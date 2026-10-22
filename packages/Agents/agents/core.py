"""Zero-dependency interfaces and conservative physics gates for cured epoxy Tg."""
from __future__ import annotations
import csv, hashlib, json, math, pathlib, re
from dataclasses import dataclass, asdict
from typing import Any

BASE = pathlib.Path(__file__).resolve().parents[1]
TRUE = {'true', 'yes', '1'}
SYSTEMS = {
    'SYS_D230_FAN': {'resin_id':'R001','hardener_id':'H006','label':'DGEBA/Jeffamine D230','file_mode':'prebuilt_network','expected_force_field':'FAN_CHARMM_REVIEWED'},
    'SYS_DDS_LIVRAGHI': {'resin_id':'R001','hardener_id':'H007','label':'DGEBA/DDS','file_mode':'cure_or_prebuilt_network','expected_force_field':'DDS_AMBER_REVIEWED'},
}
REQUIRED_TOPOLOGY = ('atoms','bonds','angles','dihedrals','impropers')

def load_json(path):
    return json.loads(pathlib.Path(path).read_text())

def read_csv(path):
    with open(path, newline='') as f: return list(csv.DictReader(f))

def is_true(v): return str(v).lower() in TRUE

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def parse_lammps_header(path):
    """Lightweight checks only: counts and presence of known sections (NOT force-field validation)."""
    txt=pathlib.Path(path).read_text(errors='replace')
    counts={}
    for name in REQUIRED_TOPOLOGY:
        m=re.search(r'^\s*(\d+)\s+'+name+r'\s*$',txt,re.M)
        counts[name]=int(m.group(1)) if m else None
    section_names=['Masses','Atoms','Bonds','Angles','Dihedrals','Impropers','Pair Coeffs','Bond Coeffs','Angle Coeffs','Dihedral Coeffs','Improper Coeffs']
    present={s:bool(re.search(r'^'+re.escape(s)+r'(?:\s*#.*)?\s*$',txt,re.M)) for s in section_names}
    # Explicitly reject missing bonds or atoms on both branches.
    return {'counts':counts,'sections':present,'nonempty':all(counts[k] is not None and counts[k]>0 for k in ('atoms','bonds')) and present['Masses'] and present['Atoms'] and present['Bonds']}

@dataclass
class Gate:
    decision: str
    reasons: list[str]
    warnings: list[str]
    stage: str
    def to_dict(self):return asdict(self)


def resolve_match(proposal,manifest):
    """Only exact system IDs and existing IDs; no fabricated conversion levels."""
    mode=proposal.get('mode')
    if mode not in ('cure_conversion','formulation'):
        return None,'INVALID_MODE'
    sid=proposal.get('system_id')
    if sid not in SYSTEMS: return None,'UNSUPPORTED_SYSTEM'
    candidates=[r for r in manifest if r.get('system_id')==sid and r.get('record_kind') in ('cured_network','generated_cured_network')]
    requested_id=proposal.get('structure_id')
    if requested_id:
        found=[r for r in candidates if r['structure_id']==requested_id]
        return (found[0],None) if len(found)==1 else (None,'STRUCTURE_NOT_FOUND_OR_NOT_CURED')
    # Formulation proposals must be matched to a discrete prebuilt network with identical formulation identifier.
    fid=proposal.get('formulation_id')
    if mode=='formulation':
        if not fid: return None,'FORMULATION_ID_REQUIRED'
        candidates=[r for r in candidates if r.get('formulation_id')==fid]
        if not candidates:return None,'NO_CURED_NETWORK_FOR_FORMULATION'
    if proposal.get('conversion_pct') is None:return None,'CONVERSION_OR_STRUCTURE_ID_REQUIRED'
    try: target=float(proposal['conversion_pct'])
    except(ValueError,TypeError):return None,'INVALID_CONVERSION'
    if not 0<=target<=100:return None,'INVALID_CONVERSION'
    candidates=[r for r in candidates if r.get('verified_cure_percent','').strip() and abs(float(r['verified_cure_percent'])-target)<1e-7]
    if proposal.get('replica_id') is not None:
        candidates=[r for r in candidates if str(r.get('replica_id'))==str(proposal['replica_id'])]
    if not candidates:return None,'NO_MATCHING_VERIFIED_CONVERSION'
    if len(candidates)>1:return None,'AMBIGUOUS_REPLICA_SPECIFY_STRUCTURE_ID'
    return candidates[0],None


def validate_proposal(proposal):
    reasons=[]
    if not isinstance(proposal,dict):return Gate('REJECT',['PROPOSAL_NOT_OBJECT'],[],'proposal')
    if proposal.get('mode') not in ('cure_conversion','formulation'):reasons.append('MODE_INVALID')
    sid=proposal.get('system_id')
    if sid not in SYSTEMS:reasons.append('SYSTEM_ID_UNSUPPORTED')
    if proposal.get('mode')=='formulation':
        # Prevent silent use of a different hardener, resin or generated chemistry.
        if not proposal.get('formulation_id'):reasons.append('FORMULATION_ID_REQUIRED')
        comps=proposal.get('components')
        if not isinstance(comps,dict):reasons.append('COMPONENTS_REQUIRED_FOR_FORMULATION')
        elif sid in SYSTEMS and (comps.get('resin_id')!=SYSTEMS[sid]['resin_id'] or comps.get('hardener_id')!=SYSTEMS[sid]['hardener_id']):reasons.append('CHEMISTRY_NOT_SUPPORTED_BY_SELECTED_SYSTEM')
    v=proposal.get('conversion_pct')
    if v is not None:
        try:
            if not 0<=float(v)<=100:reasons.append('CONVERSION_OUT_OF_RANGE')
        except (ValueError,TypeError):reasons.append('CONVERSION_NOT_NUMERIC')
    if not (proposal.get('structure_id') or v is not None):reasons.append('STRUCTURE_ID_OR_CONVERSION_REQUIRED')
    return Gate('REJECT' if reasons else 'ACCEPT',reasons,[],'proposal')


def validate_structure(row, project_root=BASE):
    """Metadata flags are checked but do not themselves prove physical accuracy."""
    if row is None:return Gate('MODIFY',['NO_STRUCTURE_RESOLVED'],[],'pre')
    reasons=[];warn=[]
    if row.get('record_kind') not in ('cured_network','generated_cured_network'):
        return Gate('REJECT',['NOT_A_CURED_NETWORK'],[],'pre')
    if row.get('system_id') not in SYSTEMS:
        return Gate('REJECT',['UNKNOWN_NETWORK_CHEMISTRY'],[],'pre')
    if row.get('verified_cure_percent','').strip()=='':reasons.append('CURE_PERCENT_UNVERIFIED')
    else:
        try:
            x=float(row['verified_cure_percent'])
            if x<0 or x>100:reasons.append('CURE_PERCENT_OUT_OF_RANGE')
        except ValueError:reasons.append('CURE_PERCENT_INVALID')
    if not row.get('replica_id','').strip():reasons.append('REPLICA_ID_MISSING')
    path=row.get('relative_path','').strip()
    if not path:reasons.append('NETWORK_PATH_MISSING');file=None
    else:
        file=(pathlib.Path(project_root)/path).resolve()
        if not file.is_relative_to(pathlib.Path(project_root).resolve()):return Gate('REJECT',['STRUCTURE_PATH_ESCAPES_PROJECT_ROOT'],[],'pre')
        if not file.is_file():reasons.append('NETWORK_FILE_MISSING')
    if file and file.exists():
        expected=row.get('file_sha256','').strip().lower()
        if not expected:reasons.append('FILE_HASH_UNVERIFIED')
        elif sha256(file)!=expected:return Gate('REJECT',['FILE_HASH_MISMATCH'],[],'pre')
        hdr=parse_lammps_header(file)
        if not hdr['nonempty']:return Gate('REJECT',['LAMMPS_HEADER_MISSING_ATOMS_BONDS_OR_SECTIONS'],[],'pre')
        if row.get('number_atoms','').strip() and int(row['number_atoms'])!=hdr['counts']['atoms']:
            return Gate('REJECT',['LAMMPS_ATOM_COUNT_MISMATCH'],[],'pre')
        warn.append('HEADER_CHECK_DOES_NOT_VALIDATE_BOND_GRAPH_OR_FORCE_FIELD')
    for k in ('topology_reviewed','force_field_reviewed','npt_equilibrated'):
        if not is_true(row.get(k,'')):reasons.append(k.upper()+'_NOT_APPROVED')
    ff=row.get('force_field_id','')
    if ff!=SYSTEMS[row['system_id']]['expected_force_field']:reasons.append('FORCE_FIELD_FAMILY_MISMATCH')
    if not row.get('equilibrated_restart_path',''):
        reasons.append('EQUILIBRATED_RESTART_MISSING')
    else:
        rst=(pathlib.Path(project_root)/row['equilibrated_restart_path']).resolve()
        if not rst.is_relative_to(pathlib.Path(project_root).resolve()):return Gate('REJECT',['RESTART_PATH_ESCAPES_PROJECT_ROOT'],[],'pre')
        if not rst.is_file():reasons.append('EQUILIBRATED_RESTART_FILE_MISSING')
    if row['system_id']=='SYS_DDS_LIVRAGHI' and row.get('record_kind')=='generated_cured_network' and not row.get('cure_protocol_id',''):
        reasons.append('DDS_CURING_PROVENANCE_MISSING')
    return Gate('MODIFY' if reasons else 'ACCEPT',reasons,warn,'pre')


def validate_protocol(protocol,system_id):
    issues=[]
    if protocol.get('engine')!='lammps' or protocol.get('ensemble')!='NPT':issues.append('ENGINE_OR_ENSEMBLE_UNSUPPORTED')
    if not protocol.get('reviewed_by_scientist'):issues.append('PROTOCOL_NOT_REVIEWED')
    if protocol.get('force_field_id')!=SYSTEMS[system_id]['expected_force_field']:issues.append('FORCE_FIELD_PROTOCOL_MISMATCH')
    if protocol.get('template_file') in (None,''):issues.append('REVIEWED_LAMMPS_TEMPLATE_REQUIRED')
    for key in ('timestep_fs','equil_steps_per_T','production_steps_per_T'):
        if not isinstance(protocol.get(key), (int,float)) or protocol[key]<=0:issues.append(key.upper()+'_INVALID')
    temps=protocol.get('temperatures_K',[])
    if len(temps)<7 or any(not isinstance(t,(int,float)) for t in temps):issues.append('TEMPERATURE_GRID_TOO_SHORT')
    elif not all(temps[i]>temps[i+1] for i in range(len(temps)-1)):issues.append('TEMPERATURE_GRID_NOT_DESCENDING')
    if protocol.get('pressure_atm') != 1.0:issues.append('PRESSURE_NOT_ONE_ATM_REVIEW_REQUIRED')
    return Gate('MODIFY' if issues else 'ACCEPT',issues,[],'protocol')


def validate_result(result):
    """Never fabricate or accept Tg from reference literature or MLIP single points."""
    issues=[]
    if result.get('source_type') != 'real_classical_MD':issues.append('NOT_REAL_CLASSICAL_MD')
    if result.get('property')!='Tg':issues.append('PROPERTY_NOT_TG')
    if not result.get('fit_accepted'):issues.append('FIT_NOT_ACCEPTED')
    if not result.get('npt_converged'):issues.append('NPT_NOT_CONVERGED')
    if not result.get('protocol_reviewed'):issues.append('PROTOCOL_NOT_REVIEWED')
    if not result.get('structure_id'):issues.append('STRUCTURE_ID_MISSING')
    for k in ('tg_K','tg_uncertainty_K'):
        try:
            v=float(result[k]);assert math.isfinite(v) and v>0
        except (KeyError,TypeError,ValueError,AssertionError):issues.append(k.upper()+'_INVALID')
    if result.get('cooling_trace_path') in (None,''):issues.append('TRACE_MISSING')
    return Gate('MODIFY' if issues else 'ACCEPT',issues,[],'post')
