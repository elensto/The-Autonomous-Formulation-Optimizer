"""DDS curing adapter: validate/prepare external, scientist-reviewed LAMMPS workflows.

No chemistry is invented here; the actual Livraghi-derived bond/react deck must be
provided by the MD team. No scheduler submission is performed.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, pathlib, re, shlex
from core import BASE, SYSTEMS, Gate, load_json, read_csv, sha256, validate_proposal

SAFE_ID=re.compile(r'^[A-Za-z0-9_-]{1,60}$')

def validate_curing_request(proposal, config, root=BASE):
    errors=[]; warnings=[]
    if not isinstance(proposal,dict):return {'decision':'REJECT','reasons':['INVALID_PROPOSAL'],'warnings':[]}
    if proposal.get('system_id')!='SYS_DDS_LIVRAGHI':errors.append('CURING_BACKEND_ONLY_SUPPORTS_DGEBA_DDS')
    if proposal.get('mode') not in ('cure_conversion','formulation'):errors.append('INVALID_MODE')
    rid=proposal.get('request_id','')
    if not isinstance(rid,str) or not SAFE_ID.fullmatch(rid):errors.append('INVALID_REQUEST_ID')
    target=proposal.get('conversion_pct')
    try:
        target=float(target)
        if not 0<target<=100:errors.append('INVALID_TARGET_CONVERSION')
    except (ValueError,TypeError):errors.append('TARGET_CONVERSION_REQUIRED');target=None
    if config.get('backend')!='external_lammps_bond_react':errors.append('UNSUPPORTED_CURING_BACKEND')
    if config.get('system_id')!='SYS_DDS_LIVRAGHI':errors.append('CONFIG_CHEMISTRY_MISMATCH')
    if config.get('force_field_family')!='AMBER':errors.append('DDS_MUST_USE_REVIEWED_AMBER_FORCE_FIELD')
    if not config.get('scientist_reviewed'):errors.append('CURING_WORKFLOW_NOT_SCIENTIST_REVIEWED')
    if not config.get('input_script_sha256'):errors.append('CURING_SCRIPT_CHECKSUM_REQUIRED')
    if not config.get('starting_structure_sha256'):errors.append('START_STRUCTURE_CHECKSUM_REQUIRED')
    root=pathlib.Path(root).resolve()
    for key,hkey in [('input_script','input_script_sha256'),('starting_structure','starting_structure_sha256')]:
        rel=config.get(key)
        if not isinstance(rel,str) or not rel:errors.append(key.upper()+'_NOT_CONFIGURED');continue
        p=(root/rel).resolve()
        if not p.is_relative_to(root):errors.append(key.upper()+'_PATH_ESCAPE');continue
        if not p.is_file():errors.append(key.upper()+'_NOT_FOUND');continue
        if config.get(hkey) and sha256(p)!=config[hkey].lower():errors.append(key.upper()+'_HASH_MISMATCH')
    # Proposal count fields are required if changing formulation. Never silently
    # assume the published 100:50 baseline represents a proposed other ratio.
    counts=(proposal.get('components') or {}).get('molecule_counts')
    if proposal.get('mode')=='formulation':
        comps=proposal.get('components') or {}
        if comps.get('resin_id')!='R001' or comps.get('hardener_id')!='H007':errors.append('FORMULATION_CHEMISTRY_UNSUPPORTED')
        if not isinstance(counts,dict):errors.append('MOLECULE_COUNTS_REQUIRED_FOR_NEW_FORMULATION')
    if not isinstance(counts,dict):
        if proposal.get('mode')=='cure_conversion':counts=config.get('starting_counts')
    if not isinstance(counts,dict):errors.append('STARTING_MOLECULE_COUNTS_UNAVAILABLE')
    else:
        a,b=counts.get('DGEBA'),counts.get('DDS')
        if type(a)!=int or type(b)!=int or a<=0 or b<=0:errors.append('INVALID_MOLECULE_COUNTS')
        else:
            ratio=4*b/(2*a)
            if ratio<=0:errors.append('INVALID_EQUIVALENCE_RATIO')
            # Chemistry changes must have independently built input and deck.
            if counts!=config.get('starting_counts'):
                errors.append('NEW_COMPOSITION_REQUIRES_REBUILT_START_AND_REVIEWED_REACTION_DECK')
            if target is not None and target>100*min(1,ratio):
                warnings.append('TARGET_EXCEEDS_EPOXY_SIDE_STOICHIOMETRIC_CEILING_IF_CONVERSION_NORMALIZED_TO_EPOXIDES')
            warnings.append(f'AMINE_H_TO_EPOXY_RATIO={ratio:.4f}')
    decision='MODIFY' if errors else 'ACCEPT'
    return {'decision':decision,'reasons':errors,'warnings':warnings,'requested_conversion_pct':target,'starting_counts':counts,'source':'validation_only'}

def prepare_curing_job(proposal,config,root=BASE,out_dir='outputs/curing_jobs'):
    root=pathlib.Path(root).resolve()
    gate=validate_curing_request(proposal,config,root)
    if gate['decision']!='ACCEPT':return gate,None
    dest=(root/out_dir/proposal['request_id']).resolve()
    if not dest.is_relative_to(root):return {'decision':'REJECT','reasons':['DESTINATION_PATH_ESCAPE']},None
    if dest.exists():return {'decision':'MODIFY','reasons':['JOB_DIRECTORY_EXISTS']},None
    script=(root/config['input_script']).resolve()
    start=(root/config['starting_structure']).resolve()
    record={
        'request_id':proposal['request_id'],'system_id':'SYS_DDS_LIVRAGHI',
        'mode':proposal['mode'],'formulation_id':proposal.get('formulation_id','DGEBA_DDS_EQUIV1'),
        'requested_conversion_pct':float(proposal['conversion_pct']),
        'achieved_conversion_pct':None,'source_type':'curing_job_request',
        'starting_counts':gate['starting_counts'],
        'input_script':str(script),'input_script_sha256':sha256(script),
        'starting_structure':str(start),'starting_structure_sha256':sha256(start),
        'force_field_family':'AMBER','cure_protocol_id':config['cure_protocol_id'],
        'status':'PREPARED_NOT_SUBMITTED',
        'notes':'Scientific input script supplied externally; output must be independently reviewed before import to cured manifest.'
    }
    return gate,(dest,record)

def write_curing_job(bundle):
    dest,record=bundle
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'curing_job.json').write_text(json.dumps(record,indent=2)+'\n')
    # Preparation never edits the bond/react script or submits scheduler jobs.
    (dest/'README.txt').write_text('Review input_script and starting_structure in curing_job.json.\nRun using the HPC team approved procedure. No command is submitted automatically.\nAfter completion, record actual conversion and validated topology before importing the network.\n')
    return dest

def assess_curing_campaign(event):
    """Interpret external campaign reports; not a live HPC query."""
    state=str(event.get('status','')).upper()
    achieved=event.get('achieved_conversion_pct')
    target=event.get('target_conversion_pct')
    reasons=[]
    if state not in ('COMPLETED','SUCCESS'):reasons.append('CURING_CAMPAIGN_NOT_COMPLETED')
    if achieved is None:reasons.append('ACHIEVED_CONVERSION_NOT_MEASURED')
    else:
        try:
            if not 0<=float(achieved)<=100:reasons.append('INVALID_ACHIEVED_CONVERSION')
        except (ValueError,TypeError):reasons.append('INVALID_ACHIEVED_CONVERSION')
    for k in ('validated_bond_graph','charge_conservation_checked','force_field_reviewed','topology_reviewed'):
        if event.get(k) is not True:reasons.append(k.upper()+'_PENDING')
    if not event.get('structure_file_sha256'):reasons.append('STRUCTURE_HASH_MISSING')
    return {'decision':'MODIFY' if reasons else 'ACCEPT','reasons':reasons,'target_conversion_pct':target,'achieved_conversion_pct':achieved,'meaning':'CURING_OUTPUT_IMPORT_GATE_ONLY_NOT_TG_MD_APPROVAL'}

def main():
    p=argparse.ArgumentParser(description='Prepare or assess DGEBA/DDS external LAMMPS curing')
    sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('prepare');a.add_argument('--proposal',required=True);a.add_argument('--config',default='config/curing_DDS_EXTERNAL.json');a.add_argument('--root',default=str(BASE))
    b=sub.add_parser('assess');b.add_argument('--event',required=True)
    args=p.parse_args()
    if args.action=='assess':result=assess_curing_campaign(load_json(args.event))
    else:
        base=pathlib.Path(args.root).resolve();result,bundle=prepare_curing_job(load_json(args.proposal),load_json(base/args.config),base)
        if bundle:result['job_dir']=str(write_curing_job(bundle))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
