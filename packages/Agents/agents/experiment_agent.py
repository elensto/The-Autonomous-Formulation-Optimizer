"""Translate BO proposal -> resolvable network -> validation -> deterministic job bundle.
Never executes LAMMPS, launches curing, or silently rewrites network files.
"""
import argparse, hashlib, json, pathlib, re, shlex
from core import BASE, SYSTEMS, Gate, load_json, read_csv, resolve_match, validate_proposal, validate_structure, validate_protocol


def make_job(proposal, manifest, protocol, project_root=BASE, out_dir='outputs/jobs'):
    pre=validate_proposal(proposal)
    if pre.decision!='ACCEPT':return {'status':'REJECT','checks':[pre.to_dict()]},None
    row,error=resolve_match(proposal,manifest)
    if error:
        # No new network is ever generated here: hand off to optional Curing Agent.
        status='MODIFY' if error in ('NO_CURED_NETWORK_FOR_FORMULATION','NO_MATCHING_VERIFIED_CONVERSION','STRUCTURE_NOT_FOUND_OR_NOT_CURED','AMBIGUOUS_REPLICA_SPECIFY_STRUCTURE_ID') else 'REJECT'
        action='CALL_CURING_AGENT_OR_IMPORT_VALIDATED_NETWORK' if proposal['system_id']=='SYS_DDS_LIVRAGHI' else 'VERIFY_ARCHIVE_AND_MAP_CONVERSION'
        return {'status':status,'reason':error,'next_action':action,'checks':[pre.to_dict()]},None
    chem=validate_structure(row,project_root)
    ff=validate_protocol(protocol,row['system_id'])
    checks=[pre.to_dict(),chem.to_dict(),ff.to_dict()]
    decision='REJECT' if 'REJECT' in (chem.decision,ff.decision) else 'MODIFY' if 'MODIFY' in (chem.decision,ff.decision) else 'ACCEPT'
    if decision!='ACCEPT':return {'status':decision,'structure_id':row['structure_id'],'checks':checks,'next_action':'REPAIR_MISSING_METADATA_FILES_OR_PROTOCOL_AND_RECHECK'},None
    rid=str(proposal.get('request_id',''))
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,60}',rid):
        return {'status':'REJECT','reason':'REQUEST_ID_UNSAFE_OR_MISSING'},None
    base=pathlib.Path(project_root).resolve(); dest=(base/out_dir/rid).resolve()
    if not dest.is_relative_to(base):return {'status':'REJECT','reason':'OUTPUT_PATH_ESCAPES_PROJECT_ROOT'},None
    if dest.exists():return {'status':'MODIFY','reason':'JOB_DIRECTORY_ALREADY_EXISTS_NO_OVERWRITE'},None
    # templates owned by sim expert; no guessing AMBER/CHARMM styles in agent.
    tmpl=(base/protocol['template_file']).resolve()
    if not tmpl.is_relative_to(base) or not tmpl.is_file():return {'status':'MODIFY','reason':'REVIEWED_TEMPLATE_FILE_MISSING'},None
    header=tmpl.read_text()
    if '{{READ_RESTART}}' not in header:return {'status':'MODIFY','reason':'TEMPLATE_MUST_CONTAIN_READ_RESTART_PLACEHOLDER'},None
    if 'REVIEW_REQUIRED' in protocol['template_file'] or 'REPLACE THIS FILE' in header:
        return {'status':'MODIFY','reason':'FORCE_FIELD_TEMPLATE_REVIEW_REQUIRED'},None
    if 'read_data ' in header.lower():return {'status':'REJECT','reason':'TEMPLATE_MUST_USE_REVIEWED_RESTART_NOT_READ_DATA'},None
    restart=(base/row['equilibrated_restart_path']).resolve()
    script=header.replace('{{READ_RESTART}}','read_restart '+str(restart))
    # Avoid velocity reseeding each plateau, preserving cooling trajectory continuity.
    stride=int(protocol.get('averaging_stride',100))
    if stride<1:return {'status':'REJECT','reason':'INVALID_AVERAGING_STRIDE'},None
    pressure=protocol['pressure_atm'];tdamp=protocol['thermostat_damp_fs'];pdamp=protocol['barostat_damp_fs']
    temps=protocol['temperatures_K'];e=int(protocol['equil_steps_per_T']);s=int(protocol['production_steps_per_T'])
    script+='\n# Cooling steps: preliminary and cooling-rate-dependent apparent Tg\n'
    for i,T in enumerate(temps):
        script+=f'''\n# Plateau {i+1}: {T} K
fix f_npt all npt temp {T} {T} {tdamp} iso {pressure} {pressure} {pdamp}
run {e}
variable currT equal temp
variable currV equal vol
fix f_avg all ave/time {stride} 1 {stride} v_currT v_currV file temp_{i:02d}_{T:.0f}K.avg
run {s}
unfix f_avg
unfix f_npt
'''
    # Require simulation settings are reviewed separately; no silent choices of FF parameterization.
    script+='write_restart final.restart\n'
    job={'request_id':rid,'mode':proposal['mode'],'system_id':row['system_id'],'formulation_id':row.get('formulation_id'),'structure_id':row['structure_id'],'replica_id':row['replica_id'],'verified_conversion_pct':float(row['verified_cure_percent']),'target_tg_C':proposal.get('target_tg_C'),'target_range_C':proposal.get('target_range_C'),'force_field_id':row['force_field_id'],'protocol_id':protocol['protocol_id'],'simulator':'LAMMPS_classical','status':'PREPARED_NOT_EXECUTED','input_structure_path':row['relative_path'],'restart_path':row['equilibrated_restart_path'],'output_dir':str(dest.relative_to(base)),'expected_data_pattern':'temp_*.avg','requires_post_physics_check':True,'source_type':'planned_simulation_job','checks':checks}
    return {'status':'ACCEPT','structure_id':row['structure_id'],'checks':checks},(job,script,dest)

def write_job(bundle):
    job,script,dest=bundle
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'in.tg_scan.lmp').write_text(script)
    (dest/'job.json').write_text(json.dumps(job,indent=2))
    (dest/'run_local.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\ncd "$(dirname "$0")"\n# Requires a verified and locally installed LAMMPS binary\nlmp -in in.tg_scan.lmp\n')
    (dest/'run_local.sh').chmod(0o755)
    return str(dest)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--proposal',required=True);p.add_argument('--manifest',default='data/network_manifest.csv');p.add_argument('--protocol',required=True);p.add_argument('--root',default=str(BASE));p.add_argument('--out',default='outputs/experiment_decision.json');a=p.parse_args()
    base=pathlib.Path(a.root);prop=load_json(a.proposal);manifest=read_csv(a.manifest);proto=load_json(a.protocol)
    report,bundle=make_job(prop,manifest,proto,base)
    if bundle:report['job_dir']=write_job(bundle)
    pathlib.Path(a.out).parent.mkdir(parents=True,exist_ok=True);pathlib.Path(a.out).write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
