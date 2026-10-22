"""One-command dual-path dispatcher. Never submits jobs or imports network data."""
import argparse, json, pathlib
from core import BASE, load_json, read_csv
from experiment_agent import make_job, write_job
from curing_agent import prepare_curing_job, write_curing_job

def route(proposal,manifest,protocol,curing_config,root=BASE):
    gate,bundle=make_job(proposal,manifest,protocol,root)
    if bundle:
        gate['job_dir']=write_job(bundle)
        gate['route']='EXISTING_CURED_NETWORK_TO_TG'
        return gate
    if gate.get('next_action')=='CALL_CURING_AGENT_OR_IMPORT_VALIDATED_NETWORK' and proposal.get('system_id')=='SYS_DDS_LIVRAGHI':
        check,cure_bundle=prepare_curing_job(proposal,curing_config,root)
        if cure_bundle:check['curing_job_dir']=str(write_curing_job(cure_bundle))
        return {'status':check['decision'],'route':'DDS_CURING_REQUIRED','experiment_gate':gate,'curing_gate':check,'note':'Never substitute requested 90% for achieved 90%.'}
    gate['route']='PREBUILT_NETWORK_LOOKUP';return gate

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--proposal',required=True)
    p.add_argument('--manifest',default='data/network_manifest.csv')
    p.add_argument('--dds-protocol',default='config/protocol_DDS.json')
    p.add_argument('--d230-protocol',default='config/protocol_D230.json')
    p.add_argument('--curing-config',default='config/curing_DDS_EXTERNAL.json')
    p.add_argument('--root',default=str(BASE))
    args=p.parse_args();root=pathlib.Path(args.root).resolve()
    prop=load_json(root/args.proposal)
    proto=load_json(root/(args.dds_protocol if prop.get('system_id')=='SYS_DDS_LIVRAGHI' else args.d230_protocol))
    report=route(prop,read_csv(root/args.manifest),proto,load_json(root/args.curing_config),root)
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
