"""Standalone before/after PhysCheck decision CLI."""
import argparse,json,pathlib
from core import BASE, load_json, read_csv, resolve_match, validate_proposal, validate_structure, validate_protocol,validate_result

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='stage',required=True)
    a=s.add_parser('pre');a.add_argument('--proposal',required=True);a.add_argument('--manifest',default='data/network_manifest.csv');a.add_argument('--protocol',required=True);a.add_argument('--root',default=str(BASE))
    b=s.add_parser('post');b.add_argument('--result',required=True)
    a=p.parse_args()
    if a.stage=='post':report=validate_result(load_json(a.result)).to_dict()
    else:
        prop=load_json(a.proposal);rows=read_csv(a.manifest);first=validate_proposal(prop)
        if first.decision=='REJECT':report=first.to_dict()
        else:
            row,error=resolve_match(prop,rows)
            if error:report={'stage':'pre','decision':'MODIFY','reasons':[error],'warnings':[]}
            else:
                ss=validate_structure(row,pathlib.Path(a.root));pp=validate_protocol(load_json(a.protocol),row['system_id'])
                all_reasons=ss.reasons+pp.reasons
                report={'stage':'pre','decision':'REJECT' if 'REJECT' in (ss.decision,pp.decision) else 'MODIFY' if all_reasons else 'ACCEPT','reasons':all_reasons,'warnings':ss.warnings+pp.warnings,'structure_id':row['structure_id']}
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
