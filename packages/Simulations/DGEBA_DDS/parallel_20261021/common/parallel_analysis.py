"""Composition-aware structural gates. Never promotes a result to BO automatically."""
from pathlib import Path
import json
import math
import sys
from analyze_reference import analyze

def main():
    root=Path(sys.argv[1]);cfg=json.loads((root/'case.json').read_text())
    result=analyze(root/'final.cure.data')
    result['epoxy_conversion_percent']=100*result['reacted_epoxies']/cfg['epoxy_groups']
    result['limiting_reagent_conversion_percent']=100*result['reacted_epoxies']/min(cfg['epoxy_groups'],cfg['amine_hydrogens'])
    result['target_reached']=(root/'target_reached.txt').read_text().strip()=='true'
    result['formulation_id']=cfg['formulation_id']
    result['structure_gate_passed']=(result['target_reached'] and result['atoms']==cfg['atoms']
        and math.isfinite(result['net_charge_e']) and abs(result['net_charge_e'])<1e-3
        and cfg['required_reactions']<=result['reacted_epoxies']<=min(cfg['epoxy_groups'],cfg['amine_hydrogens'])
        and result['largest_component_fraction']>=0.90)
    result.update({'DGEBA_molecules':cfg['DGEBA_molecules'],'DDS_molecules':cfg['DDS_molecules'],
                   'amine_epoxy_ratio':cfg['amine_epoxy_ratio'],
                   'target_epoxy_conversion_percent':cfg['target_epoxy_conversion_percent'],
                   'equilibrated':False,'tg_calculated':False,'bo_eligible':False,
                   'note':'Connectivity is a practical screen, not a finite-size/percolation validation. Cure conversion varies with composition.'})
    (root/'campaign_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
