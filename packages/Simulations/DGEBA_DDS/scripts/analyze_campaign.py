"""Check chemistry/conversion/connectivity before any property calculation."""
from pathlib import Path
import collections
import json
import math
import sys

def analyze(path):
    sections={};key=None
    for line in path.read_text().splitlines():
        s=line.split('#')[0].strip()
        if not s:continue
        if s[0].isalpha():key=s;sections[key]=[]
        elif key:sections[key].append(s.split())
    atoms={int(a[0]):a for a in sections['Atoms']}
    graph={i:set() for i in atoms}
    for row in sections['Bonds']:
        a,b=int(row[2]),int(row[3]);graph[a].add(b);graph[b].add(a)
    todo=set(atoms);sizes=[]
    while todo:
        stack=[todo.pop()];size=0
        while stack:
            n=stack.pop();size+=1;fresh=graph[n]&todo;todo-=fresh;stack.extend(fresh)
        sizes.append(size)
    types=collections.Counter(int(a[2]) for a in atoms.values())
    charge=sum(float(a[3]) for a in atoms.values())
    return {'atoms':len(atoms),'net_charge_e':charge,'reacted_epoxies':types[12],
            'epoxy_conversion_percent':100*types[12]/200,
            'connected_components':len(sizes),'largest_component_atoms':max(sizes),
            'largest_component_fraction':max(sizes)/len(atoms)}

def main():
    root=Path(sys.argv[1]); result=analyze(root/'final.cure.data')
    result['target_reached']=(root/'target_reached.txt').read_text().strip()=='true'
    result['structure_gate_passed']=(result['target_reached'] and result['atoms']==6350
        and math.isfinite(result['net_charge_e']) and abs(result['net_charge_e'])<1e-3
        and result['epoxy_conversion_percent']>=90 and result['largest_component_fraction']>=0.90)
    result.update({'equilibrated':False,'tg_calculated':False,'bo_eligible':False,
                   'note':'Connectivity gate is a practical screen, not a finite-size/percolation validation.'})
    (root/'campaign_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
