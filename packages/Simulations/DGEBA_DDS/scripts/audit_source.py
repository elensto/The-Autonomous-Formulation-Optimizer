"""Structural checks on the uncured small source, independent of MD execution."""
from pathlib import Path
import collections
import hashlib
import json
import math
import re

ROOT=Path(__file__).resolve().parents[1]
def main():
    path=ROOT/'inputs/data.initial';text=path.read_text(); sections={};section=None;box={}
    for line in text.splitlines():
        clean=line.split('#')[0].strip()
        if not clean:continue
        if re.match(r'^[A-Za-z]',clean):
            section=clean;sections[section]=[];continue
        fields=clean.split()
        if len(fields)==4 and fields[2] in ('xlo','ylo','zlo'):
            box[fields[2][0]]=float(fields[1])-float(fields[0])
        if section:sections[section].append(fields)
    atoms={int(a[0]):a for a in sections['Atoms']}
    counts=collections.Counter(int(a[2]) for a in atoms.values())
    masses={int(row[0]):float(row[1]) for row in sections['Masses']}
    distances=[]
    for row in sections['Bonds']:
        a,b=atoms[int(row[2])],atoms[int(row[3])]; d=[]
        for axis,idx in zip('xyz',(4,5,6)):
            value=float(a[idx])-float(b[idx]);value-=round(value/box[axis])*box[axis];d.append(value)
        distances.append(math.sqrt(sum(v*v for v in d)))
    for name,n in [('Bonds',2),('Angles',3),('Dihedrals',4)]:
        for row in sections[name]:assert all(int(a) in atoms for a in row[-n:])
    net=sum(float(a[3]) for a in atoms.values())
    # Original files use molecule ID zero; identify molecules by bond connectivity.
    graph={a:set() for a in atoms}
    for row in sections['Bonds']:
        a,b=int(row[2]),int(row[3]);graph[a].add(b);graph[b].add(a)
    unseen=set(atoms);sizes=collections.Counter()
    while unseen:
        start=unseen.pop();todo=[start];size=0
        while todo:
            node=todo.pop();size+=1
            fresh=graph[node]&unseen;unseen-=fresh;todo.extend(fresh)
        sizes[size]+=1
    result={'atoms':len(atoms),'bonds':len(sections['Bonds']),'angles':len(sections['Angles']),
            'dihedrals':len(sections['Dihedrals']),'molecule_atom_count_histogram':dict(sizes),
            'net_charge_e':net,'initial_density_g_cm3':sum(masses[int(a[2])] for a in atoms.values())/6.02214076e23/(math.prod(box.values())*1e-24),
            'minimum_image_bond_length_A':[min(distances),max(distances)],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'checks_passed':len(atoms)==6350 and abs(net)<1e-3 and min(distances)>0.5 and max(distances)<2.5,
            'scope':'Connectivity/reference checks only; force-field accuracy and equilibrium not established.'}
    assert sizes=={49:100,29:50},sizes
    (ROOT/'reports/source_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    assert result['checks_passed']

if __name__=='__main__':main()
