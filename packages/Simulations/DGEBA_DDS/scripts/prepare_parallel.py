"""Build an isolated 23-composition batch by removing intact source molecules.

The source box and all retained coordinates/charges/FF parameters are preserved.
Atom/topology IDs and molecule IDs are reassigned. NPT relaxes the lower density.
"""
from pathlib import Path
import collections
import argparse
import hashlib
import json
import math
import re
import shutil

ROOT=Path(__file__).resolve().parents[1]
BATCH=ROOT/'parallel_20261021'
NAMES=['Masses','Pair Coeffs','Bond Coeffs','Angle Coeffs','Dihedral Coeffs',
       'Atoms','Velocities','Bonds','Angles','Dihedrals','Impropers']

def parse(path):
    header=[];sections={};name=None
    for line in path.read_text().splitlines():
        s=line.split('#')[0].strip()
        if s in NAMES:name=s;sections[name]=[]
        elif name:
            if s:sections[name].append(s.split())
        else:header.append(line)
    return header,sections

def components(sections):
    atoms={int(a[0]):a for a in sections['Atoms']}
    graph={i:set() for i in atoms}
    for b in sections['Bonds']:
        i,j=map(int,b[2:]);graph[i].add(j);graph[j].add(i)
    todo=set(atoms);parts=[]
    while todo:
        start=min(todo);todo.remove(start);stack=[start];part={start}
        while stack:
            i=stack.pop();fresh=graph[i]&todo;todo-=fresh;part|=fresh;stack.extend(fresh)
        parts.append(sorted(part))
    return atoms,parts

def make_case(header,sections,nd,nh,index,restore=False):
    atoms,parts=components(sections)
    resin=[p for p in parts if len(p)==49];hardener=[p for p in parts if len(p)==29]
    assert len(resin)==100 and len(hardener)==50 and len(parts)==150
    selected=resin[:nd]+hardener[:nh]
    chosen=set(i for p in selected for i in p)
    remap={old:i+1 for i,old in enumerate(sorted(chosen))}
    mol={old:j+1 for j,p in enumerate(selected) for old in p}
    new={}
    for name,rows in sections.items():
        if name=='Atoms':
            new[name]=[[str(remap[int(a[0])]),str(mol[int(a[0])])]+a[2:] for a in rows if int(a[0]) in chosen]
        elif name=='Velocities':
            new[name]=[[str(remap[int(a[0])])]+a[1:] for a in rows if int(a[0]) in chosen]
        elif name in ['Bonds','Angles','Dihedrals','Impropers']:
            result=[]
            for row in rows:
                ids=list(map(int,row[2:]));flags=[i in chosen for i in ids]
                assert all(flags) or not any(flags),'Molecule selection cuts a bonded interaction'
                if all(flags):result.append([str(len(result)+1),row[1]]+[str(remap[i]) for i in ids])
            new[name]=result
        else:new[name]=rows
    count_names={'atoms':'Atoms','bonds':'Bonds','angles':'Angles','dihedrals':'Dihedrals','impropers':'Impropers'}
    new_header=[]
    for i,line in enumerate(header):
        m=re.fullmatch(r'\s*\d+\s+(atoms|bonds|angles|dihedrals|impropers)\s*',line)
        if m:line=f'{len(new.get(count_names[m[1]],[]))} {m[1]}'
        if i==0:line=f'FormulationLab DGEBA/DDS: {nd} DGEBA + {nh} DDS; whole-molecule subset of Livraghi replica 1'
        new_header.append(line)
    case_id=f'D{nd:03d}_H{nh:03d}'
    folder=BATCH/'cases'/case_id;folder.mkdir(parents=True,exist_ok=restore)
    text='\n'.join(new_header).rstrip()+'\n\n'
    for name,rows in new.items():
        title='Atoms # full' if name=='Atoms' else name
        text+=title+'\n\n'+'\n'.join(' '.join(row) for row in rows)+'\n\n'
    (folder/'data.initial').write_text(text,encoding='utf-8',newline='\n')
    masses={int(r[0]):float(r[1]) for r in new['Masses']}
    lengths=[]
    for axis in 'xyz':
        m=re.search(rf'([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+{axis}lo\s+{axis}hi',text)
        lengths.append(float(m[2])-float(m[1]))
    charge=sum(float(a[3]) for a in new['Atoms'])
    density=sum(masses[int(a[2])] for a in new['Atoms'])/6.02214076e23/(math.prod(lengths)*1e-24)
    check_atoms,check_parts=components(new)
    histogram=collections.Counter(map(len,check_parts))
    assert histogram=={49:nd,29:nh},histogram
    assert len(check_atoms)==49*nd+29*nh and abs(charge)<1e-3
    for old,newid in remap.items():
        assert check_atoms[newid][2:]==atoms[old][2:],'Retained physical atom properties changed'
    bond_lengths=[]
    for row in new['Bonds']:
        a,b=(check_atoms[int(i)] for i in row[2:]);dr=[]
        for k,L in enumerate(lengths):
            d=float(a[k+4])-float(b[k+4]);d-=round(d/L)*L;dr.append(d)
        bond_lengths.append(math.sqrt(sum(v*v for v in dr)))
    assert min(bond_lengths)>0.5 and max(bond_lengths)<2.5
    epoxy=2*nd;amine=4*nh;limit=min(epoxy,amine)
    cfg={'array_index':index,'formulation_id':case_id,'DGEBA_molecules':nd,'DDS_molecules':nh,
         'atoms':len(check_atoms),'epoxy_groups':epoxy,'amine_hydrogens':amine,
         'amine_epoxy_ratio':amine/epoxy,'maximum_epoxy_conversion_percent':100*limit/epoxy,
         'target_limiting_reagent_conversion_percent':90,
         'required_reactions':math.ceil(0.90*limit-1e-10),
         'target_epoxy_conversion_percent':90*limit/epoxy,
         'source_replica':1,'seed':4928459,'initial_net_charge_e':charge,
         'initial_density_g_cm3':density,'minimum_image_bond_length_A':[min(bond_lengths),max(bond_lengths)],
         'data_sha256':hashlib.sha256((folder/'data.initial').read_bytes()).hexdigest(),
         'retained_physical_atom_properties_unchanged':True,'whole_molecule_selection_verified':True,
         'bo_eligible':False,'note':'Single realization per composition. Conversion differs with stoichiometry; retain it alongside Tg.'}
    (folder/'case.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    return cfg

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--restore-missing-data',action='store_true');args=parser.parse_args()
    source=ROOT/'raw/lammps/data/small/data.100dgeba0_50dds_1gcm3_1'
    if BATCH.exists():
        if not args.restore_missing_data:raise SystemExit('Batch already exists; inspect it before rebuilding')
        manifest=json.loads((BATCH/'manifest.json').read_text())
        assert hashlib.sha256(source.read_bytes()).hexdigest()==manifest['source_sha256']
        header,sections=parse(source);restored=0
        for cfg in manifest['cases']:
            if not (BATCH/'cases'/cfg['formulation_id']/'data.initial').exists():
                new=make_case(header,sections,cfg['DGEBA_molecules'],cfg['DDS_molecules'],cfg['array_index'],restore=True)
                assert new['data_sha256']==cfg['data_sha256'];restored+=1
        print('Restored missing data files:',restored);return
    common=BATCH/'common';common.mkdir(parents=True)
    header,sections=parse(source)
    recipes=[(100,h) for h in range(40,50)]+[(d,50) for d in range(99,86,-1)]
    cases=[make_case(header,sections,d,h,i) for i,(d,h) in enumerate(recipes)]
    for name in ['coefficients_reaction.inc','preTop.txt','postTop.txt','2preTop.txt','2postTop.txt','map1.txt','map2.txt','in.equilibrate','in.tg_screen']:
        shutil.copy2(ROOT/'inputs'/name,common/name)
    pilot=(ROOT/'inputs/in.pilot').read_text().replace('variable heat_steps index 2000','variable epoxies index 200\nvariable heat_steps index 2000').replace('100.0*c_noh/200.0','100.0*c_noh/v_epoxies')
    (common/'in.prepare').write_text(pilot,encoding='utf-8',newline='\n')
    cure=(ROOT/'inputs/in.cure_campaign').read_text().replace('variable chunk_steps index 100000','variable epoxies index 200\nvariable chunk_steps index 100000').replace('100.0*c_noh/200.0','100.0*c_noh/v_epoxies')
    (common/'in.cure_campaign').write_text(cure,encoding='utf-8',newline='\n')
    for name in ['export_topology.py','analyze_properties.py']:
        shutil.copy2(ROOT/'scripts'/name,common/name)
    shutil.copy2(ROOT/'scripts/analyze_campaign.py',common/'analyze_reference.py')
    for name in ['parallel_analysis.py','parallel_runner.sh','parallel_array.slurm','parallel_validate.slurm']:
        shutil.copy2(ROOT/'scripts'/name,common/name)
    manifest={'batch_id':BATCH.name,'source_file':str(source.relative_to(ROOT)),
              'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'additional_cases':len(cases),'cpus_per_case':16,'maximum_array_parallelism':23,
              'source_geometry':'No insertion or coordinate scaling: remove intact molecules, renumber IDs, relax by NPT.',
              'cure_rule':'90% of the limiting reagent capacity; actual epoxy conversion retained in every result.',
              'initial_reference_job':'175498','cases':cases,
              'common_file_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in common.iterdir() if p.is_file()}}
    (BATCH/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    (BATCH/'case_ids.txt').write_text('\n'.join(c['formulation_id'] for c in cases)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'batch':str(BATCH),'cases':len(cases),'atom_range':[min(c['atoms'] for c in cases),max(c['atoms'] for c in cases)],'ratio_range':[min(c['amine_epoxy_ratio'] for c in cases),max(c['amine_epoxy_ratio'] for c in cases)]},indent=2))

if __name__=='__main__':main()
