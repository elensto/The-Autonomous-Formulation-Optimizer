"""Conservative preliminary gates; no automatic promotion to BO observations."""
from pathlib import Path
import csv
import json
import math
import random
import statistics as st
import sys

def rows(path):
    return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]

def mean(x):return st.mean(x)

def linefit(x,y):
    xm,ym=mean(x),mean(y)
    b=sum((a-xm)*(v-ym) for a,v in zip(x,y))/sum((a-xm)**2 for a in x)
    return ym-b*xm,b

def fit_break(x,y):
    best=None
    for k in range(4,len(x)-3):
        a,b=linefit(x[:k],y[:k]); c,d=linefit(x[k:],y[k:])
        if b<=0 or d<=1.05*b:continue
        tg=(a-c)/(d-b)
        if not x[k-1]<=tg<=x[k]:continue
        error=sum((v-(a+b*t if j<k else c+d*t))**2 for j,(t,v) in enumerate(zip(x,y)))
        if best is None or error<best['sse']:
            best={'Tg_K':tg,'Tg_C':tg-273.15,'sse':error,'low_intercept':a,'low_slope':b,'high_intercept':c,'high_slope':d,'split':k}
    return best

def equilibrium(root):
    data=rows(root/'equilibrium.dat')
    finite=all(math.isfinite(v) for row in data for v in row)
    blocksize=1000 # sampling each 100 fs -> 100 ps per block
    blocks=[data[i:i+blocksize] for i in range(0,len(data),blocksize) if len(data[i:i+blocksize])==blocksize]
    last=blocks[-10:]
    result={'finite':finite,'sampling_interval_fs':100,'block_ps':100,'complete_blocks':len(blocks),
            'stationarity_gate_passed':False,'bo_eligible':False}
    if len(last)==10 and finite:
        rho=[mean([r[1] for r in b]) for b in last]
        energy=[mean([r[3] for r in b]) for b in last]
        drift=abs(mean(rho[5:])-mean(rho[:5]))/mean(rho)
        edrift=abs(mean(energy[5:])-mean(energy[:5]))
        temp=mean([r[4] for b in last for r in b]);pressure=mean([r[5] for b in last for r in b])
        qmax=max(abs(r[6]) for r in data)
        passed=(drift<0.005 and edrift<0.02 and abs(temp-503)<5 and abs(pressure-1)<200 and qmax<1e-3)
        result.update({'density_g_cm3':mean(rho),'density_block_SE_g_cm3':st.stdev(rho)/math.sqrt(10),
                       'density_relative_half_window_drift':drift,'energy_half_window_drift_kcal_mol_atom':edrift,
                       'mean_temperature_K':temp,'mean_pressure_atm':pressure,'max_abs_charge_e':qmax,
                       'stationarity_gate_passed':passed,
                       'gate_definition':'Last 1 ns: rho drift <0.5%, PE drift <0.02 kcal/mol/atom, T within 5 K, mean P within 200 atm, charge <1e-3 e.',
                       'scope':'Practical stationarity screen; not proof of convergence or force-field accuracy.'})
    (root/'equilibrium_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    return result['stationarity_gate_passed']

def tg(root):
    points=[];blockvalues=[]
    for path in sorted(root.glob('tg_*.dat'),key=lambda p:int(p.stem.split('_')[1])):
        data=rows(path)
        if not data or not all(math.isfinite(v) for r in data for v in r):continue
        blocks=[mean([r[3] for r in data[i:i+100]]) for i in range(0,len(data),100) if len(data[i:i+100])==100]
        if len(blocks)<5:continue
        points.append({'T_K':mean([r[1] for r in data]),'specific_volume_cm3_g':mean([r[3] for r in data]),
                       'density_g_cm3':mean([r[2] for r in data]),'block_SE_specific_volume':st.stdev(blocks)/math.sqrt(len(blocks))})
        blockvalues.append(blocks)
    pairs=sorted(zip(points,blockvalues),key=lambda a:a[0]['T_K'])
    points=[a for a,b in pairs];blockvalues=[b for a,b in pairs]
    result={'status':'unresolved','temperature_points':len(points),'bo_eligible':False,
            'scope':'One network, rapid stepwise cooling. Bootstrap describes within-trajectory sampling only; excludes network, cooling-rate and force-field uncertainty.'}
    if len(points)==26:
        x=[r['T_K'] for r in points];y=[r['specific_volume_cm3_g'] for r in points]
        best=fit_break(x,y)
        if best:
            a,b=linefit(x,y);baseline=sum((v-a-b*t)**2 for t,v in zip(x,y))
            improvement=1-best['sse']/baseline if baseline>0 else 0
            rng=random.Random(42);boots=[]
            for _ in range(200):
                sampled=[mean(rng.choices(vals,k=len(vals))) for vals in blockvalues]
                f=fit_break(x,sampled)
                if f:boots.append(f['Tg_C'])
            boots.sort()
            interval=[boots[int(0.025*(len(boots)-1))],boots[int(0.975*(len(boots)-1))]] if boots else None
            resolved=(improvement>=0.2 and len(boots)>=160 and interval is not None and interval[1]-interval[0]<100)
            result.update({'status':'preliminary_fit' if resolved else 'unresolved','fit_candidate':best,
                           'SSE_improvement_over_single_line':improvement,'successful_bootstraps':len(boots),
                           'within_trajectory_bootstrap_95_interval_C':interval})
    if points:
        with (root/'tg_curve.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(points[0]));writer.writeheader();writer.writerows(points)
    (root/'tg_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    return True

if __name__=='__main__':
    action=sys.argv[1];root=Path(sys.argv[2])
    ok=equilibrium(root) if action=='equilibrium' else tg(root)
    raise SystemExit(0 if ok else 1)
