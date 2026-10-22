"""Summarize a bounded pilot; never promotes transient density or cure to BO data."""
from pathlib import Path
import collections
import json
import math
import re
import sys

def parse_data(path):
    counts=collections.Counter(); charge=0.0; atoms=0; state=False
    text=path.read_text()
    for line in text.splitlines():
        if line.startswith('Atoms'):
            state=True; continue
        if state:
            fields=line.split()
            if not fields: continue
            if not fields[0].isdigit(): break
            counts[int(fields[2])]+=1; charge+=float(fields[3]); atoms+=1
    return {'atoms':atoms,'net_charge_e':charge,'hydroxyl_oxygen_count':counts[12],
            'epoxy_conversion_percent':100*counts[12]/200}

def main():
    root=Path(sys.argv[1] if len(sys.argv)>1 else '.')
    log=(root/'log.lammps').read_text()
    state=parse_data(root/'pilot.data')
    timed=log.split('PILOT_REACTIVE_TIMED_BEGIN')[-1]
    matches=re.findall(r'Loop time of ([\deE.+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms',timed)
    timing=None
    if matches:
        seconds,ranks,steps,atoms=matches[-1]; seconds=float(seconds);steps=int(steps)
        layout=re.findall(r'with (\d+) MPI tasks x (\d+) OpenMP threads',timed)
        timing={'loop_seconds':seconds,'reported_total_workers':int(ranks),'steps':steps,'dt_fs':1.0,
                'mpi_ranks':int(layout[-1][0]) if layout else None,
                'openmp_threads_per_rank':int(layout[-1][1]) if layout else None,
                'ns_per_day_short_extrapolation':steps/seconds*86400/1e6}
    passed=('PILOT_COMPLETED' in log and state['atoms']==6350
            and math.isfinite(state['net_charge_e']) and abs(state['net_charge_e'])<1e-3
            and 'ERROR:' not in log)
    result={'execution_passed':passed,'reaction_observed':state['hydroxyl_oxygen_count']>0,
            'final_structure':state,'timing':timing,
            'warnings':sorted(set(re.findall(r'WARNING: ([^\n]+)',log))),
            'density_equilibrated':False,'tg_calculated':False,'bo_eligible':False,
            'scope':'Short compatibility/stability pilot. Cure conversion is not a kinetic prediction.'}
    (root/'pilot_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not passed: raise SystemExit(1)

if __name__=='__main__':main()
