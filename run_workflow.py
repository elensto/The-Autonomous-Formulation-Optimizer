#!/usr/bin/env python
"""Closed-loop autonomous formulation optimizer (DGEBA / 4,4'-DDS, target Tg window).

    python run_workflow.py                       # retrospective closed loop on real Carbono MD results
    python run_workflow.py --pool all            # full 23-formulation library (unmeasured -> NO_RESULT_YET)
    python run_workflow.py --backend slurm       # write real SLURM jobs instead of replaying results
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from workflow import ROOT
from workflow.stages import (ExperimentAgent, ObservationStore, SimulationEngine, load_candidates,
                             parse_and_check, propose_next, now)


def banner(txt):
    print('\n' + '=' * 78 + f'\n{txt}\n' + '=' * 78)


def run(args):
    cfg = json.loads(Path(args.config).read_text())
    if args.tg_min is not None: cfg['target']['tg_min_C'] = args.tg_min
    if args.tg_max is not None: cfg['target']['tg_max_C'] = args.tg_max
    if args.iterations: cfg['optimizer']['max_iterations'] = args.iterations
    if args.batch_size: cfg['optimizer']['batch_size'] = args.batch_size
    backend = args.backend or cfg['simulation']['backend']
    run_dir = ROOT / cfg['store_dir'] / args.run_name

    store = ObservationStore(run_dir, fresh=not args.resume)
    sim = SimulationEngine(cfg, backend, run_dir)
    agent = ExperimentAgent(cfg)
    pool = args.pool or ('replay-available' if backend == 'replay' else 'all')
    candidates = load_candidates(cfg, pool, list(sim.results.index))
    t = cfg['target']
    banner(f"TARGET  Tg in [{t['tg_min_C']}, {t['tg_max_C']}] C | backend={backend} | pool={pool} "
           f"({len(candidates)} formulations) | run={run_dir.relative_to(ROOT)}")
    store.log('config', config=cfg, backend=backend, pool=pool, n_candidates=len(candidates))

    hit = None
    for it in range(1, cfg['optimizer']['max_iterations'] + 1):
        banner(f'ITERATION {it}')
        # ---- 1. Bayesian optimizer
        proposals, meta = propose_next(candidates, store.observations, cfg)
        print(f"[1 BO] method={meta.get('method', meta['status'])} | trained on "
              f"{meta.get('n_approved_observations', 0)} accepted MD observations")
        store.log('bo', iteration=it, meta=meta, proposals=proposals.to_dict('records'))
        if meta['status'] != 'OK':
            print(f"[1 BO] stop: {meta['status']} - {meta.get('message')}")
            break
        for _, p in proposals.iterrows():
            pred = '' if np.isnan(p.predicted_tg_C) else \
                f" | GP Tg={p.predicted_tg_C:.0f}+/-{p.tg_predictive_sd_C:.0f} C, P(in window)={p.p_tg_in_window:.2f}"
            print(f"[1 BO] propose {p.candidate_id}: DGEBA={p.n_DGEBA} DDS={p.n_DDS} r={p.ratio:.3f} "
                  f"cure={p.conversion_pct:.1f}%{pred}")
            # ---- 2. Experiment agent + PhysCheck
            check, req = agent.prepare(p, it)
            store.log('physcheck_pre', iteration=it, candidate=p.candidate_id, check=check, request=req)
            print(f"[2 AGENT] pre-PhysCheck {check['decision']} {check['reasons'] or ''}")
            idx = candidates.index[candidates.candidate_id == p.candidate_id]
            if check['decision'] != 'ACCEPT':
                candidates = candidates.drop(idx).reset_index(drop=True)
                continue
            # ---- 3. Simulation engine
            res = sim.run(req)
            store.log('simulation', iteration=it, result=res)
            print(f"[3 SIM] {req['job_id']} -> {res['status']} ({res['job_dir'].split(str(ROOT) + '/')[-1]})")
            if res['status'] != 'COMPLETED':
                candidates.loc[idx, 'job_status'] = 'pending'
                continue
            # ---- 4. Result parser + post-PhysCheck
            obs, post = parse_and_check(res, req, cfg)
            obs['iteration'] = it
            store.log('physcheck_post', iteration=it, candidate=p.candidate_id, check=post, observation=obs)
            inwin = t['tg_min_C'] <= obs['tg_C'] <= t['tg_max_C']
            print(f"[4 PARSE] Tg={obs['tg_C']:.1f} C (95% CI {obs['tg_ci95_low_C']:.0f}-{obs['tg_ci95_high_C']:.0f}), "
                  f"cure achieved={obs['conversion_pct']:.1f}% | post-PhysCheck {post['decision']} "
                  f"{post['reasons'] or ''} | {'IN TARGET WINDOW' if inwin else 'outside target'}")
            # ---- 5. Observation store -> fed back to BO next iteration
            store.add(obs)
            candidates.loc[idx, 'job_status'] = 'completed'
            print(f"[5 STORE] {len(store.observations)} observations ({len(store.trainable())} BO-trainable)")
            if inwin and post['decision'] == 'ACCEPT' and hit is None:
                hit = {**obs, 'found_at_iteration': it}
        if hit and args.stop_on_target:
            print(f"\nTarget reached at iteration {it}: {hit['candidate_id']} Tg={hit['tg_C']} C -> stopping.")
            break

    summary = write_report(store, cfg, hit, len(candidates), backend, pool, run_dir)
    banner('SUMMARY')
    print(json.dumps(summary, indent=2))
    return summary


def write_report(store, cfg, hit, n_pool, backend, pool, run_dir):
    obs = store.observations
    t = cfg['target']
    n_hits = int(((obs.tg_C.astype(float) >= t['tg_min_C']) & (obs.tg_C.astype(float) <= t['tg_max_C'])).sum()) if len(obs) else 0
    summary = {'finished_utc': now(), 'backend': backend, 'pool': pool, 'target_C': [t['tg_min_C'], t['tg_max_C']],
               'n_simulations_used': int(len(obs)), 'n_bo_trainable': int(len(store.trainable())),
               'n_in_window': n_hits, 'best': hit,
               'note': ('Retrospective closed loop: the simulation stage replays REAL Carbono LAMMPS results '
                        '(1 network realization each, cooling-rate-dependent preliminary Tg).' if backend == 'replay'
                        else 'SLURM job folders written; results enter the loop when parsed.')}
    if backend == 'replay' and pool == 'replay-available':
        res = pd.read_csv(ROOT / cfg['simulation']['replay_results'])
        k = int(((res.tg_C >= t['tg_min_C']) & (res.tg_C <= t['tg_max_C'])).sum())
        if k:
            summary['random_search_expected_evals_to_first_hit'] = round((len(res) + 1) / (k + 1), 2)
    (run_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    if len(obs):
        try:
            import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
            o = obs.copy(); o['tg_C'] = o.tg_C.astype(float)
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            ax[0].axhspan(t['tg_min_C'], t['tg_max_C'], color='tab:green', alpha=.15, label='target window')
            ax[0].errorbar(range(1, len(o) + 1), o.tg_C, yerr=[o.tg_C - o.tg_ci95_low_C.astype(float),
                           o.tg_ci95_high_C.astype(float) - o.tg_C], fmt='o-', capsize=3)
            for i, r in enumerate(o.itertuples(), 1): ax[0].annotate(r.candidate_id, (i, r.tg_C), fontsize=7, xytext=(4, 4), textcoords='offset points')
            ax[0].set(xlabel='simulation #', ylabel='Tg (°C)', title='Closed-loop trajectory'); ax[0].legend()
            sc = ax[1].scatter(o.ratio.astype(float), o.conversion_pct.astype(float), c=o.tg_C, cmap='viridis', s=80)
            fig.colorbar(sc, ax=ax[1], label='Tg (°C)')
            ax[1].set(xlabel='amine-H / epoxy ratio', ylabel='achieved epoxy conversion (%)', title='Explored design space')
            fig.tight_layout(); fig.savefig(run_dir / 'trajectory.png', dpi=130); plt.close(fig)
        except Exception as e:  # plotting is optional
            print('plot skipped:', e)
    return summary


def cli():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config', default=str(ROOT / 'config/experiment.json'))
    p.add_argument('--backend', choices=['replay', 'slurm'])
    p.add_argument('--pool', choices=['replay-available', 'all'])
    p.add_argument('--iterations', type=int)
    p.add_argument('--batch-size', type=int)
    p.add_argument('--tg-min', type=float)
    p.add_argument('--tg-max', type=float)
    p.add_argument('--run-name', default='demo')
    p.add_argument('--resume', action='store_true', help='keep existing observations in runs/<run-name>')
    p.add_argument('--stop-on-target', action='store_true')
    return p


if __name__ == '__main__':
    run(cli().parse_args())
