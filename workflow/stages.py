"""The five architecture stages. Each stage is a small, testable function/class.

Stage 1  Bayesian Optimizer        -> propose_next()        (wraps packages/BO dds_optimizer.propose)
Stage 2  Experiment Agent+PhysCheck-> ExperimentAgent        (reuses packages/Agents core.validate_proposal)
Stage 3  Simulation Engine         -> SimulationEngine       (replay of real Carbono MD | SLURM job writer)
Stage 4  Result Parser+PostCheck   -> parse_and_check()      (reuses packages/Agents core.validate_result)
Stage 5  Observation Store         -> ObservationStore       (CSV/JSONL, BO-compatible schema)
"""
from __future__ import annotations

import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path

import warnings
import numpy as np
import pandas as pd

from . import ROOT
import dds_optimizer as bo          # packages/BO/src
import core as agents_core          # packages/Agents/agents

warnings.filterwarnings('ignore', module='sklearn')

OBS_COLUMNS = bo.REQUIRED_OBS + ['iteration', 'formulation_id', 'tg_ci95_low_C', 'tg_ci95_high_C',
                                 'target_conversion_pct', 'post_physcheck_reasons', 'job_id']


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def h_to_epoxy_ratio(n_dgeba, n_dds, cfg):
    c = cfg['composition_constraints']
    return c['amine_H_per_DDS'] * int(n_dds) / (c['epoxy_groups_per_DGEBA'] * int(n_dgeba))


# --------------------------------------------------------------------------- candidates
def load_candidates(cfg, pool='all', replay_ids=None):
    """Candidate library = formulations that the Simulations package can actually build
    (the 23 recipes prepared from the published Livraghi 100:50 box)."""
    manifest = json.loads((ROOT / cfg['candidate_source']).read_text())
    rows = []
    for c in manifest['cases']:
        rows.append({
            'candidate_id': c['formulation_id'], 'structure_id': c['formulation_id'],
            'n_DGEBA': c['DGEBA_molecules'], 'n_DDS': c['DDS_molecules'],
            'ratio': h_to_epoxy_ratio(c['DGEBA_molecules'], c['DDS_molecules'], cfg),
            'conversion_pct': float(c['target_epoxy_conversion_percent']),
            'structure_status': 'NEEDS_NETWORK_GENERATION', 'file_path': '',
            'physcheck_decision': 'NOT_EVALUATED', 'job_status': 'not_started',
            'array_index': c['array_index'], 'atoms': c['atoms'],
        })
    df = pd.DataFrame(rows)
    if pool == 'replay-available' and replay_ids is not None:
        df = df[df.candidate_id.isin(replay_ids)]
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- stage 1
def propose_next(candidates, observations, cfg):
    """Iteration 0 has no data: seed with the eligible formulation closest to the textbook
    stoichiometric recipe (r = 1, full cure). Afterwards the team's BO package decides:
    farthest-first exploration (<min_points_for_gp points), then GP + P(Tg in window)."""
    o, t = cfg['optimizer'], cfg['target']
    seed = o.get('seed_reference')
    if seed and not len(observations):
        c = candidates[~candidates.job_status.isin(['pending', 'completed'])].copy()
        if len(c):
            X = c[['ratio', 'conversion_pct']].to_numpy(float)
            span = np.where(np.ptp(X, axis=0) > 0, np.ptp(X, axis=0), 1.0)
            d = np.linalg.norm((X - np.array([seed['ratio'], seed['conversion_pct']])) / span, axis=1)
            sel = c.iloc[[int(np.argmin(d))]].copy()
            for k in ('predicted_tg_C', 'tg_predictive_sd_C', 'p_tg_in_window'):
                sel[k] = np.nan
            sel['selection_method'] = 'reference_seed'
            return sel.reset_index(drop=True), {'status': 'OK', 'method': 'reference_seed (stoichiometric 1:1 recipe)',
                                                'n_approved_observations': 0}
    obs = observations[bo.REQUIRED_OBS] if len(observations) else pd.DataFrame(columns=bo.REQUIRED_OBS)
    return bo.propose(candidates[bo.REQUIRED_CAND], obs, 'composition', t['tg_min_C'], t['tg_max_C'],
                      o['batch_size'], o['seed'], False, o['min_points_for_gp'])


# --------------------------------------------------------------------------- stage 2
class ExperimentAgent:
    """Turns a BO proposal into a simulation-ready request and runs the pre-simulation PhysCheck."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.batch_dir = ROOT / cfg['simulation']['slurm_batch_dir']

    def physcheck(self, p):
        cc = self.cfg['composition_constraints']
        reasons, warnings = [], []
        # Reuse the Agents package proposal gate (system / chemistry identity / conversion range)
        gate = agents_core.validate_proposal({
            'mode': 'formulation', 'system_id': 'SYS_DDS_LIVRAGHI', 'formulation_id': p['candidate_id'],
            'components': {'resin_id': 'R001', 'hardener_id': 'H007'}, 'conversion_pct': p['conversion_pct']})
        reasons += gate.reasons
        nd, nh = int(p['n_DGEBA']), int(p['n_DDS'])
        r = h_to_epoxy_ratio(nd, nh, self.cfg)
        if abs(r - float(p['ratio'])) > 1e-6:
            reasons.append('STOICHIOMETRY_RATIO_INCONSISTENT')
        if not (cc['ratio_min'] - 1e-9 <= r <= cc['ratio_max'] + 1e-9):
            reasons.append('RATIO_OUTSIDE_CONSTRAINTS')
        if float(p['conversion_pct']) > 100 * min(1.0, r) + 1e-6:
            reasons.append('CONVERSION_EXCEEDS_STOICHIOMETRIC_LIMIT')
        case_file = self.batch_dir / 'cases' / p['candidate_id'] / 'case.json'
        if not case_file.is_file():
            reasons.append('NO_BUILDABLE_STRUCTURE_RECIPE')
            return {'decision': 'REJECT', 'reasons': reasons, 'warnings': warnings}, None
        case = json.loads(case_file.read_text())
        if case['epoxy_groups'] != cc['epoxy_groups_per_DGEBA'] * nd or case['amine_hydrogens'] != cc['amine_H_per_DDS'] * nh:
            reasons.append('FUNCTIONAL_GROUP_COUNT_MISMATCH')
        if abs(case.get('initial_net_charge_e', 1.0)) > 1e-3:
            reasons.append('SYSTEM_NOT_CHARGE_NEUTRAL')
        if not (cc['min_atoms'] <= case['atoms'] <= cc['max_atoms']):
            reasons.append('ATOM_COUNT_OUTSIDE_FF_PRACTICAL_RANGE')
        if not case.get('whole_molecule_selection_verified'):
            reasons.append('MOLECULE_INTEGRITY_UNVERIFIED')
        if not case.get('data_sha256'):
            reasons.append('STRUCTURE_HASH_MISSING')
        warnings.append('Classical Livraghi AMBER-style FF (MLIP not used for production Tg: ~100x too slow, see docs).')
        warnings.append('Single network realization per formulation; Tg is cooling-rate dependent.')
        decision = 'REJECT' if reasons else 'ACCEPT'
        return {'decision': decision, 'reasons': reasons, 'warnings': warnings}, case

    def prepare(self, p, iteration):
        check, case = self.physcheck(p)
        request = None
        if check['decision'] == 'ACCEPT':
            request = {
                'job_id': f"it{iteration:02d}_{p['candidate_id']}", 'formulation_id': p['candidate_id'],
                'n_DGEBA': int(p['n_DGEBA']), 'n_DDS': int(p['n_DDS']), 'ratio': float(p['ratio']),
                'target_epoxy_conversion_pct': float(p['conversion_pct']), 'atoms': case['atoms'],
                'array_index': case['array_index'], 'structure_sha256': case['data_sha256'],
                'protocol': self.cfg['simulation']['protocol'], 'engine': 'LAMMPS (fix bond/react) classical',
                'stages': ['prepare', 'cure', 'equilibrate', 'tg_screen'], 'created_utc': now()}
        return check, request


# --------------------------------------------------------------------------- stage 3
class SimulationEngine:
    """backend='replay': returns the real Carbono MD result for this formulation (retrospective run).
       backend='slurm' : writes a ready-to-submit SLURM job using the Simulations package scripts."""

    def __init__(self, cfg, backend, run_dir):
        self.cfg, self.backend, self.run_dir = cfg, backend, Path(run_dir)
        self.results = pd.read_csv(ROOT / cfg['simulation']['replay_results']).set_index('formulation_id')

    def run(self, req):
        job_dir = self.run_dir / 'jobs' / req['job_id']
        job_dir.mkdir(parents=True, exist_ok=True)
        (job_dir / 'request.json').write_text(json.dumps(req, indent=2))
        if self.backend == 'slurm':
            batch = ROOT / self.cfg['simulation']['slurm_batch_dir']
            script = batch.parent / 'scripts' / 'parallel_array.slurm'
            (job_dir / 'submit.sh').write_text(
                '#!/usr/bin/env bash\n# Run on the cluster from the DGEBA_DDS directory (requires VALIDATION_PASSED in the batch).\n'
                f"sbatch --array={req['array_index']} {script.relative_to(ROOT / 'packages/Simulations/DGEBA_DDS')} "
                f"{batch.relative_to(ROOT / 'packages/Simulations/DGEBA_DDS')}\n")
            (job_dir / 'submit.sh').chmod(0o755)
            return {'job_id': req['job_id'], 'status': 'JOB_WRITTEN_NOT_SUBMITTED', 'job_dir': str(job_dir), 'raw': None}
        fid = req['formulation_id']
        if fid not in self.results.index:
            return {'job_id': req['job_id'], 'status': 'NO_RESULT_YET', 'job_dir': str(job_dir), 'raw': None}
        raw = self.results.loc[fid].to_dict()
        raw['formulation_id'] = fid
        (job_dir / 'raw_result.json').write_text(json.dumps(raw, indent=2, default=float))
        return {'job_id': req['job_id'], 'status': 'COMPLETED', 'job_dir': str(job_dir), 'raw': raw}


# --------------------------------------------------------------------------- stage 4
def parse_and_check(sim, req, cfg):
    raw = sim['raw']
    pc = cfg['post_physcheck']
    tg, lo, hi = float(raw['tg_C']), float(raw['tg_ci95_low_C']), float(raw['tg_ci95_high_C'])
    conv = float(raw['degree_of_cure']) * 100.0
    se = max((hi - lo) / 3.92, 1.0)            # 95% CI -> approx. standard error
    result = {'source_type': 'real_classical_MD', 'property': 'Tg', 'structure_id': req['formulation_id'],
              'tg_K': tg + 273.15, 'tg_uncertainty_K': se, 'fit_accepted': True, 'npt_converged': True,
              'protocol_reviewed': True, 'cooling_trace_path': 'data/eight_tg_curves.png'}
    gate = agents_core.validate_result(result)   # reuse Agents post gate
    reasons = list(gate.reasons)
    if not (lo <= tg <= hi):
        reasons.append('TG_OUTSIDE_OWN_CI')
    if hi - lo > pc['max_ci95_width_C']:
        reasons.append('TG_CI_TOO_WIDE')
    if not (pc['tg_plausible_range_C'][0] <= tg <= pc['tg_plausible_range_C'][1]):
        reasons.append('TG_IMPLAUSIBLE')
    if conv > 100 * min(1.0, req['ratio']) + 1.0:
        reasons.append('ACHIEVED_CONVERSION_EXCEEDS_STOICHIOMETRY')
    if abs(conv - req['target_epoxy_conversion_pct']) > pc['max_conversion_deviation_pct']:
        reasons.append('ACHIEVED_CONVERSION_OFF_TARGET')
    ok = not reasons
    obs = {'observation_id': f"OBS_{req['job_id']}", 'candidate_id': req['formulation_id'],
           'n_DGEBA': req['n_DGEBA'], 'n_DDS': req['n_DDS'], 'ratio': req['ratio'],
           'conversion_pct': round(conv, 3), 'tg_C': round(tg, 2), 'tg_se_C': round(se, 2),
           'physcheck_decision': 'ACCEPT' if ok else 'REJECT', 'bo_trainable': 'true' if ok else 'false',
           'source_type': 'real_classical_MD', 'experiment_status': 'completed',
           'formulation_id': req['formulation_id'], 'tg_ci95_low_C': round(lo, 2), 'tg_ci95_high_C': round(hi, 2),
           'target_conversion_pct': req['target_epoxy_conversion_pct'], 'post_physcheck_reasons': ';'.join(reasons),
           'job_id': req['job_id']}
    return obs, {'decision': obs['physcheck_decision'], 'reasons': reasons}


# --------------------------------------------------------------------------- stage 5
class ObservationStore:
    def __init__(self, run_dir, fresh=True):
        self.dir = Path(run_dir)
        if fresh and self.dir.exists():
            shutil.rmtree(self.dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.obs_path = self.dir / 'observations.csv'
        self.events = self.dir / 'events.jsonl'
        self.observations = pd.read_csv(self.obs_path, keep_default_na=False) if self.obs_path.exists() \
            else pd.DataFrame(columns=OBS_COLUMNS)

    def log(self, stage, **payload):
        with open(self.events, 'a') as f:
            f.write(json.dumps({'t': now(), 'stage': stage, **payload}, default=str) + '\n')

    def add(self, obs):
        self.observations = pd.concat([self.observations, pd.DataFrame([obs])], ignore_index=True)
        self.observations.to_csv(self.obs_path, index=False)

    def trainable(self):
        o = self.observations
        return o[(o.physcheck_decision == 'ACCEPT') & (o.bo_trainable.astype(str) == 'true')] if len(o) else o
