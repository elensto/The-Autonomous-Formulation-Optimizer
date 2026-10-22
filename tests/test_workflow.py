import json, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pandas as pd
from workflow.stages import ExperimentAgent, load_candidates, parse_and_check
import run_workflow

CFG = json.loads((ROOT / 'config/experiment.json').read_text())


class TestStages(unittest.TestCase):
    def test_library_has_23_buildable_formulations(self):
        c = load_candidates(CFG)
        self.assertEqual(len(c), 23)
        self.assertTrue(((c.ratio >= 0.8) & (c.ratio <= 1.2)).all())

    def test_pre_physcheck_accepts_valid_and_rejects_bad(self):
        agent = ExperimentAgent(CFG)
        good = load_candidates(CFG).iloc[0].copy()
        self.assertEqual(agent.physcheck(good)[0]['decision'], 'ACCEPT')
        bad = good.copy(); bad['ratio'] = 1.5
        self.assertIn('STOICHIOMETRY_RATIO_INCONSISTENT', agent.physcheck(bad)[0]['reasons'])
        over = good.copy(); over['conversion_pct'] = 95.0       # r=0.8 -> max 80 %
        self.assertIn('CONVERSION_EXCEEDS_STOICHIOMETRIC_LIMIT', agent.physcheck(over)[0]['reasons'])
        ghost = good.copy(); ghost['candidate_id'] = 'D050_H050'
        self.assertEqual(agent.physcheck(ghost)[0]['decision'], 'REJECT')

    def test_post_physcheck_rejects_wide_ci(self):
        req = {'job_id': 'x', 'formulation_id': 'X', 'n_DGEBA': 100, 'n_DDS': 45, 'ratio': 0.9, 'target_epoxy_conversion_pct': 81}
        raw = {'tg_C': 200, 'tg_ci95_low_C': 100, 'tg_ci95_high_C': 300, 'degree_of_cure': 0.81}
        obs, chk = parse_and_check({'raw': raw}, req, CFG)
        self.assertEqual(chk['decision'], 'REJECT'); self.assertEqual(obs['bo_trainable'], 'false')

    def test_closed_loop_finds_target(self):
        args = run_workflow.cli().parse_args(['--run-name', 'unittest', '--stop-on-target'])
        s = run_workflow.run(args)
        self.assertIsNotNone(s['best'])
        self.assertTrue(150 <= s['best']['tg_C'] <= 220)
        self.assertLessEqual(s['n_simulations_used'], 8)


if __name__ == '__main__':
    unittest.main()
