import csv
import json
import pathlib
import sys
import unittest
import numpy as np
import pandas as pd
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from dds_optimizer import ratio_from_counts,build_composition_pool,propose,read_csv,REQUIRED_CAND,REQUIRED_OBS,accepted_observations

class DDSBOTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.empty=read_csv(ROOT/'data/md_tg_observations.csv',REQUIRED_OBS)
        cls.fixture=ROOT/'examples/DEMO_SYNTHETIC_ONLY'
        cls.demo_obs=read_csv(cls.fixture/'DEMO_observations.csv',REQUIRED_OBS)
        cls.demo_pool=read_csv(cls.fixture/'DEMO_candidates.csv',REQUIRED_CAND)

    def test_ratio_and_limit(self):
        self.assertEqual(ratio_from_counts(100,50),1.0)
        self.assertEqual(ratio_from_counts(100,40),0.8)
        pool=build_composition_pool()
        self.assertTrue((pool.conversion_pct<=np.minimum(1,pool.ratio)*100).all())
        self.assertFalse(((pool.ratio==0.8)&(pool.conversion_pct==90)).any())
        self.assertEqual(len(pool),18)

    def test_real_data_empty(self):
        real=read_csv(ROOT/'data/dds_network_manifest.csv',REQUIRED_CAND)
        out,m=propose(real,self.empty,'conversion')
        self.assertEqual(m['status'],'NO_ELIGIBLE_CANDIDATES')
        self.assertEqual(len(out),0)

    def test_real_composition_is_exploration(self):
        pool=read_csv(ROOT/'data/dds_composition_candidate_pool.csv',REQUIRED_CAND)
        out,m=propose(pool,self.empty,'composition')
        self.assertEqual(m['method'],'diversity_exploration')
        self.assertEqual(len(out),3)
        self.assertTrue(out.proposed_job_type.eq('generate_and_cure_network').all())
        self.assertTrue(out.predicted_tg_C.isna().all())

    def test_synthetic_not_approved_without_explicit_flag(self):
        self.assertEqual(len(accepted_observations(self.demo_obs)),0)
        out,m=propose(self.demo_pool,self.demo_obs,'conversion')
        self.assertEqual(m['n_approved_observations'],0)
        self.assertEqual(m['method'],'diversity_exploration')

    def test_synthetic_gp_only_when_opted_in(self):
        out,m=propose(self.demo_pool,self.demo_obs,'conversion',allow_synthetic_demo=True)
        self.assertEqual(m['method'],'GP_target_interval_probability')
        self.assertEqual(set(out.conversion_pct),{85,90})
        self.assertTrue(np.isfinite(out.predicted_tg_C).all())
        self.assertTrue(((out.p_tg_in_window>=0)&(out.p_tg_in_window<=1)).all())

    def test_duplicate_pending_excluded(self):
        pool=read_csv(ROOT/'data/dds_composition_candidate_pool.csv',REQUIRED_CAND)
        pool.loc[pool.index[0],'job_status']='running'
        out,m=propose(pool,self.empty,'composition',batch_size=18)
        self.assertNotIn(pool.iloc[0].candidate_id,set(out.candidate_id))

    def test_invalid_ratio_rejected(self):
        pool=read_csv(ROOT/'data/dds_composition_candidate_pool.csv',REQUIRED_CAND)
        pool.loc[0,'ratio']=0.99
        with self.assertRaises(ValueError):propose(pool,self.empty,'composition')

if __name__=='__main__':unittest.main()
