import json, pathlib, sys, tempfile, hashlib, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'agents'))
from curing_agent import validate_curing_request, prepare_curing_job,write_curing_job,assess_curing_campaign

class CuringAgentTests(unittest.TestCase):
 def setUp(self):
  self.base=json.loads((ROOT/'config/curing_DDS_EXTERNAL.json').read_text())
  self.proposal=json.loads((ROOT/'examples/bo_proposal_DDS_cure_request.json').read_text())
 def test_default_is_blocked(self):
  d=validate_curing_request(self.proposal,self.base,ROOT)
  self.assertEqual(d['decision'],'MODIFY');self.assertIn('CURING_WORKFLOW_NOT_SCIENTIST_REVIEWED',d['reasons'])
 def test_report_does_not_count_target_as_achieved(self):
  e=json.loads((ROOT/'data/dds_campaign_status_2026-10-21.json').read_text())
  d=assess_curing_campaign(e);self.assertEqual(d['decision'],'MODIFY')
  self.assertIsNone(d['achieved_conversion_pct'])
 def test_new_ratio_requires_new_workflow(self):
  prop=json.loads((ROOT/'examples/bo_proposal_DDS_new_ratio.json').read_text())
  d=validate_curing_request(prop,self.base,ROOT)
  self.assertIn('NEW_COMPOSITION_REQUIRES_REBUILT_START_AND_REVIEWED_REACTION_DECK',d['reasons'])
 def test_d230_cannot_be_cured_with_dds_script(self):
  prop={**self.proposal,'system_id':'SYS_D230_FAN'}
  self.assertIn('CURING_BACKEND_ONLY_SUPPORTS_DGEBA_DDS',validate_curing_request(prop,self.base,ROOT)['reasons'])
 def test_verified_fixture_prepares_job_only(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=pathlib.Path(tmp);(base/'raw').mkdir()
   for name in ('start.data','cure.lmp'):(base/'raw'/name).write_text('TEST FILE NOT A VALID LAMMPS INPUT')
   cfg={**self.base,'scientist_reviewed':True,'input_script':'raw/cure.lmp','starting_structure':'raw/start.data'}
   cfg['input_script_sha256']=hashlib.sha256((base/'raw/cure.lmp').read_bytes()).hexdigest()
   cfg['starting_structure_sha256']=hashlib.sha256((base/'raw/start.data').read_bytes()).hexdigest()
   gate,bundle=prepare_curing_job(self.proposal,cfg,base)
   self.assertEqual(gate['decision'],'ACCEPT',gate)
   dest=write_curing_job(bundle)
   self.assertTrue((dest/'curing_job.json').exists())
   job=json.loads((dest/'curing_job.json').read_text())
   self.assertEqual(job['status'],'PREPARED_NOT_SUBMITTED')
   self.assertIsNone(job['achieved_conversion_pct'])
   self.assertEqual(prepare_curing_job(self.proposal,cfg,base)[0]['decision'],'MODIFY')
if __name__=='__main__':unittest.main()
