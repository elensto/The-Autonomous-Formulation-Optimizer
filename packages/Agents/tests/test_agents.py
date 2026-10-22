import csv, hashlib, json, pathlib, sys, tempfile, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'agents'))
from core import read_csv,resolve_match,validate_structure,validate_result,validate_proposal
from experiment_agent import make_job,write_job
class AgentsTest(unittest.TestCase):
 def test_uncured_dds_is_not_runnable(self):
  row=read_csv(ROOT/'data/network_manifest.csv')[-6]
  self.assertEqual(row['record_kind'],'uncured_start')
  self.assertEqual(validate_structure(row,ROOT).decision,'REJECT')
 def test_real_manifest_no_verified_structure(self):
  m=read_csv(ROOT/'data/network_manifest.csv')
  self.assertEqual(len(m),33)
  for file in ['bo_proposal_D230.json','bo_proposal_DDS.json','bo_proposal_DDS_composition.json']:
   p=json.loads((ROOT/'examples'/file).read_text())
   row,error=resolve_match(p,m)
   self.assertIsNone(row);self.assertIsNotNone(error)
 def test_no_literature_as_observation(self):
  r=json.loads((ROOT/'examples/md_result_PENDING.json').read_text())
  self.assertEqual(validate_result(r).decision,'MODIFY')
  r.update(source_type='literature',property='Tg',fit_accepted=True,npt_converged=True,protocol_reviewed=True,structure_id='X',tg_K=455,tg_uncertainty_K=5,cooling_trace_path='a')
  self.assertIn('NOT_REAL_CLASSICAL_MD',validate_result(r).reasons)
 def test_proposal_validation(self):
  p=json.loads((ROOT/'examples/bo_proposal_DDS_composition.json').read_text())
  self.assertEqual(validate_proposal(p).decision,'ACCEPT')
  p['components']['hardener_id']='H006'
  self.assertEqual(validate_proposal(p).decision,'REJECT')
 def test_verified_fixture_can_prepare_job(self):
  with tempfile.TemporaryDirectory() as d:
   b=pathlib.Path(d);(b/'raw').mkdir();(b/'config').mkdir()
   data=b/'raw/network.data'
   data.write_text('''LAMMPS data test\n2 atoms\n1 bonds\n0 angles\n0 dihedrals\n0 impropers\n\nMasses\n\n1 12.0\n\nAtoms # full\n\n1 1 1 0.0 1 2 3\n2 1 1 0.0 2 3 4\n\nBonds\n\n1 1 1 2\n''')
   restart=b/'raw/equilibrated.restart';restart.write_text('MOCK RESTART NOT VALID LAMMPS: TEST ONLY')
   (b/'config/tmpl.lmp').write_text('units real\natom_style full\n{{READ_RESTART}}\ntimestep 1.0\n')
   row={'structure_id':'F1','system_id':'SYS_D230_FAN','formulation_id':'DGEBA_D230_NOMINAL','record_kind':'cured_network','verified_cure_percent':'80','replica_id':'1','relative_path':'raw/network.data','file_sha256':hashlib.sha256(data.read_bytes()).hexdigest(),'number_atoms':'2','topology_reviewed':'true','force_field_reviewed':'true','npt_equilibrated':'true','force_field_id':'FAN_CHARMM_REVIEWED','equilibrated_restart_path':'raw/equilibrated.restart'}
   proto={'protocol_id':'DEMO_TEST_ONLY','engine':'lammps','ensemble':'NPT','force_field_id':'FAN_CHARMM_REVIEWED','reviewed_by_scientist':True,'template_file':'config/tmpl.lmp','temperatures_K':[450,420,390,360,330,300,270],'timestep_fs':1.0,'pressure_atm':1.0,'thermostat_damp_fs':100,'barostat_damp_fs':1000,'equil_steps_per_T':100,'production_steps_per_T':200,'averaging_stride':10}
   proposal={'request_id':'TEST_01','mode':'cure_conversion','system_id':'SYS_D230_FAN','conversion_pct':80,'replica_id':'1','target_tg_C':85}
   report,bundle=make_job(proposal,[row],proto,b)
   self.assertEqual(report['status'],'ACCEPT',report)
   path=pathlib.Path(write_job(bundle))
   self.assertTrue((path/'in.tg_scan.lmp').is_file())
   self.assertEqual(make_job(proposal,[row],proto,b)[0]['status'],'MODIFY')
   # Not a real LAMMPS restart: this test verifies orchestration only.
 def test_missing_metadata_blocks(self):
  manifest=read_csv(ROOT/'data/network_manifest.csv')
  p=json.loads((ROOT/'examples/bo_proposal_D230.json').read_text())
  proto=json.loads((ROOT/'config/protocol_D230.json').read_text())
  report,bundle=make_job(p,manifest,proto,ROOT)
  self.assertEqual(report['status'],'MODIFY');self.assertIsNone(bundle)
if __name__=='__main__':unittest.main()
