#!/usr/bin/env python3
"""Check the recovered-configuration gate and profile campaign protocol."""
import ast
import contextlib
import copy
import importlib.util
import io
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('n640_profiles',str(HERE/'campaign.py'))
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)


class FakeCheck:
    RAW_COLUMNS=['configuration','tau','piter','axis','R','long','tr1','tr2','S_RAW','Nsrc']
    SU3_COLUMNS=['configuration','tau','piter','unit_max']
    avg=None
    def read_measurement(self,path,cfg,tau,piter):
        scan={(-1,4,i,0,0):dict(S_RAW='1.0',Nsrc='786432') for i in range(332)}
        return {},scan,[],dict(configuration=cfg,tau=tau,piter=piter,unit_max='1e-15')
    def raw_rows(self,cfg,tau,piter,scan):
        return [dict(configuration=cfg,tau=tau,piter=piter,axis=k[0],R=k[1],long=k[2],tr1=k[3],tr2=k[4],**v) for k,v in scan.items()]


class FakeTrace:
    TRACE_FIELDS=['T_TRACE','LT_RE','LT_IM']
    def trace_records(self,text,scan,avg):
        return {k:dict(T_TRACE='0.1',LT_RE='0.2',LT_IM='0') for k in scan}


class Tests(unittest.TestCase):
    def setUp(self):
        output=contextlib.redirect_stdout(io.StringIO());output.__enter__();self.addCleanup(output.__exit__,None,None,None)

    def test_reviewed_recovery_receipt(self):
        spec=c.read(HERE/'INPUTS.json');summary=c.read(HERE/'inputs/recovery_summary.json')
        for name,digest in spec['input_pins'].items():self.assertEqual(c.digest(HERE/name),digest)
        self.assertTrue(summary['complete_eigenpair_set_valid'])
        self.assertEqual((summary['configuration'],summary['times'],summary['audited_modes']),(640,192,19200))
        self.assertEqual(len(c.read(HERE/'inputs/recovery_saved_files.json')),193)
        self.assertFalse(summary['profile_measurements_added'])

    def test_exact_case_plan(self):
        self.assertEqual([r['tau'] for r in c.cases()],[2,4,6])
        self.assertEqual(sum(len(r['piters']) for r in c.cases()),12)
        self.assertTrue(all(r['configuration']==640 and r['piters']==[4,10,20,40] for r in c.cases()))

    def test_raw_trace_join_and_invariance(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for p in (4,10,20,40):(out/('p%d.out'%p)).write_text('fixture\n')
            trace=FakeTrace();check=FakeCheck();entry=c.cases()[0]
            raw,windows=c.read_case(out,entry,trace,check)
            self.assertEqual((len(raw),len(windows)),(1328,4))
            self.assertIn('T_TRACE',raw[0]);self.assertEqual(raw[0]['configuration'],'640')
            original=check.read_measurement
            def changed(path,cfg,tau,piter):
                a,scan,b,w=original(path,cfg,tau,piter)
                if piter==10:scan[(-1,4,0,0,0)]['S_RAW']='1.1'
                return a,scan,b,w
            check.read_measurement=changed
            with self.assertRaisesRegex(ValueError,'clover/source'):c.read_case(out,entry,trace,check)

    def test_recovery_gate_rejects_unreviewed_result(self):
        expected=c.read(HERE/'inputs/recovery_summary.json');altered=copy.deepcopy(expected);altered['audited_modes']=19100
        recovery=SimpleNamespace(prepared=lambda name:(Path('/tmp'),{}),final_replay=lambda f,p:altered)
        with patch.object(c,'dependencies',return_value=(dict(recovery_name='full192'),None,None,None)),patch.object(c,'module',return_value=recovery):
            with self.assertRaisesRegex(ValueError,'reviewed result'):c.recovery()

    def test_prepare_submission_and_immutability(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/profiles';here.mkdir(parents=True)
            exe=root/'existing_exe';exe.write_bytes(b'validated executable')
            receipt=c.read(HERE/'inputs/recovery_summary.json')
            spec=dict(settings={'Nv':100},scope='fixture')
            trace=FakeTrace();trace.EXE=exe;check=FakeCheck();args=SimpleNamespace(name='matched')
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'dependencies',return_value=(spec,trace,check,'manifest')),patch.object(c,'recovery',return_value=receipt):
                c.prepare(args);target,plan,_,_=c.prepared('matched')
                self.assertEqual(plan['cases'],c.cases())
                with patch.object(c.subprocess,'run',return_value=SimpleNamespace(returncode=0,stdout='12345\n',stderr='')) as run:
                    c.submit(args);self.assertIn('--array=0-2%3',run.call_args[0][0])
                    with self.assertRaises(ValueError):c.submit(args)
                (target/'PREPARED.json').write_text('{}\n')
                with self.assertRaisesRegex(ValueError,'plan changed'):c.prepared('matched')

    def test_run_all_cases_and_collect_protocol(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/profiles';here.mkdir(parents=True)
            target=here/'campaigns/matched';target.mkdir(parents=True)
            trace=FakeTrace();trace.EXE=root/'existing_exe';check=FakeCheck()
            receipt=c.read(HERE/'inputs/recovery_summary.json')
            plan=dict(cases=c.cases(),recovery=receipt,scope='fixture')
            (target/'PREPARED.json').write_text('{}\n');args=SimpleNamespace(name='matched')
            def fake_mpi(command,**kwargs):
                self.assertEqual(command[:5],['mpiexec','-n','64','-ppn','8'])
                env=kwargs['env'];self.assertEqual(env['QX_CFG'],'640');self.assertEqual(env['EVROOT'],receipt['eigenvector_directory'])
                self.assertEqual(int(env['T0'])+int(env['TAU_LENGTH'])//2,98)
                kwargs['stdout'].write('fixture\n');return SimpleNamespace(returncode=0)
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'prepared',return_value=(target,plan,trace,check)),patch.object(c.subprocess,'run',side_effect=fake_mpi):
                for slot in range(3):
                    with patch.dict(os.environ,{'SLURM_NTASKS':'64','SLURM_JOB_NUM_NODES':'8','SLURM_ARRAY_TASK_ID':str(slot)}):c.run(args)
                c.collect(args);summary=c.read(target/'SUMMARY.json')
                self.assertTrue(summary['success']);self.assertEqual(summary['raw_records'],3984)
                self.assertEqual(summary['window_records'],12);self.assertFalse(summary['ensemble_71_recomputed'])
                self.assertTrue((root/'qx_n640_profiles_matched_review.zip').is_file())
                (target/'case_000/p4.out').write_text('changed')
                with self.assertRaisesRegex(ValueError,'output changed'):c.launch_records(target/'case_000',c.cases()[0])

    def test_incomplete_cases_export_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/profiles';here.mkdir(parents=True)
            target=here/'campaigns/matched';target.mkdir(parents=True)
            (target/'failure.err').write_text('failed attempt\n')
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'prepared',side_effect=ValueError('case incomplete')):
                c.collect(SimpleNamespace(name='matched'))
            self.assertFalse(c.read(target/'SUMMARY.json')['success'])
            import zipfile
            with zipfile.ZipFile(root/'qx_n640_profiles_matched_review.zip') as z:self.assertIn('campaign/failure.err',z.namelist())

    def test_name_and_syntax(self):
        for name in ('../x','x/y','x\n',''):
            with self.assertRaises(ValueError):c.folder(name)
        ast.parse((HERE/'campaign.py').read_text())
        subprocess.run(['bash','-n',str(HERE/'submit_micro.sh')],check=True)


if __name__=='__main__':unittest.main()
