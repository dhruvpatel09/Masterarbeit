#!/usr/bin/env python3
"""Native contractions, complete synthetic replay and campaign failure tests."""
import argparse
import ast
import contextlib
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import campaign as c
import check_outputs as check
import generate_sources as helper
from matrix_plan import (PILOT,CUTOFFS,RADII,TAUS,AUDIT_TIMES,
                         expected_keys,direct_keys,legacy_keys,block_keys)
sys.path.insert(0,str(c.HERE/'tests'))
import synthetic_outputs as synthetic


class OutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs=synthetic.references([650])
        cls.text=synthetic.scan_fixture(650,10)

    def replay(self,text,refs=None,piter=10):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'scan.out';path.write_text(text)
            return check.scan(path,650,piter,self.refs if refs is None else refs)

    def test_complete_matrix_retains_ordered_pairs_and_signed_complex_values(self):
        values,blocks,summary=self.replay(self.text)
        self.assertEqual((len(expected_keys()),len(direct_keys()),len(legacy_keys()),len(block_keys())*64),(8448,1536,384,12288))
        self.assertEqual((len(values),len(blocks),summary['diagonal_reference_rows']),(8448,12288,2112))
        self.assertEqual(summary['equal_time_gram_checks'],16)
        self.assertTrue(any(float(r['L_RE'])<0 for r in values))
        self.assertTrue(any(float(r['L_IM'])!=0 for r in values))
        self.assertNotEqual(synthetic.value(650,40,100,4,12,0),synthetic.value(650,100,40,4,12,0).conjugate())
        high,_,_=self.replay(synthetic.scan_fixture(650,40),piter=40)
        self.assertEqual(check.projection_pair(values,high),0.)
        bad=[dict(r) for r in high];bad[100]['L_IM']='1'
        with self.assertRaises(ValueError):check.projection_pair(values,bad)

    def test_missing_duplicate_matrix_pair_and_wrong_orientation_rejected(self):
        first=next(x for x in self.text.splitlines(True) if x.startswith('MATRIX_RT '))
        variants=[self.text.replace(first,'',1),self.text+first,
            self.text.replace('source_cutoffs=40,60,80,100','source_cutoffs=40,80,60,100',1),
            synthetic.replace_record(self.text,'MATRIX_RT',{'Nv_src':'60'})]
        for text in variants:
            with self.subTest(text_length=len(text)):
                with self.assertRaises(ValueError):self.replay(text)

    def test_axis_average_count_and_diagonal_reference_checks(self):
        variants=[synthetic.replace_record(self.text,'MATRIX_RT',{'Nsrc':'1'}),
            synthetic.replace_record(self.text,'MATRIX_RT',{'L_RE':'0.5'},lambda r:r['axis']=='-1'),
            synthetic.replace_record(self.text,'MATRIX_RT',{'L_RE':'nan'})]
        for text in variants:
            with self.assertRaises(ValueError):self.replay(text)
        refs=dict(self.refs);key=(650,100,32,32,-1);refs[key]=dict(refs[key],L_RE='1')
        with self.assertRaises(ValueError):self.replay(self.text,refs)

    def test_direct_checks_are_recomputed_and_cannot_inflate_tolerance(self):
        direct=next(x for x in self.text.splitlines(True) if x.startswith('MATRIX_DIRECT '))
        legacy=next(x for x in self.text.splitlines(True) if x.startswith('MATRIX_LEGACY '))
        variants=[self.text.replace(direct,'',1),self.text.replace(legacy,'',1),
            synthetic.replace_record(self.text,'MATRIX_DIRECT',{'tolerance':'1'}),
            synthetic.replace_record(self.text,'MATRIX_DIRECT',{'direct_re':'1','abs_error':'0'}),
            synthetic.replace_record(self.text,'MATRIX_DIRECT',{'rank':'1'}),
            direct.replace('status=PASS','status=PASS status=PASS')+self.text]
        for text in variants:
            with self.assertRaises(ValueError):self.replay(text)

    def test_projection_roundtrip_and_rank_block_failure_detection(self):
        block=next(x for x in self.text.splitlines(True) if x.startswith('MATRIX_BLOCK '))
        variants=[self.text.replace(block,'',1),
            synthetic.replace_record(self.text,'MATRIX_BLOCK',{'rank':'1'}),
            synthetic.replace_record(self.text,'MATRIX_BLOCK',{'px':'7'}),
            synthetic.replace_record(self.text,'MATRIX_BLOCK',{'L_RE':'1'}),
            synthetic.replace_record(self.text,'SU3_TIME',{'unit_max':'2e-12'}),
            synthetic.replace_record(self.text,'SU3_WINDOW',{'t1':'121'}),
            synthetic.replace_record(self.text,'PERIODIC_ROUNDTRIP',{'status':'FAIL'}),
            synthetic.replace_record(self.text,'MATRIX_RT',{'L_IM':'1e-5'},lambda r:r['R']=='32')]
        for text in variants:
            with self.assertRaises(ValueError):self.replay(text)

    def test_equal_time_gram_accepts_singular_psd_and_rejects_invalid_matrices(self):
        for matrix in ([[complex(1) for _ in range(4)] for _ in range(4)],[[0j]*4 for _ in range(4)]):
            check.equal_time_matrix(matrix)
        negative=[[complex(int(i==j)) for j in range(4)] for i in range(4)];negative[2][2]=-1+0j
        nonhermitian=[[complex(int(i==j)) for j in range(4)] for i in range(4)];nonhermitian[1][0]=.2j
        zero_pivot=[[0j]*4 for _ in range(4)];zero_pivot[0][1]=zero_pivot[1][0]=1+0j
        for matrix in (negative,nonhermitian,zero_pivot):
            with self.assertRaises(ValueError):check.equal_time_matrix(matrix)

    def test_audit_keeps_strict_thresholds_saved_order_and_all33_times(self):
        text=synthetic.audit_fixture(640)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'audit.out';path.write_text(text)
            self.assertEqual(check.audit(path,640)['modes'],3300)
            first=next(x for x in text.splitlines(True) if x.startswith('MODE ') and 't=122 ' in x)
            variants=[text.replace('residual_relative=5e-12','residual_relative=1.01e-9',1),
                text.replace(first,'',1),text.replace('lambda=0.01 norm=1','lambda=0.5 norm=1',1),
                text.replace('gram_offdiag=1e-15','gram_offdiag=1.01e-9',1)]
            for changed in variants:
                path.write_text(changed)
                with self.assertRaises(ValueError):check.audit(path,640)


class SourceTests(unittest.TestCase):
    def test_dependency_pins_match_published_parent_manifest(self):
        current=json.loads((c.HERE/'INPUTS.json').read_text())['files']
        parent_name='laplace_static/validation/qx_static_extension/INPUTS.json'
        self.assertEqual(c.digest(c.ROOT/parent_name),current[parent_name])
        published=json.loads((c.ROOT/parent_name).read_text())['files']
        shared=set(current)&set(published)
        self.assertGreaterEqual(len(shared),35)
        for name in sorted(shared):
            self.assertEqual(current[name],published[name],name)

    def test_python36_shell_and_source_pins(self):
        for path in c.HERE.rglob('*.py'):
            if 'build' in path.parts or 'campaigns' in path.parts:continue
            options=dict(feature_version=(3,6)) if sys.version_info>=(3,8) else {}
            ast.parse(path.read_text(),filename=str(path),**options)
        for name in ('build.sh','submit_micro.sh'):
            subprocess.run(['bash','-n',str(c.HERE/name)],check=True)
        c.verify_sources()
        settings=json.loads((c.HERE/'INPUTS.json').read_text())
        self.assertEqual(settings['production_target_configurations'],list(range(100,801,10)))
        self.assertEqual(settings['pilot_configurations'],list(PILOT))

    def test_all_published_pilot_diagonals_exist(self):
        refs=check.load_references(c.REFERENCE,PILOT)
        self.assertEqual(len(refs),6336)
        self.assertEqual({k[0] for k in refs},set(PILOT))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'refs.csv';synthetic.write_references(path,[650])
            self.assertEqual(len(check.load_references(path,[650])),2112)
            with path.open('a') as f:f.write(path.read_text().splitlines()[1]+'\n')
            with self.assertRaises(ValueError):check.load_references(path,[650])

    @unittest.skipUnless(shutil.which('gcc'),'gcc is required for native synthetic tests')
    def test_native_all_rectangles_and_existing_periodic_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp)
            tests=c.HERE.parent/'qx_static_extension/tests'
            for modes in (4,100):
                exe=Path(tmp)/('mixed%d'%modes)
                subprocess.run(['gcc','-std=c11','-O2','-DNV='+str(modes),'-I'+str(tests),'-I'+tmp,
                    '-I'+str(c.HERE),str(c.HERE/'tests/mixed_contractions.c'),'-lm','-o',str(exe)],check=True)
                subprocess.run([str(exe)],check=True)
            exe=Path(tmp)/'periodic'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(tests),'-I'+tmp,
                str(tests/'kernels_and_shifts.c'),'-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is required for captured-interface syntax checks')
    def test_driver_syntax_and_unchanged_io_payload_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp)
            include=[c.HERE.parent/'qx_static_extension/tests',Path(tmp),c.ROOT/'solver/include',
                c.ROOT/'laplace_static/validation/qx_n640_diagnostic/inputs/qcd_include',
                c.ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230',c.HERE]
            subprocess.run(['gcc','-std=c11','-fsyntax-only','-Werror=implicit-function-declaration',
                '-Werror=format']+['-I'+str(p) for p in include]+[str(c.HERE/'scan.c')],check=True)
            audit=(Path(tmp)/'audit_selected.c').read_text()
            self.assertIn('static const int selected_times[]={'+','.join(map(str,AUDIT_TIMES))+'};',audit)
            self.assertNotIn('1.1e-9',audit)
            with tempfile.TemporaryDirectory() as old:
                helper.load_parent().generate(old)
                for name in ('io_helpers.h','mpi_helpers.h','audit_selected.c','metrics.h'):
                    self.assertEqual((Path(tmp)/name).read_bytes(),(Path(old)/name).read_bytes())


class CampaignTests(unittest.TestCase):
    def isolated(self):
        stack=contextlib.ExitStack();self.addCleanup(stack.close)
        root=Path(stack.enter_context(tempfile.TemporaryDirectory()))
        here=root/'laplace_static/validation/qx_static_matrix';here.mkdir(parents=True)
        (here/'README.md').write_text('Synthetic lifecycle test only.\n')
        build=here/'build';build.mkdir();(build/'MANIFEST.json').write_text('{}\n')
        ref=root/'reference.csv';synthetic.write_references(ref)
        physical=root/'physical';physical.mkdir()
        inputs={}
        for cfg in PILOT:
            folder=physical/str(cfg);folder.mkdir();ev=folder/'vectors';ev.mkdir()
            gauge=folder/('qX_p2n%d'%cfg);gauge.write_bytes(b'fixture gauge')
            eigenvalues=folder/'values.h5';eigenvalues.write_bytes(b'fixture values')
            vectors=[]
            for t in AUDIT_TIMES:
                p=ev/('qX_p2n%d_evec_t%d.h5'%(cfg,t));p.write_bytes(b'fixture vector');vectors.append(c.fingerprint(p))
            inputs[cfg]=dict(configuration=cfg,eigenvectors=str(ev),eigenvalues=str(eigenvalues),
                eigenvalue_sha256=c.digest(eigenvalues),gauge=c.fingerprint(gauge),vectors=vectors)
        stack.enter_context(mock.patch.multiple(c,ROOT=root,HERE=here,BUILD=build,REFERENCE=ref))
        stack.enter_context(mock.patch.object(c,'verify_sources'))
        stack.enter_context(mock.patch.object(c,'build_manifest',return_value='synthetic-build-hash'))
        stack.enter_context(mock.patch.object(c,'input_record',side_effect=lambda cfg:inputs[cfg]))
        stack.enter_context(mock.patch.dict(os.environ,dict(USER='fixture',SLURM_NTASKS='64',SLURM_JOB_NUM_NODES='8',SLURM_JOB_ID='123456',SLURM_ARRAY_TASK_ID='0')))
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        args=argparse.Namespace(name='pilot',concurrent=3)
        c.prepare(args)
        return root,here,args,inputs

    def fake_mpi(self,command,**kwargs):
        cfg=int(kwargs['env']['QX_CFG']);piter=int(kwargs['env']['HYP_PITER'])
        text=synthetic.audit_fixture(cfg) if Path(command[3]).name=='audit_selected' else synthetic.scan_fixture(cfg,piter)
        kwargs['stdout'].write(text)
        return subprocess.CompletedProcess(command,0)

    def test_three_case_lifecycle_failure_retry_partial_export_and_tamper(self):
        root,here,args,inputs=self.isolated();target=here/'campaigns/pilot'
        with self.assertRaises(SystemExit) as stopped:c.collect(args)
        self.assertEqual(stopped.exception.code,1)
        self.assertFalse(json.loads((target/'SUMMARY.json').read_text())['success'])
        with mock.patch.object(c.subprocess,'run',side_effect=self.fake_mpi):
            for index in range(3):
                os.environ['SLURM_ARRAY_TASK_ID']=str(index);c.run(args)
            c.run(args)  # The valid marker replays rather than launching again.
        c.collect(args);summary=json.loads((target/'SUMMARY.json').read_text())
        self.assertTrue(summary['success'])
        self.assertEqual((summary['matrix_rows'],summary['rank_block_rows']),(33792,49152))
        archive=root/'qx_static_matrix_pilot_review.zip'
        with zipfile.ZipFile(archive) as z:
            self.assertIsNone(z.testzip());self.assertTrue(any(n.endswith('/p40.err') for n in z.namelist()))
            self.assertFalse(any('/physical/' in n for n in z.namelist()))
        entry=json.loads((target/'PREPARED.json').read_text())['cases'][2]
        out=c.complete_case(target,entry)
        path=out/'matrix_correlators.csv';saved=path.read_bytes();path.write_bytes(saved+b'changed')
        with self.assertRaises(ValueError):c.collect(args)
        path.write_bytes(saved);(out/'p40.err').write_text('changed')
        with self.assertRaises(ValueError):c.collect(args)

    def test_failed_launch_has_no_success_marker_and_retry_preserves_attempt(self):
        _,here,args,_=self.isolated();target=here/'campaigns/pilot';os.environ['SLURM_ARRAY_TASK_ID']='2'
        def failed(command,**kwargs):
            if kwargs['env']['HYP_PITER']=='40':return subprocess.CompletedProcess(command,2)
            return self.fake_mpi(command,**kwargs)
        with mock.patch.object(c.subprocess,'run',side_effect=failed):
            with self.assertRaises(ValueError):c.run(args)
        case=target/'case_002'
        self.assertFalse((case/'LATEST.json').exists());self.assertTrue((case/'attempt_0001/FAILED.json').exists())
        with mock.patch.object(c.subprocess,'run',side_effect=self.fake_mpi):c.run(args)
        self.assertEqual(json.loads((case/'LATEST.json').read_text())['attempt'],'attempt_0002')
        self.assertTrue((case/'attempt_0001/FAILED.json').exists())

    def test_submission_blocks_active_ambiguous_or_changed_physical_inputs(self):
        _,here,args,inputs=self.isolated();target=here/'campaigns/pilot'
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'123456\n','')) as run:
            c.submit(args);self.assertEqual(run.call_count,1)
            self.assertIn('--array=0,1,2%3',run.call_args[0][0])
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'123456 123456\n','')) as run:
            with self.assertRaises(ValueError):c.submit(args)
            self.assertEqual(run.call_count,1);self.assertEqual(run.call_args[0][0][0],'squeue')
        c.write(target/'SUBMISSION_002.json',dict(state='submitting'))
        with self.assertRaises(ValueError):c.submit(args)
        Path(inputs[170]['vectors'][0]['path']).write_bytes(b'changed physical input')
        with self.assertRaises(ValueError):c.prepared('pilot',physical=True)

    def test_plan_checksum_geometry_and_build_binding(self):
        _,here,_,_=self.isolated();path=here/'campaigns/pilot/PREPARED.json'
        plan=json.loads(path.read_text());plan['geometry']['Nv_sink']=[40,100];c.write(path,plan)
        with self.assertRaises(ValueError):c.prepared('pilot')
        (path.parent/'prepared.sha256').write_text(c.digest(path)+'\n')
        with self.assertRaises(ValueError):c.prepared('pilot')


if __name__=='__main__':unittest.main()
