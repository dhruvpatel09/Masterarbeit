#!/usr/bin/env python3
"""Full71 orchestration tests; MPI and physical storage are mocked explicitly."""
import argparse
import ast
import csv
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
import zipfile
from pathlib import Path
from unittest import mock

import campaign as c


class PublishedSourceTests(unittest.TestCase):
    def test_exact_source_and_reference_pins(self):
        c.verify_sources()
        base=c.load_base()
        base.verify_sources()
        self.assertEqual(base.geometry(),c.geometry())
        self.assertEqual(len(base.check.expected_keys()),8448)
        self.assertEqual(len(base.check.block_keys())*64,12288)

    def test_base_import_restores_generic_module_names(self):
        sentinel=types.ModuleType('check_outputs')
        original=c.BASE;c.BASE=None
        try:
            with mock.patch.dict(sys.modules,{'check_outputs':sentinel}):
                base=c.load_base()
                self.assertIs(sys.modules['check_outputs'],sentinel)
                self.assertIsNot(base.check,sentinel)
        finally:c.BASE=original

    def test_python36_and_shell_syntax(self):
        for p in c.HERE.glob('*.py'):
            source=p.read_text()
            if sys.version_info>=(3,8):ast.parse(source,filename=str(p),feature_version=(3,6))
            else:ast.parse(source,filename=str(p))
        subprocess.run(['bash','-n',str(c.HERE/'submit_micro.sh')],check=True)

    def test_every_target_is_counted_once_and_slots_do_not_include_reuse(self):
        indices=[i for slot in c.slots() for i in slot]
        self.assertEqual(indices,list(range(68)))
        self.assertEqual(len(c.slots()),34)
        self.assertTrue(all(len(s)==2 for s in c.slots()))
        self.assertEqual(sorted(c.NEW+c.PILOT),c.TARGET)
        self.assertFalse(set(c.NEW)&set(c.PILOT))

    def test_actual_published_pilot_stdout_and_csv_replay_without_physical_storage(self):
        base=c.load_base();source=c.read(c.HERE/'INPUTS.json')
        archive=c.ROOT/source['published_pilot_archive']
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            with zipfile.ZipFile(archive) as z:
                names=z.namelist()
                self.assertEqual(len(names),len(set(names)))
                for item in z.infolist():
                    path=Path(item.filename)
                    if path.is_absolute() or '..' in path.parts or (item.external_attr>>16)&0o170000==0o120000:
                        raise ValueError('Unsafe archived pilot path')
                z.extractall(root)
            target=root/'laplace_static/validation/qx_static_matrix/campaigns/pilot'
            pp=c.read(target/'PREPARED.json')
            self.assertEqual(c.digest(target/'PREPARED.json'),(target/'prepared.sha256').read_text().strip())
            build=root/'laplace_static/validation/qx_static_matrix/build/MANIFEST.json'
            self.assertEqual(c.digest(build),pp['build_manifest_sha256'])
            def offline_prepared(name,physical=False):
                if name!='pilot' or physical:raise ValueError('Offline replay has no physical storage')
                return target,pp
            with mock.patch.object(c,'ROOT',root),mock.patch.object(c,'load_base',return_value=base):
                with mock.patch.object(base,'prepared',side_effect=offline_prepared):
                    _,_,receipts=c.admit_pilot('pilot',reparse=True)
            self.assertEqual(receipts,source['pilot_receipts']['cases'])


class LifecycleTests(unittest.TestCase):
    """Reuse genuine case-marker/CSV helpers; replace only HPC-specific actions."""
    @classmethod
    def setUpClass(cls):
        cls.native=c.load_base()

    def setUp(self):
        silent=mock.patch('sys.stdout',new=io.StringIO())
        silent.start();self.addCleanup(silent.stop)
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.here=self.root/'full71';self.here.mkdir()
        self.pilot=self.root/'native/campaigns/pilot';self.pilot.mkdir(parents=True)
        self.build=self.root/'native/build';self.build.mkdir()
        (self.build/'MANIFEST.json').write_text(json.dumps(dict(files={})))
        self.reference=self.root/'reference.csv';self.reference.write_text('published diagonal reference\n')
        (self.here/'INPUTS.json').write_text(json.dumps(dict(files={})))
        (self.here/'FILE_MANIFEST.json').write_text(json.dumps(dict(files={},dependencies={})))
        self.build_hash=c.digest(self.build/'MANIFEST.json')
        self.starts=[];self.physical=[]
        self.fake=types.SimpleNamespace(HERE=self.root/'native',BUILD=self.build,REFERENCE=self.reference,
            complete_case=self.native.complete_case,csv_write=self.native.csv_write,
            build_manifest=lambda:self.build_hash,prepared=self.pilot_prepared,
            input_record=lambda cfg:dict(configuration=cfg,admitted=True),
            validate_inputs=self.validate_inputs,replay_case=self.replay_case,
            launch=self.launch,live_previous=self.native.live_previous,
            check=types.SimpleNamespace(load_references=self.references,audit=self.audit))
        for name,value in [('ROOT',self.root),('HERE',self.here),('MATRIX_ROWS',1),('BLOCK_ROWS',1)]:
            patch=mock.patch.object(c,name,value);patch.start();self.addCleanup(patch.stop)
        patch=mock.patch.object(c,'load_base',return_value=self.fake);patch.start();self.addCleanup(patch.stop)
        cases=[dict(index=i,configuration=cfg,piters=[10,40] if cfg==650 else [10],
                    input=dict(configuration=cfg,admitted=True)) for i,cfg in enumerate(c.PILOT)]
        self.pp=dict(stage='technical_pilot',geometry=c.geometry(),cases=cases,
                     build_manifest_sha256=self.build_hash)
        (self.pilot/'PREPARED.json').write_text(json.dumps(self.pp))
        (self.pilot/'prepared.sha256').write_text(c.digest(self.pilot/'PREPARED.json')+'\n')
        for entry in cases:self.seed_case(self.pilot,entry)
        receipts=[]
        for entry in cases:
            out=self.native.complete_case(self.pilot,entry)
            receipts.append(dict(configuration=entry['configuration'],attempt=out.relative_to(self.root).as_posix(),
                summary_sha256=c.digest(out/'SUMMARY.json'),case_sha256=c.digest(out/'CASE.json'),
                launches_sha256=c.digest(out/'LAUNCHES.json')))
        self.published_pilot=self.root/'published_pilot.zip';self.published_pilot.write_text('published pilot bytes')
        (self.here/'INPUTS.json').write_text(json.dumps(dict(files={},published_pilot_archive='published_pilot.zip',
            pilot_receipts=dict(name='pilot',prepared_sha256=c.digest(self.pilot/'PREPARED.json'),cases=receipts))))
        self.starts.clear()

    def references(self,path,configs):
        if Path(path)!=self.reference or not set(configs).issubset(c.TARGET):
            raise ValueError('Reference admission')
        return set(configs)

    def validate_inputs(self,entry):
        if not entry.get('admitted') or entry['configuration'] not in c.TARGET:
            raise ValueError('Physical input changed')
        self.physical.append(entry['configuration'])

    def pilot_prepared(self,name,physical=False):
        if name!='pilot':raise ValueError('Unknown pilot')
        if c.digest(self.pilot/'PREPARED.json')!=(self.pilot/'prepared.sha256').read_text().strip():
            raise ValueError('Pilot plan checksum')
        pp=c.read(self.pilot/'PREPARED.json')
        if physical:
            for e in pp['cases']:self.validate_inputs(e['input'])
        return self.pilot,pp

    def audit(self,path,cfg):
        if c.read(path)!=dict(configuration=cfg,kind='audit'):raise ValueError('Audit failed')
        return dict(configuration=cfg,modes=3300)

    def launch(self,out,kind,piter,entry,launches):
        cfg=entry['configuration'];stem='audit' if kind=='audit' else 'p%d'%piter
        (out/(stem+'.out')).write_text(json.dumps(dict(configuration=cfg,kind=kind)))
        (out/(stem+'.err')).write_text('')
        launches.append(dict(kind=kind,piter=piter,configuration=cfg,stem=stem,seconds=1.,exit_code=0,
            out_sha256=c.digest(out/(stem+'.out')),err_sha256=c.digest(out/(stem+'.err'))))
        c.write(out/'LAUNCHES.json',launches)
        self.starts.append((cfg,kind,piter))

    def replay_case(self,out,entry,references):
        cfg=entry['configuration'];launches=c.read(out/'LAUNCHES.json')
        if cfg not in references:raise ValueError('Missing diagonal reference')
        expected=[('audit',0)]+[('scan',p) for p in entry['piters']]
        if [(r['kind'],r['piter']) for r in launches]!=expected:raise ValueError('Launch inventory')
        for r in launches:
            if r['exit_code']!=0:raise ValueError('Failed native launch')
            for suffix in ('out','err'):
                if c.digest(out/(r['stem']+'.'+suffix))!=r[suffix+'_sha256']:
                    raise ValueError('Native output changed')
        self.audit(out/'audit.out',cfg)
        rows=[dict(cfg=str(cfg),piter=str(p),Nv_src='40',Nv_sink='40',R='4',tau='0',
                   axis='-1',L_RE=str(cfg),L_IM='0') for p in entry['piters']]
        blocks=[dict(cfg=str(cfg),piter=str(p),rank='0',L_RE=str(cfg)) for p in entry['piters']]
        return rows,blocks,dict(configuration=cfg,audit=dict(modes=3300),scans=list(entry['piters']))

    def seed_case(self,target,entry):
        case=target/('case_%03d'%entry['index']);case.mkdir()
        out=case/'attempt_0001';out.mkdir();launches=[]
        c.write(out/'CASE.json',dict(entry,prepared_sha256=c.digest(target/'PREPARED.json')))
        self.launch(out,'audit',0,entry,launches)
        for p in entry['piters']:self.launch(out,'scan',p,entry,launches)
        rows,blocks,summary=self.replay_case(out,entry,set(c.TARGET))
        self.native.csv_write(out/'matrix_correlators.csv',rows);self.native.csv_write(out/'rank_blocks.csv',blocks)
        c.write(out/'SUMMARY.json',dict(summary,success=True,launches_sha256=c.digest(out/'LAUNCHES.json'),
            case_sha256=c.digest(out/'CASE.json'),
            outputs_sha256={n:c.digest(out/n) for n in ('matrix_correlators.csv','rank_blocks.csv')}))
        c.write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=c.digest(out/'SUMMARY.json')))
        return out

    def prepare(self):
        c.prepare(argparse.Namespace(name='full71',pilot='pilot'))
        return c.prepared('full71')

    def test_prepare_and_replay_reuses_three_and_preserves_pilot_controls(self):
        target,plan=self.prepare()
        self.assertEqual(len(plan['cases']),68);self.assertEqual(len(plan['slots']),34)
        self.assertEqual(plan['reused_configurations'],[170,640,650])
        self.assertEqual(self.starts,[])
        self.assertTrue(all(e['piters']==[10] for e in plan['cases']))

    def test_pilot_receipt_tamper_is_rejected(self):
        target,plan=self.prepare()
        path=self.root/plan['pilot']['cases'][0]['attempt']/'SUMMARY.json'
        path.write_text(path.read_text()+' ')
        with self.assertRaises(ValueError):c.prepared('full71')

    def test_pilot_stdout_change_is_rejected_before_submission(self):
        target,plan=self.prepare()
        path=self.root/plan['pilot']['cases'][0]['attempt']/'p10.out'
        path.write_text(path.read_text()+' ')
        with self.assertRaises(ValueError):c.prepared('full71')

    def test_pilot_geometry_changes_are_not_reused(self):
        pp=c.read(self.pilot/'PREPARED.json');pp['geometry']['R']=[4,12,23,24,32]
        (self.pilot/'PREPARED.json').write_text(json.dumps(pp))
        (self.pilot/'prepared.sha256').write_text(c.digest(self.pilot/'PREPARED.json'))
        with self.assertRaises(ValueError):c.prepare(argparse.Namespace(name='full71',pilot='pilot'))

    def test_case_assembly_corruption_is_rejected_even_with_updated_csv_receipt(self):
        e=self.pp['cases'][0];out=self.native.complete_case(self.pilot,e)
        rows=c.csv_rows(out/'matrix_correlators.csv');rows[0]['L_RE']='wrong'
        self.native.csv_write(out/'matrix_correlators.csv',rows)
        summary=c.read(out/'SUMMARY.json');summary['outputs_sha256']['matrix_correlators.csv']=c.digest(out/'matrix_correlators.csv')
        c.write(out/'SUMMARY.json',summary)
        c.write(out.parent/'LATEST.json',dict(attempt=out.name,summary_sha256=c.digest(out/'SUMMARY.json')))
        with self.assertRaises(ValueError):c.replay(out,e,set(c.TARGET))

    def test_changed_slot_plan_and_workflow_are_rejected(self):
        target,plan=self.prepare()
        plan['slots'][0]=[0,0]
        (target/'PREPARED.json').write_text(json.dumps(plan))
        (target/'prepared.sha256').write_text(c.digest(target/'PREPARED.json'))
        with self.assertRaises(ValueError):c.prepared('full71')

    def test_changed_workflow_manifest_and_native_build_are_rejected(self):
        target,plan=self.prepare()
        original=(self.here/'FILE_MANIFEST.json').read_text()
        (self.here/'FILE_MANIFEST.json').write_text(original+' ')
        with self.assertRaises(ValueError):c.prepared('full71')
        (self.here/'FILE_MANIFEST.json').write_text(original)
        self.fake.build_manifest=lambda:'different native binary manifest'
        with self.assertRaises(ValueError):c.prepared('full71')

    def test_new_case_cannot_claim_another_configurations_inputs(self):
        target,plan=self.prepare()
        plan['cases'][0]['input']['configuration']=110
        (target/'PREPARED.json').write_text(json.dumps(plan))
        (target/'prepared.sha256').write_text(c.digest(target/'PREPARED.json'))
        with self.assertRaises(ValueError):c.prepared('full71')

    def test_retry_skips_validated_case_and_preserves_failed_attempt(self):
        target,plan=self.prepare();entry=plan['cases'][0]
        original=self.fake.launch
        def broken(out,kind,piter,entry,launches):
            original(out,kind,piter,entry,launches)
            if kind=='scan':raise ValueError('simulated MPI failure')
        self.fake.launch=broken
        with self.assertRaises(ValueError):c.run_case(target,entry,set(c.TARGET))
        case=target/'case_000'
        self.assertTrue((case/'attempt_0001/FAILED.json').is_file())
        self.assertFalse((case/'LATEST.json').exists());self.assertFalse((case/'RUNNING.json').exists())
        self.fake.launch=original;c.run_case(target,entry,set(c.TARGET))
        self.assertEqual(c.read(case/'LATEST.json')['attempt'],'attempt_0002')
        starts=len(self.starts);c.run_case(target,entry,set(c.TARGET))
        self.assertEqual(len(self.starts),starts)

    def test_collect_requires_the_full_cohort_before_writing_tables(self):
        target,plan=self.prepare()
        with self.assertRaises(ValueError):c.collect(argparse.Namespace(name='full71'))
        self.assertFalse((target/'SUMMARY.json').exists())
        self.assertFalse((target/'matrix_correlators.csv').exists())

    def test_full71_streaming_collect_excludes_p40_and_archives_pilot(self):
        target,plan=self.prepare()
        for e in plan['cases']:c.run_case(target,e,set(c.TARGET))
        c.collect(argparse.Namespace(name='full71'))
        rows=c.csv_rows(target/'matrix_correlators.csv')
        self.assertEqual([int(r['cfg']) for r in rows],c.TARGET)
        self.assertEqual(set(r['piter'] for r in rows),{'10'})
        summary=c.read(target/'SUMMARY.json')
        self.assertEqual(summary['matrix_rows'],71)
        self.assertEqual(sum(r['reused_from_pilot'] for r in summary['cases']),3)
        archive=self.root/'qx_static_matrix_full71_review.zip'
        with zipfile.ZipFile(archive) as z:
            self.assertIsNone(z.testzip())
            self.assertIn('native/campaigns/pilot/case_002/attempt_0001/p40.out',z.namelist())
            self.assertIn('full71/campaigns/full71/ARCHIVE_FILES.json',z.namelist())

    def test_collect_failure_cleans_partial_streams_without_success_summary(self):
        target,plan=self.prepare()
        for e in plan['cases']:self.seed_case(target,e)
        entry=plan['cases'][9]
        out=self.fake.complete_case(target,entry)
        (out/'p10.out').write_text('changed stdout')
        with self.assertRaises(ValueError):c.collect(argparse.Namespace(name='full71'))
        self.assertFalse((target/'matrix_correlators.csv.pending').exists())
        self.assertFalse((target/'rank_blocks.csv.pending').exists())
        self.assertFalse((target/'SUMMARY.json').exists())
        self.assertFalse((target/'matrix_correlators.csv').exists())

    def test_live_and_ambiguous_submission_receipts_block_duplicate_submission(self):
        target,plan=self.prepare()
        c.write(target/'SUBMISSION_001.json',dict(returncode=0,stdout='12345\n'))
        with mock.patch.dict(os.environ,{'USER':'test'}):
            with mock.patch('subprocess.run',return_value=subprocess.CompletedProcess([],0,'12345 12345\n','')):
                with self.assertRaises(ValueError):self.fake.live_previous(target)
        c.write(target/'SUBMISSION_001.json',dict(state='submitting'))
        with self.assertRaises(ValueError):self.fake.live_previous(target)

    def test_submit_array_has_34_slots_and_four_task_cap(self):
        target,plan=self.prepare()
        self.fake.live_previous=lambda p:None
        with mock.patch('subprocess.run',return_value=subprocess.CompletedProcess([],0,'98765\n','')) as run:
            c.submit(argparse.Namespace(name='full71',start=0,count=None,concurrent=4))
        cmd=run.call_args[0][0]
        self.assertIn('--array='+','.join(map(str,range(34)))+'%4',cmd)
        self.assertEqual(cmd[-1],'full71')

    def test_submit_skips_a_completed_slot_and_blocks_running_case(self):
        target,plan=self.prepare();self.fake.live_previous=lambda p:None
        for e in plan['cases'][:2]:self.seed_case(target,e)
        with mock.patch('subprocess.run',return_value=subprocess.CompletedProcess([],0,'98765\n','')) as run:
            c.submit(argparse.Namespace(name='full71',start=0,count=2,concurrent=4))
        self.assertIn('--array=1%4',run.call_args[0][0])
        case=target/'case_002';case.mkdir();(case/'RUNNING.json').write_text('{}')
        with self.assertRaises(ValueError):c.submit(argparse.Namespace(name='full71',start=1,count=1,concurrent=4))

    def test_slot_ranges_concurrency_and_allocation_are_checked(self):
        target,plan=self.prepare()
        for start,count in [(-1,1),(34,1),(0,35),(0,0)]:
            with self.assertRaises(ValueError):c.selected_slots(argparse.Namespace(start=start,count=count),plan)
        with self.assertRaises(ValueError):c.submit(argparse.Namespace(name='full71',start=0,count=None,concurrent=5))
        with mock.patch.dict(os.environ,{'SLURM_NTASKS':'1','SLURM_JOB_NUM_NODES':'1','SLURM_ARRAY_TASK_ID':'0'}):
            with self.assertRaises(ValueError):c.run(argparse.Namespace(name='full71'))

    def test_one_array_worker_runs_exactly_two_new_configurations(self):
        target,plan=self.prepare()
        with mock.patch.dict(os.environ,{'SLURM_NTASKS':'64','SLURM_JOB_NUM_NODES':'8','SLURM_ARRAY_TASK_ID':'0'}):
            c.run(argparse.Namespace(name='full71'))
        self.assertEqual(self.starts,[(100,'audit',0),(100,'scan',10),(110,'audit',0),(110,'scan',10)])
        self.assertIsNotNone(self.fake.complete_case(target,plan['cases'][0]))
        self.assertIsNotNone(self.fake.complete_case(target,plan['cases'][1]))

    def test_case_lock_is_exclusive(self):
        case=self.root/'lock_case';case.mkdir()
        with c.case_lock(case):
            with self.assertRaises(ValueError):
                with c.case_lock(case):pass
        self.assertFalse((case/'RUNNING.json').exists())


if __name__=='__main__':unittest.main()
