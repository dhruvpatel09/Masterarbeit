import contextlib
import csv
import io
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

BASE=Path(__file__).resolve().parents[3]
HELPER=BASE/'laplace_static/validation/hyp_projection_configs_qx_200_230_queue'
sys.path.insert(0,str(HELPER))
import queue_common as q
import run_remaining as r
import collect_mixed as collector


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.root=Path(self.temporary.name)
        case=self.root/'laplace_static/validation/hyp_projection_configs_qx_200_230'
        helper=case.with_name('hyp_projection_configs_qx_200_230_queue')
        shutil.copytree(BASE/'laplace_static/validation/hyp_projection_configs_qx_200_230',case)
        shutil.copytree(HELPER,helper,ignore=shutil.ignore_patterns('__pycache__'))
        (case/'build').mkdir()
        (case/'build/manifest.sha256').write_text('fixture build provenance\n')
        (case/'prepared.sha256').write_text('fixture preparation provenance\n')
        self.patches=[patch.object(q.c,'ROOT',self.root),patch.object(q.c,'CASE',case),
                      patch.object(q,'HERE',helper),patch.object(q,'LEDGER',helper/'RUNS.json'),
                      patch.object(q.c,'prepared',return_value={}),
                      patch.object(q.checks,'check',side_effect=self.fake_scientific_check),
                      patch.object(r.time,'sleep',return_value=None)]
        for p in self.patches:p.start()
        self.checks=[];self.submissions=[];self.active=None;self.stage=0
        self.create_outputs(0,'5918419')
        self.silent=contextlib.redirect_stdout(io.StringIO())
        self.silent.__enter__()

    def tearDown(self):
        self.silent.__exit__(None,None,None)
        for p in reversed(self.patches):p.stop()
        self.temporary.cleanup()

    def create_outputs(self,index,job):
        selected=q.c.case_for(index)
        folder=q.c.CASE/f'results_{job}_{index}'
        folder.mkdir()
        (folder/'CASE.json').write_text(json.dumps(dict(selected,array_job_id=job)))
        for name in q.PROVENANCE:
            (folder/name.replace('/','_')).write_bytes((q.c.CASE/name).read_bytes())
        launches=[dict(configuration=selected['configuration'],tau=selected['tau'],
                       piter=p,label=f'p{p}',exit_code=0,seconds=5) for p in (4,40)]
        if index==0:
            launches += [dict(configuration=100,tau=4,piter=p,label=f'n100_p{p}',
                              exit_code=0,seconds=5) for p in (4,40)]
        (folder/'LAUNCHES.json').write_text(json.dumps(launches))
        summary=dict(configuration=selected['configuration'],tau=selected['tau'],success=True,
                     p4_window=dict(unit_max='1e-6',det_max='5e-7'),
                     max_relative_difference=dict(L_RE=1e-8,C_combined=1e-5))
        if index==0:summary['n100_regression']=dict(success=True)
        (folder/'SUMMARY.json').write_text(json.dumps(summary))
        for name,columns,count in [('profiles_raw.csv',q.checks.RAW_COLUMNS,664),
                                   ('profiles_paired.csv',q.checks.PAIRED_COLUMNS,332),
                                   ('su3_times.csv',q.checks.SU3_COLUMNS,2*selected['tau']),
                                   ('su3_windows.csv',q.checks.SU3_COLUMNS,2)]:
            row={k:'0' for k in columns}
            row.update(configuration=str(selected['configuration']),tau=str(selected['tau']))
            q.checks.write_csv(folder/name,columns,[row.copy() for _ in range(count)])
        return folder

    def fake_scientific_check(self,folder,cfg,tau,regression=False):
        # These tests exercise new admission/provenance/coverage logic only.
        # Original numerical checkers are unchanged and run on SuperMUC logs.
        s=json.loads((folder/'SUMMARY.json').read_text())
        self.assertEqual((s['configuration'],s['tau']),(cfg,tau))
        if regression:self.assertEqual(cfg,200)
        self.checks.append(str(folder))
        return s

    def fake_slurm(self,command,**kwargs):
        if command[0]=='sbatch':
            self.assertIsNone(self.active,'Two cases admitted at once')
            self.assertEqual(command[1],'--parsable')
            index=int(command[2].split('=')[1])
            self.assertTrue(1<=index<=11)
            job=str(710000+index)
            self.create_outputs(index,job)
            self.active=(index,job);self.stage=0;self.submissions.append(index)
            return subprocess.CompletedProcess(command,0,job+'\n','')
        self.assertEqual(command[0],'squeue')
        self.assertIn('--array',command)
        self.assertNotIn('--jobs',command)
        if self.active is None:return subprocess.CompletedProcess(command,0,'','')
        index,job=self.active
        if self.stage<2:
            state=('PENDING','RUNNING')[self.stage]
            self.stage+=1
            return subprocess.CompletedProcess(command,0,f'{job}_{index}|{state}\n','')
        self.active=None
        return subprocess.CompletedProcess(command,0,'','')

    def test_sequential_campaign_collect_and_resume(self):
        q.verify_environment()
        data=r.initialize('5918419')
        with patch.object(r.subprocess,'run',side_effect=self.fake_slurm):
            r.run(data,30)
            self.assertEqual(self.submissions,list(range(1,12)))
            saved=q.load_ledger()
            self.assertTrue(all(e['state']=='completed' for e in saved['cases']))
            archive=Path(saved['review_archive'])
            self.assertTrue(archive.is_file())
            collection=q.HERE/'collection_mixed_5918419'
            summary=json.loads((collection/'SUMMARY.json').read_text())
            self.assertTrue(summary['success'])
            self.assertEqual(len(set(summary['job_mapping'].values())),12)
            for name,count in [('profiles_raw.csv',7968),('profiles_paired.csv',3984),
                               ('su3_times.csv',96),('su3_windows.csv',24)]:
                with (collection/name).open() as f:self.assertEqual(len(list(csv.DictReader(f))),count)
            before=q.c.digest(archive)
            r.run(saved,30)
            self.assertEqual(self.submissions,list(range(1,12)))
            self.assertEqual(q.c.digest(archive),before)

    def test_initial_job_can_still_be_pending_then_running(self):
        shutil.rmtree(q.c.CASE/'results_5918419_0')
        data=r.initialize('5918419')
        self.assertEqual(data['cases'][0]['state'],'submitted')
        states=iter(['PENDING','RUNNING',None])
        def initial_queue(entry):
            state=next(states)
            if state is None:self.create_outputs(0,'5918419')
            return state
        with patch.object(r,'active_job',side_effect=initial_queue):
            result=r.wait_for_completion(data['cases'][0],30)
        self.assertTrue(result['success'])
        self.assertEqual(self.submissions,[])

    def test_failed_initial_mpi_stops_before_remaining_submissions(self):
        data=r.initialize('5918419')
        p=q.c.CASE/'results_5918419_0/LAUNCHES.json'
        launches=json.loads(p.read_text());launches[0]['exit_code']=1;p.write_text(json.dumps(launches))
        with patch.object(r.subprocess,'run',side_effect=self.fake_slurm):
            with self.assertRaisesRegex(RuntimeError,'MPI launch'):r.run(data,30)
        self.assertEqual(self.submissions,[])

    def test_submit_limit_retains_unsubmitted_case(self):
        data=r.initialize('5918419');entry=data['cases'][1]
        response=subprocess.CompletedProcess([],1,'','QOSMaxSubmitJobPerUserLimit')
        with patch.object(r.subprocess,'run',return_value=response):
            with self.assertRaisesRegex(RuntimeError,'QoS submit limit'):r.submit_one(data,entry)
        saved=q.load_ledger()['cases'][1]
        self.assertEqual(saved['state'],'planned');self.assertIsNone(saved['array_job_id'])

    def test_timeout_discards_stale_reply_and_remains_ambiguous(self):
        data=r.initialize('5918419');entry=data['cases'][1]
        entry['sbatch_reply']=dict(returncode=1,stdout='',stderr='QOSMaxSubmitJobPerUserLimit')
        with patch.object(r.subprocess,'run',side_effect=subprocess.TimeoutExpired('sbatch',60)):
            with self.assertRaises(subprocess.TimeoutExpired):r.submit_one(data,entry)
        saved=q.load_ledger();entry=saved['cases'][1]
        self.assertEqual(entry['state'],'submitting');self.assertNotIn('sbatch_reply',entry)
        with self.assertRaisesRegex(RuntimeError,'ambiguous'):r.recover_reply(saved,entry)

    def test_recover_recorded_accepted_reply_without_resubmit(self):
        data=r.initialize('5918419');entry=data['cases'][1]
        entry.update(state='submitting',sbatch_reply=dict(returncode=0,stdout='710001;local\n',stderr=''))
        with patch.object(r.subprocess,'run') as call:r.recover_reply(data,entry)
        call.assert_not_called()
        self.assertEqual(entry['array_job_id'],'710001');self.assertEqual(entry['state'],'submitted')

    def test_queue_error_is_not_completion(self):
        data=r.initialize('5918419');entry=data['cases'][1]
        entry.update(state='submitted',array_job_id='710001')
        response=subprocess.CompletedProcess([],1,'','controller not reachable')
        before=len(self.checks)
        with patch.object(r.subprocess,'run',return_value=response):
            with self.assertRaisesRegex(RuntimeError,'squeue failed'):r.wait_for_completion(entry,30)
        self.assertEqual(len(self.checks),before)

    def test_failed_mpi_stops_further_allocation(self):
        data=r.initialize('5918419')
        def failed_slurm(command,**kwargs):
            result=self.fake_slurm(command,**kwargs)
            if command[0]=='sbatch' and self.submissions[-1]==2:
                path=q.c.CASE/'results_710002_2/LAUNCHES.json'
                launches=json.loads(path.read_text());launches[0]['exit_code']=1
                path.write_text(json.dumps(launches))
            return result
        with patch.object(r.subprocess,'run',side_effect=failed_slurm):
            with self.assertRaisesRegex(RuntimeError,'MPI launch'):r.run(data,30)
        self.assertEqual(self.submissions,[1,2]);self.assertFalse((self.root/'hyp_projection_configs_200_230_review_mixed_5918419.zip').exists())

    def test_wrong_case_and_changed_provenance_rejected(self):
        data=r.initialize('5918419');first=data['cases'][0]
        path=q.output_directory(first)/'CASE.json'
        case=json.loads(path.read_text());case['configuration']=210;path.write_text(json.dumps(case))
        with self.assertRaisesRegex(RuntimeError,'Wrong case'):q.validate_case(first)
        path.write_text(json.dumps(dict(q.c.case_for(0),array_job_id='5918419')))
        (q.output_directory(first)/'package.sha256').write_text('changed provenance')
        with self.assertRaisesRegex(RuntimeError,'provenance changed'):q.validate_case(first)

    def test_missing_cases_full_refused_partial_labelled(self):
        data=r.initialize('5918419')
        with self.assertRaisesRegex(RuntimeError,'Incomplete/failed'):collector.collect(data)
        self.assertFalse((q.HERE/'collection_mixed_5918419').exists())
        archive=collector.collect(data,partial=True)
        self.assertIn('partial_',archive.name)
        folder=next(q.HERE.glob('collection_partial_*'))
        summary=json.loads((folder/'SUMMARY.json').read_text())
        self.assertFalse(summary['success']);self.assertEqual(summary['validated_cases'],1)
        self.assertEqual(len(summary['errors']),11)

    def test_duplicate_job_mapping_rejected(self):
        data=r.initialize('5918419')
        data['cases'][1].update(state='submitted',array_job_id='5918419')
        with self.assertRaisesRegex(RuntimeError,'distinct'):q.validate_ledger(data)


if __name__=='__main__':unittest.main(verbosity=2)
