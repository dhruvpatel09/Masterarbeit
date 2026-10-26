#!/usr/bin/env python3
"""Scientific regression checks for completing the matched qX cohort."""
import csv
import importlib.util
import json
import shutil
import tempfile
import unittest
from decimal import Decimal as D, localcontext
from pathlib import Path

import reproduce as a
import verify_review


class CompleteEnsembleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trace = a.module('qx71_test_trace', a.TRACE70 / 'reproduce.py')
        cls.temporal = a.module('qx71_test_temporal', a.TIME70 / 'reproduce.py')
        cls.binning = a.module('qx71_test_binning', a.BIN70 / 'reproduce.py')

    def test_ratio_of_means_and_trace_identity(self):
        values = {(c,p):tuple(map(D,row)) for c,row in
                  [(100,('1','3','1','.1','.2')),(110,('3','3','1','.1','.2'))]
                  for p in (4,40)}
        result = self.trace.paired_statistics(values,[100,110])
        self.assertEqual(result['C_RAW_p40'], D('.5'))
        self.assertNotEqual(result['C_RAW_p40'], (D(3)/1+D(3)/3)/2-1)
        self.assertEqual(result['C_RAW_p40']-result['C_TL_p40'],result['C_TRACE_p40'])

    def test_paired_errors_preserve_exact_common_fluctuations(self):
        values = {(c,p):tuple(map(D,row)) for c,row in
                  [(100,('1','3','1','.1','.2')),(110,('3','3','1','.1','.2'))]
                  for p in (4,40)}
        result=self.trace.paired_statistics(values,[100,110])
        self.assertGreater(result['C_RAW_p40_jackknife_se'],0)
        self.assertEqual(result['delta_C_RAW_p4_minus_p40_jackknife_se'],0)

    def primary_rows(self):
        return [dict(configuration=c,piter=p,Nv=100,tau=2,t0=97,t1=99,tmid=98,
                     axis=-1,R=4,long=2,tr1=0,tr2=0,long_from_midpoint=0,
                     Nsrc=3*64**3,L_RE='1',L_IM='0',LS_RE='2',LS_IM='0',S_RAW='1')
                for c in (100,110) for p in (4,40)]

    def group_fixture(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'data.csv'
            a.write_rows(path,rows)
            return a.grouped([(path,(4,40))],[100,110],(4,40),{(2,-1,4,2,0,0)})

    def test_missing_and_duplicate_configurations_rejected(self):
        rows=self.primary_rows()
        self.assertEqual(len(self.group_fixture(rows)),1)
        with self.assertRaisesRegex(ValueError,'Incomplete matched'):
            self.group_fixture(rows[:-1])
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            self.group_fixture(rows+[rows[0]])

    def test_changed_original_clover_and_nonfinite_primary_rejected(self):
        rows=self.primary_rows();rows[-1]['S_RAW']='1.1'
        with self.assertRaisesRegex(ValueError,'Original-field'):
            self.group_fixture(rows)
        rows=self.primary_rows();rows[-1]['LS_RE']='NaN'
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            self.group_fixture(rows)

    def test_paired_covariance_and_temporal_difference(self):
        samples=[[D(1),D(2)],[D(2),D(4)],[D(3),D(6)]]
        cov=self.temporal.jackknife_covariance(samples)
        with localcontext() as context:
            context.prec=50
            expected=D(4)/3
            self.assertEqual(cov[0][0],expected)
            self.assertEqual(cov[0][1],D(8)/3)
            direct=self.temporal.jackknife_covariance([[v[1]-v[0]] for v in samples])[0][0]
            self.assertLess(abs(direct-(cov[1][1]+cov[0][0]-2*cov[0][1])),D('1e-48'))

    def test_complete_cohort_bin_support_is_rebuilt(self):
        expected={1:(71,71,[]),2:(70,35,[800]),4:(68,17,[780,790,800]),
                  8:(64,8,[740,750,760,770,780,790,800])}
        for width,(n,b,dropped) in expected.items():
            plan=self.binning.make_plan(a.TARGET,width)
            self.assertEqual((plan['configurations'],plan['bins'],plan['dropped_incomplete_bin_configurations']),(n,b,dropped))
            self.assertEqual(plan['gaps'],[])
            self.assertIn(640,plan['selected_configurations'])
        old=self.binning.make_plan(a.BASE,2)
        self.assertEqual(old['configurations'],70)
        self.assertNotEqual(old['blocks'],self.binning.make_plan(a.TARGET,2)['blocks'])

    def test_preserved_coverage_and_matched_trace_scope(self):
        s=json.loads((a.RESULTS/'SUMMARY.json').read_text())
        self.assertEqual((s['configuration_count'],s['cases'],s['raw_primary_rows'],s['trace_primary_rows'],s['central_window_rows']),
                         (71,213,282864,141432,852))
        self.assertEqual(s['raw_piters'],[4,10,20,40])
        self.assertEqual(s['traceless_matched_piters'],[4,40])
        self.assertEqual(s['excluded_configurations'],[])
        coverage=list(a.read_rows(a.RESULTS/'configuration_coverage_71configs.csv'))
        self.assertEqual([int(r['configuration']) for r in coverage],a.TARGET)

    def test_baseline_regressions_cover_every_archived_value(self):
        v=json.loads((a.RESULTS/'VALIDATION.json').read_text())
        self.assertEqual(v['baseline70'],dict(raw_fields_exact=21912,trace_fields_exact=31872,temporal_fields_exact=19920))
        self.assertTrue(v['temporal_marginals_reproduced_exactly'])
        self.assertEqual(v['n640_review']['output_stream_hashes'],24)

    def test_tampered_stdout_is_rejected_before_analysis(self):
        source=a.HERE/'inputs/job5946008/campaign/case_000'
        workflow=a.module('qx71_test_archived_workflow',a.HERE/'inputs/job5946008/workflow/campaign.py')
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)/'case_000';shutil.copytree(source,out)
            with (out/'p4.out').open('a') as stream:stream.write('changed\n')
            entry=dict(configuration=640,tau=2,piters=[4,10,20,40])
            with self.assertRaisesRegex(ValueError,'output changed'):
                workflow.launch_records(out,entry)


if __name__=='__main__':
    unittest.main()
