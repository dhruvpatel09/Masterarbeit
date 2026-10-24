#!/usr/bin/env python3
"""Check nonlinear averaging, paired covariance and archived raw regressions."""
import importlib.util
import math
import unittest
from decimal import Decimal as D
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('clover_full70',str(HERE/'reproduce.py'))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class AnalysisTests(unittest.TestCase):
    def synthetic(self):
        values = {}
        for cfg, entries in enumerate(((1,2,1,'.1','.3'),(2,8,1,'.2','.6'),(3,18,1,'.3','.9'))):
            for p in a.PITERS:
                values[cfg,p] = tuple(D(str(x)) for x in entries)
        return values

    def test_ratio_of_means_not_mean_of_ratios(self):
        result = a.paired_statistics(self.synthetic(),[0,1,2])
        self.assertAlmostEqual(float(result['C_RAW_p40']),28/6-1,places=14)
        self.assertNotAlmostEqual(float(result['C_RAW_p40']),(1+3+5)/3,places=10)

    def test_identical_prescriptions_have_zero_paired_error(self):
        result = a.paired_statistics(self.synthetic(),[0,1,2])
        self.assertGreater(result['C_RAW_p40_jackknife_se'],0)
        for definition in ('RAW','TRACE','TL'):
            self.assertEqual(result['delta_C_%s_p4_minus_p40' % definition],0)
            self.assertEqual(result['delta_C_%s_p4_minus_p40_jackknife_se' % definition],0)

    def test_paired_trace_and_fraction_errors_against_float_reference(self):
        values = self.synthetic()
        result = a.paired_statistics(values,[0,1,2])
        traces, fractions = [], []
        for omitted in range(3):
            selected = [tuple(float(x) for x in values[c,40]) for c in range(3) if c != omitted]
            L,LS,S,T,LT = [sum(row[j] for row in selected) for j in range(5)]
            raw, trace = LS/L-S/2,LT/L-T/2
            traces.append(trace)
            fractions.append(100*trace/raw)
        for field,samples in [('C_TRACE_p40',traces),('trace_fraction_percent_p40',fractions)]:
            center = sum(samples)/3
            expected = math.sqrt(2/3*sum((x-center)**2 for x in samples))
            self.assertAlmostEqual(float(result[field+'_jackknife_se']),expected,places=13)

    def test_duplicate_missing_and_excluded_configuration_rejected(self):
        row = a.read_csv(HERE/'inputs/profiles_raw_and_trace_70configs.csv')[0]
        with self.assertRaises(ValueError):
            a.group_primaries([row,row])
        with self.assertRaises(ValueError):
            a.group_primaries([])
        with self.assertRaises(ValueError):
            a.group_primaries([dict(row,configuration='640')])
        with self.assertRaises(ValueError):
            a.group_primaries([dict(row,T_TRACE='NaN')])

    def test_raw_estimates_and_paired_errors_match_previous_full70(self):
        previous = a.read_csv(a.ROOT/'laplace_static/analysis/qx_projection_full_ensemble/results/profiles_ensemble_70configs.csv')
        current = a.read_csv(HERE/'results/profiles_ensemble_70configs.csv')
        old = {tuple(int(r[k]) for k in a.KEYS):r for r in previous}
        mapping = [('C_RAW_p4','C_p4'),('C_RAW_p40','C_p40'),
                   ('C_RAW_p4_jackknife_se','C_p4_jackknife_se'),
                   ('C_RAW_p40_jackknife_se','C_p40_jackknife_se'),
                   ('delta_C_RAW_p4_minus_p40','delta_C_p4_minus_p40'),
                   ('delta_C_RAW_p4_minus_p40_jackknife_se','delta_C_paired_jackknife_se')]
        self.assertEqual(len(current),996)
        for row in current:
            ref = old[tuple(int(row[k]) for k in a.KEYS)]
            for new, prior in mapping:
                self.assertLess(abs(D(row[new])-D(ref[prior])),D('1e-40'))


if __name__ == '__main__':
    unittest.main()
