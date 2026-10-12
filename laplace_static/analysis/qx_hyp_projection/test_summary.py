#!/usr/bin/env python3
"""Scientific estimator/cohort checks and regression to the preserved analysis."""
import copy
import unittest
from decimal import Decimal

import make_plots as m


class SummaryChecks(unittest.TestCase):
    def fixture(self):
        # Equal raw weights give C=(2+8)/(1+3)-1=1.5.
        # Mean of connected configuration ratios would incorrectly give 4/3.
        return [dict(configuration=c, piter=p, tau=2, axis=-1, R=4,
                     long=2, tr1=0, tr2=0, L_RE=str(loop), LS_RE=str(probe), S_RAW='1')
                for c, loop, probe in ((110, 1, 2), (120, 3, 8)) for p in (4, 40)]

    def test_ratio_of_raw_means(self):
        rows, configs, piters = m.aggregate(self.fixture())
        self.assertEqual(rows[0]['C_p40'], Decimal('1.5'))
        self.assertEqual(rows[0]['delta_C_p4_minus_p40'], 0)
        self.assertGreater(rows[0]['C_p40_jackknife_se'], 0)
        self.assertEqual(rows[0]['delta_C_paired_jackknife_se'], 0)
        self.assertEqual((configs, piters), ([110, 120], [4, 40]))

    def test_incomplete_and_duplicate_cohorts(self):
        rows = self.fixture()
        with self.assertRaisesRegex(ValueError, 'Incomplete matched cohort'):
            m.aggregate(rows[:-1])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            m.aggregate(rows+[copy.deepcopy(rows[0])])

    def test_unmeasured_intermediate_curve_rejected(self):
        rows = self.fixture()
        extra = dict(rows[0], piter=20)
        with self.assertRaisesRegex(ValueError, 'Incomplete matched cohort'):
            m.aggregate(rows+[extra])

    def test_clover_and_control_scope(self):
        rows = self.fixture()
        rows[0]['S_RAW'] = '1.01'
        with self.assertRaisesRegex(ValueError, 'clover changed'):
            m.aggregate(rows)
        with self.assertRaisesRegex(ValueError, 'regression control'):
            m.aggregate([dict(r, configuration=100) for r in self.fixture()])

    def test_all_preserved_ensemble_values(self):
        sources = m.check_sources()
        raw = [r for s in sources['default_profile_tables'] for r in m.read_csv(m.REPO/s['path'])]
        rows, configs, piters = m.aggregate(raw)
        lookup = {tuple(r[k] for k in m.KEYS): r for r in rows}
        original = m.read_csv(m.REPO/'laplace_static/analysis/hyp_projection_13configs/results/ensemble_13_combined.csv')
        self.assertEqual((len(rows), len(original)), (996, 249))
        self.assertEqual((len(configs), piters), (13, [4, 40]))
        for old in original:
            new = lookup[tuple(int(old[k]) for k in m.KEYS)]
            for field in ('C_p4', 'C_p40', 'delta_C_p4_minus_p40',
                          'C_p4_jackknife_se', 'C_p40_jackknife_se', 'delta_C_paired_jackknife_se'):
                self.assertEqual(new[field], Decimal(old[field]))


if __name__ == '__main__':
    unittest.main()
