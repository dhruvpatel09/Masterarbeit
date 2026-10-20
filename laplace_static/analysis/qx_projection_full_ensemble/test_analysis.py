#!/usr/bin/env python3
"""Check estimator continuity and explicit full-cohort inclusion rules."""
import unittest
from decimal import Decimal
from pathlib import Path
import make_plots as a


class FullAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base = a.module('qx_full_test29', a.ROOT/'laplace_static/analysis/qx_projection_spectrum_29configs/make_plots.py')
        _, cls.raw, _, cls.reference, _ = base.load_data()

    def test_preserved_29_config_estimates_and_paired_errors_are_identical(self):
        records, configs, piters = a.aggregate_records(self.raw)
        self.assertEqual(records, self.reference)
        self.assertEqual(len(configs), 29)
        self.assertEqual(piters, a.PITERS)

    def test_n100_is_a_measured_member_in_the_new_explicit_cohort(self):
        first = self.raw[0]
        key = tuple(first[k] for k in a.KEYS)
        selected = [dict(r) for r in self.raw if int(r['configuration']) == 110 and
                    tuple(r[k] for k in a.KEYS) == key]
        self.assertEqual(len(selected), 4)
        added = [dict(r, configuration='100') for r in selected]
        records, configs, _ = a.aggregate_records(selected+added)
        self.assertEqual(configs, [100, 110])
        self.assertEqual(records[0]['configurations'], 2)
        self.assertEqual(records[0]['C_p4_jackknife_se'], Decimal(0))
        with self.assertRaises(ValueError):
            a.core.aggregate(selected+added)
        with self.assertRaises(ValueError):
            a.aggregate_records(selected+added[:-1])

    def test_n640_cannot_be_silently_added(self):
        with self.assertRaises(ValueError):
            a.aggregate_records([dict(self.raw[0], configuration='640')])

    def test_partial_ensemble_cannot_produce_full70_results(self):
        with self.assertRaisesRegex(ValueError, 'Incomplete full ensemble'):
            a.load_data([], [])


if __name__ == '__main__':
    unittest.main()
