#!/usr/bin/env python3
"""Independent equal-bin jackknife and gap-handling checks."""
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction
import reproduce as analysis

D = Decimal


class BinningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporal, cls.source, _ = analysis.load_dependencies()

    def values(self, series):
        return {t: {(100 + i * 10, p): (D(1), D(x), D(0), D(0), D(0))
                    for i, x in enumerate(series[t]) for p in (4, 40)} for t in (2, 4, 6)}

    def test_full70_bin_support_and_gap_preservation(self):
        labels = [c for c in range(100, 801, 10) if c != 640]
        expected = {1: (70, 70), 2: (70, 35), 4: (68, 17), 8: (64, 8)}
        for width, counts in expected.items():
            plan = analysis.make_plan(labels, width)
            self.assertEqual((plan['configurations'], plan['bins']), counts)
            self.assertEqual(plan['gaps'], [{'left': 630, 'right': 650}])
            self.assertEqual(sorted(plan['selected_configurations'] + plan['dropped_incomplete_bin_configurations']), labels)
            for block in plan['blocks']:
                self.assertTrue(all(b - a == 10 for a, b in zip(block, block[1:])))
        self.assertEqual(analysis.make_plan(labels, 4)['dropped_incomplete_bin_configurations'], [620, 630])
        self.assertEqual(analysis.make_plan(labels, 8)['dropped_incomplete_bin_configurations'], list(range(580, 640, 10)))

    def test_equal_bin_mean_error_matches_sample_variance_of_bin_means(self):
        values = self.values({2: [1, 2, 4, 6, 7, 10], 4: [2, 3, 5, 7, 8, 11], 6: [3, 4, 6, 8, 9, 12]})
        blocks = [[100, 110], [120, 130], [140, 150]]
        point, loo = analysis.block_estimates(values, blocks, self.source)
        result, _ = self.temporal.channel_statistics(point, loo, 'RAW', 40)
        with localcontext() as context:
            context.prec = 50
            means = [D('1.5'), D(5), D('8.5')]
            center = sum(means) / 3
            expected_variance = sum((x - center) ** 2 for x in means) / 6
            self.assertLess(abs(result['C_tau2_jackknife_se'] ** 2 - expected_variance), D('1e-45'))
            self.assertEqual(result['delta_tau6_minus_tau4'], 1)
            self.assertEqual(result['delta_tau6_minus_tau4_jackknife_se'], 0)

    def test_nonlinear_ratio_and_error_match_rational_delete_bin_calculation(self):
        labels, L, LS = [100, 110, 120, 130, 140, 150], [1, 2, 4, 8, 2, 3], [2, 1, 5, 6, 8, 4]
        values = {t: {(c, p): (D(L[i]), D(LS[i]), D(i), D(0), D(0))
                      for i, c in enumerate(labels) for p in (4, 40)} for t in (2, 4, 6)}
        blocks = [[100, 110], [120, 130], [140, 150]]
        point, loo = analysis.block_estimates(values, blocks, self.source)
        stats, _ = self.temporal.channel_statistics(point, loo, 'RAW', 40)
        expected_point = Fraction(sum(LS), sum(L)) - Fraction(15, 6)
        samples = []
        for i in (0, 2, 4):
            indices = [j for j in range(6) if j not in (i, i + 1)]
            samples.append(Fraction(sum(LS[j] for j in indices), sum(L[j] for j in indices))
                           - Fraction(sum(indices), 4))
        center = sum(samples) / 3
        variance = Fraction(2, 3) * sum((x - center) ** 2 for x in samples)
        with localcontext() as context:
            context.prec = 50
            dec = lambda x: D(x.numerator) / D(x.denominator)
            self.assertLess(abs(stats['C_tau2'] - dec(expected_point)), D('1e-45'))
            self.assertLess(abs(stats['C_tau2_jackknife_se'] ** 2 - dec(variance)), D('1e-45'))

    def test_unit_bins_reproduce_delete_one_configuration(self):
        values = self.values({2: [1, 2, 4, 6], 4: [2, 4, 1, 5], 6: [4, 3, 2, 1]})
        labels = [100, 110, 120, 130]
        point, loo = analysis.block_estimates(values, [[c] for c in labels], self.source)
        reference, reference_loo = self.temporal.paired_estimates(values, labels, self.source)
        self.assertEqual((point, loo), (reference, reference_loo))

    def test_duplicate_or_unequal_bins_rejected(self):
        values = self.values({2: [1, 2, 4, 6], 4: [2, 4, 1, 5], 6: [4, 3, 2, 1]})
        for blocks in ([[100, 110], [110, 120]], [[100, 110], [120]]):
            with self.assertRaises(ValueError):
                analysis.block_estimates(values, blocks, self.source)

    def test_missing_matching_tau_or_piter_rejected(self):
        values = self.values({2: [1, 2, 4, 6], 4: [2, 4, 1, 5], 6: [4, 3, 2, 1]})
        del values[4][100, 40]
        with self.assertRaisesRegex(ValueError, 'Incomplete matched'):
            analysis.block_estimates(values, [[100, 110], [120, 130]], self.source)

    def test_bin_construction_requires_order_and_enough_bins(self):
        for labels, width in [([110, 100], 1), ([100, 100], 1), ([100, 110], 4),
                              ([100, 110], 0), ([100, 110], '2'), ([100, 110], True)]:
            with self.assertRaises(ValueError):
                analysis.make_plan(labels, width)


if __name__ == '__main__':
    unittest.main()
