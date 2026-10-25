#!/usr/bin/env python3
"""Independent covariance and ratio tests for the paired temporal comparison."""
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction

import reproduce as analysis

D = Decimal


class TemporalStatisticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source, _ = analysis.verified_upstream()

    def values(self, series):
        return {t: {(c, p): (D(1), D(series[t][i]), D(0), D(0), D(0))
                    for i, c in enumerate((100, 110, 120, 130)) for p in (4, 40)}
                for t in (2, 4, 6)}

    def test_linear_mean_covariance_and_correlated_contrasts(self):
        series = {2: [1, 3, 4, 8], 4: [8, 6, 5, 1], 6: [2, 4, 5, 9]}
        point, loo = analysis.paired_estimates(self.values(series), [100, 110, 120, 130], self.source)
        result, cov = analysis.channel_statistics(point, loo, 'RAW', 40)
        # For constant denominators the estimator is a mean: covariance is sample covariance / N.
        with localcontext() as context:
            context.prec = 50
            for i, t in enumerate((2, 4, 6)):
                for j, u in enumerate((2, 4, 6)):
                    mt, mu = sum(map(D, series[t])) / 4, sum(map(D, series[u])) / 4
                    exact = sum(((D(a) - mt) * (D(b) - mu) for a, b in zip(series[t], series[u])), D(0)) / 12
                    self.assertLess(abs(cov[i][j] - exact), D('1e-46'))
            # tau6 = tau2 + 1 on every configuration: the paired difference has zero error.
            self.assertEqual(result['delta_tau6_minus_tau2'], D(1))
            self.assertEqual(result['delta_tau6_minus_tau2_jackknife_se'], D(0))
            self.assertGreater(result['delta_tau6_minus_tau2_se_ignoring_tau_covariance'], 0)
            # tau4 is perfectly anticorrelated with tau2: the difference variance is four times one marginal.
            self.assertLess(abs(result['delta_tau4_minus_tau2_jackknife_se'] ** 2 - 4 * cov[0][0]), D('1e-46'))

    def test_nonlinear_ratio_against_independent_rational_arithmetic(self):
        labels, denominators = [100, 110, 120, 130], [1, 2, 4, 8]
        numerators = {2: [2, 3, 7, 4], 4: [4, 1, 2, 9], 6: [1, 5, 3, 2]}
        values = {t: {(c, p): (D(denominators[i]), D(numerators[t][i]), D(i), D(0), D(0))
                      for i, c in enumerate(labels) for p in (4, 40)} for t in (2, 4, 6)}
        point, loo = analysis.paired_estimates(values, labels, self.source)
        stats, _ = analysis.channel_statistics(point, loo, 'RAW', 40)
        full = {t: Fraction(sum(numerators[t]), sum(denominators)) - Fraction(6, 4) for t in (2, 4, 6)}
        rational_differences = []
        for i in range(4):
            estimate = {t: Fraction(sum(numerators[t]) - numerators[t][i], sum(denominators) - denominators[i])
                           - Fraction(6 - i, 3) for t in (2, 4, 6)}
            rational_differences.append(estimate[6] - estimate[4])
        center = sum(rational_differences) / 4
        exact_var = Fraction(3, 4) * sum((x - center) ** 2 for x in rational_differences)
        with localcontext() as context:
            context.prec = 50
            convert = lambda x: D(x.numerator) / D(x.denominator)
            self.assertLess(abs(stats['C_tau2'] - convert(full[2])), D('1e-46'))
            self.assertLess(abs(stats['delta_tau6_minus_tau4'] - convert(full[6] - full[4])), D('1e-46'))
            self.assertLess(abs(stats['delta_tau6_minus_tau4_jackknife_se'] ** 2 - convert(exact_var)), D('1e-46'))
            incorrect_mean_of_ratios = sum((Fraction(a, b) for a, b in zip(numerators[2], denominators))) / 4 - Fraction(6, 4)
            self.assertNotEqual(full[2], incorrect_mean_of_ratios)

    def test_unmatched_tau_cohort_rejected(self):
        values = self.values({2: [1, 2, 3, 4], 4: [2, 3, 4, 5], 6: [3, 4, 5, 6]})
        del values[4][110, 40]
        with self.assertRaisesRegex(ValueError, 'identical matched cohort'):
            analysis.paired_estimates(values, [100, 110, 120, 130], self.source)

    def test_zero_leave_one_out_denominator_rejected(self):
        values = self.values({2: [1, 2, 3, 4], 4: [2, 3, 4, 5], 6: [3, 4, 5, 6]})
        for p in (4, 40):
            for c, L in zip((100, 110, 120, 130), (1, 1, -2, 3)):
                values[2][c, p] = (D(L), D(1), D(0), D(0), D(0))
        with self.assertRaisesRegex(ValueError, 'loop ratio'):
            analysis.paired_estimates(values, [100, 110, 120, 130], self.source)

    def test_identical_projection_controls_cancel_with_zero_error(self):
        values = self.values({2: [1, 3, 4, 8], 4: [8, 6, 5, 1], 6: [2, 4, 5, 9]})
        point, loo = analysis.paired_estimates(values, [100, 110, 120, 130], self.source)
        for definition in ('RAW', 'TL'):
            for row in analysis.projection_controls(point, loo, definition):
                self.assertEqual(row['contrast_p4_minus_p40'], 0)
                self.assertEqual(row['contrast_p4_minus_p40_jackknife_se'], 0)

    def test_projection_contrast_preserves_four_way_pairing(self):
        values = self.values({2: [1, 3, 4, 8], 4: [8, 6, 5, 1], 6: [2, 4, 5, 9]})
        labels, shifts, coefficients = (100, 110, 120, 130), (2, 5, 8, 9), {2: 1, 4: 3, 6: 4}
        for t in (2, 4, 6):
            for c, shift in zip(labels, shifts):
                old = values[t][c, 4]
                values[t][c, 4] = (old[0], old[1] + coefficients[t] * D(shift)) + old[2:]
        point, loo = analysis.paired_estimates(values, list(labels), self.source)
        controls = analysis.projection_controls(point, loo, 'RAW')
        with localcontext() as context:
            context.prec = 50
            for row in controls:
                coefficient = D(coefficients[row['tau_high']] - coefficients[row['tau_low']])
                self.assertEqual(row['contrast_p4_minus_p40'], 6 * coefficient)
                expected_variance = D('2.5') * coefficient ** 2
                self.assertLess(abs(row['contrast_p4_minus_p40_jackknife_se'] ** 2 - expected_variance), D('1e-46'))

    def test_nonfinite_covariance_input_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            analysis.jackknife_covariance([[D(1)], [D('NaN')]])


if __name__ == '__main__':
    unittest.main()
