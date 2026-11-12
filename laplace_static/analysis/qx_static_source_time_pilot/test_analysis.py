"""Meaningful guards for frozen-source pilot analysis and archive replay."""
import ast
import math
import sys
import tempfile
import unittest
from pathlib import Path

from analysis import summarize_projection
from common import HERE, normalized_matrix, delta_matrix, frobenius, phase, ratio, csv_write
from replay import close_summary
from reproduce import compare_table


class ProjectionTests(unittest.TestCase):
    def test_common_normalization_retains_amplitude_change(self):
        row = summarize_projection(2 - 3j, 1 + 2j, 8 + 0j, 4 + 0j)
        self.assertEqual(row['common_t90_normalized_re'], .5)
        self.assertEqual(row['own_equal_time_normalized_re'], .25)
        self.assertEqual(row['delta_Q_RE'], 1)
        self.assertEqual(row['delta_Q_IM'], -5)

    def test_signed_values_are_not_filtered(self):
        row = summarize_projection(-2 + 0j, 1 + 0j, 1 + 0j, 1 + 0j)
        self.assertEqual(row['real_sign'], -1)
        self.assertEqual(row['common_t90_normalized_re'], -2)
        self.assertAlmostEqual(abs(row['phase']), math.pi)

    def test_zero_reference_is_explicitly_undefined(self):
        row = summarize_projection(1j, 0j, 1 + 0j, 1 + 0j)
        self.assertFalse(row['relative_change_defined'])
        self.assertIsNone(row['complex_relative_change'])
        self.assertIsNone(phase(0j))
        self.assertIsNone(ratio(1, 0))

    def test_equal_time_scale_must_be_positive(self):
        with self.assertRaises(ValueError):
            summarize_projection(1j, 0j, -1 + 0j, 1 + 0j)

    def test_relative_matrix_norm_keeps_complex_orientation(self):
        A = [[complex(i + 1, j + 2) for j in range(4)] for i in range(4)]
        B = [[2 * z for z in row] for row in A]
        normalized = normalized_matrix(A, [1., 4., 9., 16.])
        self.assertEqual(normalized[0][1], A[0][1] / 2)
        self.assertAlmostEqual(frobenius(delta_matrix(B, A)), frobenius(A))


class ReplayGuards(unittest.TestCase):
    def test_only_named_pivot_gets_roundoff_comparison(self):
        close_summary({'minimum_equal_time_relative_pivot': .01 + 1e-16},
                      {'minimum_equal_time_relative_pivot': .01})
        with self.assertRaises(ValueError):
            close_summary({'rows': 8447}, {'rows': 8448})
        with self.assertRaises(ValueError):
            close_summary({'minimum_equal_time_relative_pivot': .02},
                          {'minimum_equal_time_relative_pivot': .01})

    def test_table_masks_and_tiny_control_values_are_not_swallowed(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'control.csv'
            csv_write(path, [dict(delta=1e-24, missing=None)])
            compare_table(path, [dict(delta=1e-24, missing=None)])
            with self.assertRaises(ValueError):
                compare_table(path, [dict(delta=2e-24, missing=None)])
            with self.assertRaises(ValueError):
                compare_table(path, [dict(delta=1e-24, missing=0.)])

    def test_analysis_sources_support_python36_syntax(self):
        for path in HERE.glob('*.py'):
            options = {'feature_version': 6} if sys.version_info >= (3, 8) else {}
            ast.parse(path.read_text(), filename=str(path), **options)

    def test_production_commands_are_absent_from_analysis(self):
        for name in ('analysis.py', 'run_analysis.py', 'render_plots.py'):
            text = (HERE / name).read_text()
            for token in ('sbatch', 'mpiexec', 'campaign.py submit'):
                self.assertNotIn(token, text)


if __name__ == '__main__':
    unittest.main()
