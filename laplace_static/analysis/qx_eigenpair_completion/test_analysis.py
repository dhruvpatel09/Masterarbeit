#!/usr/bin/env python3
"""Reparse measured exports and check coverage and sampled-quality reductions."""
import unittest
import analyze


class AnalysisTests(unittest.TestCase):
    def test_real_archived_outputs(self):
        summary,generation,audits=analyze.verify_and_reduce()
        self.assertEqual(summary['generated_count'],70)
        self.assertEqual(summary['newly_generated_count'],13)
        self.assertEqual(summary['freshly_audited_count'],16)
        self.assertEqual(summary['audited_mode_records'],20800)
        self.assertEqual(summary['audited_slice_records'],208)
        self.assertEqual(summary['paired_profile_count'],13)
        self.assertEqual(summary['target_profile_count'],69)
        self.assertNotIn(640,summary['generated_configurations'])
        self.assertEqual(summary['modes_above_requested_1e_minus9'],0)
        self.assertEqual(max(audits,key=lambda r:r['max_relative_residual'])['configuration'],680)
        self.assertLess(summary['max_central_relative_residual'],1.2e-11)


if __name__=='__main__':unittest.main()
