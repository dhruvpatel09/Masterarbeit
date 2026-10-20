#!/usr/bin/env python3
"""Check preserved numerics, input views and the n650-before-n100 audit gate."""
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import audit as a


class AuditTests(unittest.TestCase):
    def test_only_process_layout_differs_from_preserved_c_checker(self):
        expected = (a.PRESERVED/'saved_mode_audit.c').read_text()
        for before, after in [('ranks!=384', 'ranks!=64'), ('need 384 ranks', 'need 64 ranks'),
                              ('P[4]={6,4,4,4}', 'P[4]={1,4,4,4}'), ('grid=6x4x4x4', 'grid=1x4x4x4')]:
            self.assertEqual(expected.count(before), 1)
            expected = expected.replace(before, after)
        self.assertEqual((a.HERE/'test_saved_modes.c').read_text(), expected)
        a.verify_sources()

    def test_real_preserved_control_passes_and_wrong_spectrum_is_rejected(self):
        passed, modes, slices = a.parser.parse(a.PRESERVED/'job5925408/n650.out', 650)
        self.assertTrue(passed)
        self.assertEqual((len(modes), len(slices)), (1300, 13))
        a.control_regression(modes)
        altered = [dict(row) for row in modes]
        altered[0]['lambda'] += 1e-3
        with self.assertRaises(ValueError):
            a.control_regression(altered)

    def make_inputs(self, root):
        gauges = root/'gauges'
        gauges.mkdir()
        rows = {}
        for cfg, times in [(100, range(192)), (650, a.TIMES)]:
            run = root/('n%d' % cfg)
            (run/'eigenvectors').mkdir(parents=True)
            (run/'eigenvalues').mkdir()
            (gauges/('qX_p2n%d' % cfg)).write_bytes(b'input inventory fixture')
            for t in times:
                (run/'eigenvectors'/('qX_p2n%d_evec_t%d.h5' % (cfg, t))).write_bytes(b'inventory fixture')
            (run/'eigenvalues'/'values.h5').write_bytes(b'inventory fixture')
            if cfg == 650:
                marker = run/'COMPLETED.json'
                marker.write_text('{"configuration":650,"timeslices":192,"eigenvector_files":192}')
                rows[650] = dict(configuration=650, run=str(run), completion_marker_sha256=a.digest(marker))
        return gauges, rows

    def test_preparation_preserves_original_n100_files_and_creates_only_a_view(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            gauges, rows = self.make_inputs(root)
            results = root/'results_n100'
            with patch.object(a, 'RESULTS', results), patch.object(a, 'verify_build', return_value='build'), \
                    patch.object(a.legacy, 'build_manifest', return_value='observable'), \
                    patch.object(a.paths, 'N100_RUN', root/'n100'), \
                    patch.object(a.legacy, 'GAUGES', gauges), patch.object(a.legacy, 'choose_runs', return_value=rows):
                a.prepare(SimpleNamespace(eigenvalues=None))
            self.assertFalse((root/'n100/COMPLETED.json').exists())
            self.assertTrue((results/'input_run/eigenvectors').is_symlink())
            self.assertEqual((results/'input_run/eigenvectors').resolve(), root/'n100/eigenvectors')
            marker = json.loads((results/'input_run/COMPLETED.json').read_text())
            self.assertEqual((marker['configuration'], marker['timeslices'], marker['eigenvector_files']), (100, 192, 192))
            self.assertFalse((results/'SUMMARY.json').exists())
            prepared = json.loads((results/'PREPARED.json').read_text())
            self.assertEqual(prepared['configurations'], [650, 100])
            self.assertEqual(prepared['times'], a.TIMES)

    def test_missing_n100_file_stops_before_creating_results(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            gauges, rows = self.make_inputs(root)
            (root/'n100/eigenvectors/qX_p2n100_evec_t101.h5').unlink()
            results = root/'results_n100'
            with patch.object(a, 'RESULTS', results), patch.object(a, 'verify_build', return_value='build'), \
                    patch.object(a.legacy, 'build_manifest', return_value='observable'), \
                    patch.object(a.paths, 'N100_RUN', root/'n100'), \
                    patch.object(a.legacy, 'GAUGES', gauges), patch.object(a.legacy, 'choose_runs', return_value=rows):
                with self.assertRaises(FileNotFoundError):
                    a.prepare(SimpleNamespace(eigenvalues=None))
            self.assertFalse(results.exists())

    def test_failed_layout_control_prevents_launching_n100(self):
        with tempfile.TemporaryDirectory() as name:
            results = Path(name)
            data = dict(runs=[dict(configuration=cfg, gauge='/unused', eigenvectors='/unused', eigenvalues='/unused')
                              for cfg in (650, 100)])
            with patch.object(a, 'RESULTS', results), patch.object(a, 'prepared', return_value=data), \
                    patch.dict(os.environ, SLURM_NTASKS='64', SLURM_JOB_NUM_NODES='8'), \
                    patch.object(a.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0)) as launch, \
                    patch.object(a.parser, 'parse', return_value=(False, [], [])):
                with self.assertRaises(ValueError):
                    a.run(SimpleNamespace())
                self.assertEqual(launch.call_count, 1)
                self.assertFalse((results/'n100.out').exists())
                self.assertFalse((results/'SUMMARY.json').exists())
                self.assertTrue((results/'FAILED.json').exists())

    def test_ambiguous_submission_is_not_retried(self):
        with tempfile.TemporaryDirectory() as name:
            results = Path(name)
            (results/'SUBMISSION.json').write_text('{"state":"submitting"}')
            with patch.object(a, 'RESULTS', results), patch.object(a, 'prepared'), \
                    patch.object(a.subprocess, 'run') as launch:
                with self.assertRaises(FileExistsError):
                    a.submit(SimpleNamespace())
                launch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
