#!/usr/bin/env python3
"""Check full-cohort planning, admission, grouped execution and submission."""
import io
import json
import os
import subprocess
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch
import campaign as c


class FullCampaignTests(unittest.TestCase):
    def test_one_full_cohort_and_no_remeasurement_of_baseline(self):
        self.assertEqual(len(c.TARGET), 70)
        self.assertEqual(len(c.BASELINE), 29)
        self.assertEqual(c.MISSING, [100]+list(range(240, 631, 10)))
        self.assertNotIn(640, c.TARGET)
        self.assertEqual(set(c.TARGET), set(c.BASELINE)|set(c.MISSING))
        self.assertFalse(set(c.BASELINE)&set(c.MISSING))
        cases = c.case_plan()
        self.assertEqual(len(cases), 123)
        slots = c.slot_plan(cases)
        self.assertEqual(len(slots), 31)
        self.assertEqual([i for slot in slots for i in slot], list(range(123)))
        self.assertTrue(all(1 <= len(slot) <= 4 for slot in slots))
        self.assertEqual(sum(len(x['piters']) for x in cases), 492)

    def test_missing_n100_audit_does_not_create_a_partial_campaign(self):
        known = {n: {} for n in c.MISSING if n != 100}
        with tempfile.TemporaryDirectory() as name:
            target = Path(name)/'full70'
            with patch.object(c.legacy, 'choose_runs', return_value=known), patch.object(c, 'folder', return_value=target):
                with redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                    c.prepare(SimpleNamespace(name='full70'))
                self.assertFalse(target.exists())

    def test_old_audit_directory_is_discovered_without_bypassing_quality(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            audit = root/'mental/validation/earlier_saved_mode_audit/results'
            audit.mkdir(parents=True)
            row = dict(configuration=100, run='/unused/n100', completion_marker_sha256='recorded-digest')
            (audit/'PREPARED.json').write_text(json.dumps(dict(Nv=100, times=list(range(95, 102)), runs=[row])))
            (audit/'SUMMARY.json').write_text(json.dumps(dict(success=True, runs=[dict(
                configuration=100, success=True, mpi_exit_code=0, maximum_relative_residual=1e-11)])))
            with patch.object(c, 'ROOT', root), patch.object(c, 'verify_sources'), \
                    patch.object(c.legacy, 'choose_runs', return_value={}):
                report, known = c.inventory()
                self.assertEqual(known[100], row)
                self.assertNotIn(100, report['missing_selected_central_time_audits'])
                (audit/'SUMMARY.json').write_text(json.dumps(dict(success=True, runs=[dict(
                    configuration=100, success=True, mpi_exit_code=0, maximum_relative_residual=2e-8)])))
                _, known = c.inventory()
                self.assertNotIn(100, known)

    def test_targeted_n100_audit_requires_successful_layout_control(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            result = root/'laplace_static/validation/qx_saved_modes_n100/results_n100'
            result.mkdir(parents=True)
            row = dict(configuration=100, run='/audited/n100', completion_marker_sha256='recorded-digest')
            prep = result/'PREPARED.json'
            prep.write_text(json.dumps(dict(Nv=100, times=list(range(95, 102)), runs=[row])))
            (result/'prepared.sha256').write_text(c.digest(prep)+'\n')
            summary = dict(success=True, n650_layout_control_passed=False, runs=[dict(
                configuration=100, success=True, mpi_exit_code=0, maximum_relative_residual=1e-11)])
            with patch.object(c, 'ROOT', root), patch.object(c, 'verify_sources'), \
                    patch.object(c.legacy, 'choose_runs', return_value={}):
                (result/'SUMMARY.json').write_text(json.dumps(summary))
                _, known = c.inventory()
                self.assertNotIn(100, known)
                summary['n650_layout_control_passed'] = True
                (result/'SUMMARY.json').write_text(json.dumps(summary))
                _, known = c.inventory()
                self.assertEqual(known[100], row)
                (result/'FAILED.json').write_text('{}')
                _, known = c.inventory()
                self.assertNotIn(100, known)

    def test_submit_uses_one_array_with_global_concurrency_limit(self):
        data = dict(cases=c.case_plan(), slots=c.slot_plan(c.case_plan()))
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)
            result = subprocess.CompletedProcess([], 0, '12345\n', '')
            with patch.object(c, 'prepared', return_value=(path, data)), patch.object(c.legacy, 'pilot_gate'), \
                    patch.object(c.subprocess, 'run', return_value=result) as launch, redirect_stdout(io.StringIO()):
                c.submit(SimpleNamespace(name='full70', pilot='pilot', concurrent=4))
                command = launch.call_args[0][0]
                self.assertIn('--array=0-30%4', command)
                self.assertEqual(launch.call_args[1]['universal_newlines'], True)
                self.assertEqual(json.loads((path/'SUBMISSION.json').read_text())['stdout'], '12345\n')

    def test_uncertain_submission_is_never_automatically_retried(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)
            (path/'SUBMISSION.json').write_text('{"state":"submitting"}')
            with patch.object(c, 'prepared', return_value=(path, {'slots': [[]]})), \
                    patch.object(c.legacy, 'pilot_gate'), patch.object(c.subprocess, 'run') as launch:
                with self.assertRaises(FileExistsError):
                    c.submit(SimpleNamespace(name='full70', pilot='pilot', concurrent=4))
                launch.assert_not_called()

    def test_grouped_worker_executes_every_case_and_prescription_once(self):
        cases = c.case_plan()
        data = dict(cases=cases, slots=c.slot_plan(cases), runs=[
            dict(configuration=100, run='/unused/n100'), dict(configuration=240, run='/unused/n240')])
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)
            (path/'PREPARED.json').write_text('{}')
            measured = []
            def launch(command, env, stdout, stderr, check):
                measured.append(tuple(int(env[k]) for k in ('QX_CFG', 'TAU_LENGTH', 'HYP_PITER')))
                return subprocess.CompletedProcess(command, 0)
            with patch.object(c, 'prepared', return_value=(path, data)), \
                    patch.dict(os.environ, SLURM_NTASKS='64', SLURM_JOB_NUM_NODES='8', SLURM_ARRAY_TASK_ID='0', SLURM_JOB_ID='12345'), \
                    patch.object(c.subprocess, 'run', side_effect=launch), \
                    patch.object(c, 'read_case', side_effect=lambda out, entry, data: {'success': True}), \
                    redirect_stdout(io.StringIO()):
                c.run(SimpleNamespace(name='full70'))
            wanted = [(entry['configuration'], entry['tau'], p) for entry in cases[:4] for p in c.PITERS]
            self.assertEqual(measured, wanted)
            self.assertEqual(len(list(path.glob('case_*/SUMMARY.json'))), 4)

    def test_n100_has_matched_measurements_and_retains_tau4_regression(self):
        key = (0, 1, 0, 0, 0)
        check = SimpleNamespace(
            read_measurement=Mock(return_value=({}, {key: {'S_RAW': '0.5'}}, [], {})),
            raw_rows=Mock(return_value=[{}]), write_csv=Mock(), RAW_COLUMNS=[], SU3_COLUMNS=[],
            regression_check=Mock())
        with patch.object(c.legacy, 'checker', return_value=check):
            for tau in (2, 4, 6):
                result = c.read_case(Path('/unused'), dict(configuration=100, tau=tau), {'stage': 'full70'})
                self.assertEqual(result['window_records'], 4)
                self.assertEqual(result['piters'], c.PITERS)
        self.assertEqual(check.read_measurement.call_count, 12)
        check.regression_check.assert_called_once()

    def test_original_measurement_package_remains_pinned(self):
        c.verify_sources()
        self.assertEqual(c.digest(c.legacy.CASE/'package.sha256'),
                         '1ec6167c70ac4fe0cc12c5d7e5fbb4162fb86d8afbd1f5e29eae916f52e3706b')


if __name__ == '__main__':
    unittest.main()
