#!/usr/bin/env python3
"""Check pilot reuse, raw regression, complete cohorts and recovery markers."""
import copy
import importlib.util
import json
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('qx_trace_ensemble_test', str(HERE/'campaign.py'))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CampaignTests(unittest.TestCase):
    def test_sources_and_validated_seed(self):
        c.verify_sources()
        raw, windows = c.seed_data('full70')
        self.assertEqual((len(raw), len(windows)), (664, 2))
        self.assertEqual({int(r['configuration']) for r in raw}, {100})
        self.assertEqual({int(r['tau']) for r in raw}, {4})
        self.assertEqual({int(r['piter']) for r in raw}, {4, 40})

    def test_full_plan_has_418_new_launches_and_35_slots(self):
        cases, slots = c.case_plan('full70'), c.slot_plan('full70')
        self.assertEqual(len(c.TARGET), 70)
        self.assertNotIn(640, c.TARGET)
        self.assertEqual(len(cases), 209)
        self.assertEqual(sum(len(r['piters']) for r in cases), 418)
        self.assertEqual(len(slots), 35)
        self.assertEqual([i for slot in slots for i in slot], list(range(209)))
        self.assertTrue(all(1 <= len(slot) <= 6 for slot in slots))
        self.assertNotIn((100, 4), {(r['configuration'], r['tau']) for r in cases})

    def test_spectrum_plan_is_explicitly_a_selected_central_diagnostic(self):
        cases = c.case_plan('spectrum')
        self.assertEqual({r['configuration'] for r in cases}, {100, 170, 670, 740})
        self.assertEqual({r['tau'] for r in cases}, {2, 4, 6})
        self.assertEqual(sum(len(r['piters']) for r in cases), 36)
        self.assertTrue(all(r['piters'] == [6, 8, 80] for r in cases))
        self.assertEqual(len(c.slot_plan('spectrum')), 12)

    def test_pilot_replay_and_raw_reference_rejection(self):
        entry = dict(configuration=100, tau=4, piters=[4, 40])
        refs = c.references([entry])
        raw, windows = c.read_case(c.analysis.INPUT/'outputs', entry, refs)
        expected, _ = c.seed_data('full70')
        self.assertEqual(raw, expected)
        self.assertEqual(len(windows), 2)
        corrupt = copy.deepcopy(refs)
        key = next(k for k in corrupt if k[2] == 4)
        corrupt[key]['L_RE'] = str(float(corrupt[key]['L_RE'])*1.001)
        with self.assertRaises(RuntimeError):
            c.read_case(c.analysis.INPUT/'outputs', entry, corrupt)

    def test_partial_seed_cannot_be_labelled_full_ensemble(self):
        raw, windows = c.seed_data('full70')
        with self.assertRaises(ValueError):
            c.assert_cohort(raw, windows, 'full70')

    def test_complete_spectrum_cohort_and_duplicate_detection(self):
        seed, _ = c.seed_data('full70')
        geometry = [r for r in seed if r['piter'] == '40']
        raw, windows = [], []
        for entry in c.case_plan('spectrum'):
            for p in entry['piters']:
                group = dict(configuration=str(entry['configuration']), tau=str(entry['tau']), piter=str(p))
                raw += [dict(row, **group) for row in geometry]
                windows.append(group)
        c.assert_cohort(raw, windows, 'spectrum')
        with self.assertRaises(ValueError):
            c.assert_cohort(raw+[raw[0]], windows, 'spectrum')
        with self.assertRaises(ValueError):
            c.assert_cohort(raw, windows[:-1], 'spectrum')

    def test_only_validated_attempts_are_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)
            case = target/'case_000'
            attempt = case/'attempt_0001'
            attempt.mkdir(parents=True)
            (attempt/'FAILED.json').write_text('{"success": false}\n')
            self.assertIsNone(c.current_attempt(target, 0))
            (attempt/'SUMMARY.json').write_text('{"success": true}\n')
            (case/'LATEST.json').write_text(json.dumps(dict(attempt=attempt.name,
                                                          summary_sha256=c.digest(attempt/'SUMMARY.json'))))
            self.assertEqual(c.current_attempt(target, 0), attempt)
            (attempt/'SUMMARY.json').write_text('{"success": false}\n')
            with self.assertRaises(ValueError):
                c.current_attempt(target, 0)

    def test_invalid_names_and_stage_rejected(self):
        with self.assertRaises(ValueError):
            c.folder('../escape')
        with self.assertRaises(ValueError):
            c.case_plan('unknown')

    def test_failed_attempt_is_preserved_then_resumed_without_repeating_success(self):
        entry = dict(index=0, configuration=100, tau=4, piters=[4, 40])
        data = dict(cases=[entry], slots=[[0]], runs=[dict(configuration=100, run='/unused/audited_run')])
        environment = dict(SLURM_NTASKS='64', SLURM_JOB_NUM_NODES='8',
                           SLURM_ARRAY_TASK_ID='0', SLURM_JOB_ID='5942805')
        def launcher(fail):
            def invoke(command, env, stdout, stderr, check):
                p = int(env['HYP_PITER'])
                stdout.write((c.analysis.INPUT/'outputs'/('p%d.out' % p)).read_text())
                return SimpleNamespace(returncode=1 if fail and p == 40 else 0)
            return invoke
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)
            (target/'PREPARED.json').write_text('{}\n')
            args = SimpleNamespace(name='test')
            with mock.patch.object(c, 'prepared', return_value=(target, data)), mock.patch.dict(c.os.environ, environment):
                with mock.patch.object(c.subprocess, 'run', side_effect=launcher(True)):
                    with self.assertRaises(ValueError):
                        c.run(args)
                self.assertIsNone(c.current_attempt(target, 0))
                self.assertTrue((target/'case_000/attempt_0001/FAILED.json').is_file())
                with mock.patch.object(c.subprocess, 'run', side_effect=launcher(False)):
                    c.run(args)
                self.assertEqual(c.current_attempt(target, 0).name, 'attempt_0002')
                with mock.patch.object(c.subprocess, 'run', side_effect=AssertionError('Completed launches repeated')):
                    c.run(args)


if __name__ == '__main__':
    unittest.main()
