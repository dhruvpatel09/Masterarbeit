#!/usr/bin/env python3
"""Admission and submission guards; these tests do not run MPI measurements."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import campaign as c


class CampaignTests(unittest.TestCase):
    def test_controls_and_failed_configuration_excluded(self):
        for text in ('100','640','110,110','111','810'):
            with self.assertRaises(ValueError):c.configs(text)
        selected=c.configs(','.join(str(x) for x in range(240,801,10) if x!=640))
        plan=c.case_plan(selected,'extended')
        self.assertEqual(len(selected),56)
        self.assertEqual(len(plan),168)
        self.assertTrue(all(r['piters']==[4,10,20,40] for r in plan))

    def test_matched_intermediate_cohort(self):
        selected=list(range(110,231,10))
        plan=c.case_plan(selected,'intermediate13')
        self.assertEqual(len(plan),39)
        self.assertTrue(all(r['piters']==[10,20] for r in plan))
        with self.assertRaises(ValueError):c.case_plan(selected[:-1],'intermediate13')
        with self.assertRaises(ValueError):c.case_plan(selected,'pilot')

    def test_no_submission_after_ambiguous_intent(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            (folder/'submission_000.json').write_text('{"state":"submitting"}')
            args=type('Args',(),dict(name='pilot',start=0))()
            with patch.object(c,'prepared',return_value=(folder,{'stage':'pilot','cases':[{}]})),patch.object(c.subprocess,'run') as launch:
                with self.assertRaises(FileExistsError):c.submit(args)
                launch.assert_not_called()

    def test_pilot_requires_all_prescriptions_and_regression(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);out=folder/'case_000';out.mkdir()
            (out/'SUMMARY.json').write_text('{"success":true}')
            (out/'LAUNCHES.json').write_text(json.dumps([dict(configuration=110,tau=4,piter=p,exit_code=0) for p in c.PITERS]))
            data={'stage':'pilot','cases':[{'configuration':110}]}
            with patch.object(c,'prepared',return_value=(folder,data)),patch.object(c,'read_case') as science:
                with self.assertRaises(ValueError):c.pilot_gate('pilot')
                science.assert_not_called()

    def test_source_pins_match_repository(self):
        c.verify_sources()
        chosen=c.choose_runs()
        self.assertTrue(set(range(110,231,10)).issubset(chosen))
        self.assertTrue(set(range(650,801,10)).issubset(chosen))


if __name__=='__main__':unittest.main()
