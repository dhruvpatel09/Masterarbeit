#!/usr/bin/env python3
"""Synthetic statistical checks and published-data regression."""
import ast
import copy
import json
import tempfile
import unittest
import sys
from pathlib import Path

import numpy as np
import analysis
import approach_model
import data
import reproduce
from statistics_tools import (assert_close, constant_fit, covariance, energy,
                              standard_error, whiten)
from decimal import Decimal as D


class StatisticalTests(unittest.TestCase):
    def test_log_after_averaging_preserves_signed_samples(self):
        x=[D('-1'),D('4'),D('6')]
        y=[D('2'),D('3'),D('4')]
        result,delta,loo,_=energy(x,y)
        self.assertAlmostEqual(float(result[0]),0.)
        self.assertEqual(len(delta),3)
        self.assertEqual(len(loo),3)
        self.assertAlmostEqual(float(result[1]),standard_error(delta))

    def test_nonpositive_mean_is_undefined(self):
        result,delta,loo,_=energy([D('-2'),D('1')],[D('1'),D('2')])
        self.assertIsNone(result)
        self.assertIsNone(delta)
        self.assertIsNone(loo)

    def test_correlated_gls_matches_direct_linear_algebra(self):
        rng=np.random.default_rng(28)
        d=rng.normal(size=(71,6))@np.array([[2,.8,.4,.1,0,0],[0,1,.3,0,0,0],
            [0,0,1,.5,0,0],[0,0,0,1,.6,0],[0,0,0,0,1,.5],[0,0,0,0,0,1]])
        y=np.array([.2,.25,.19,.28,.21,.24])
        fitted,w=constant_fit(y,d);C=covariance(d)
        expected=np.linalg.solve(C,np.ones(6));expected/=expected.sum()
        np.testing.assert_allclose(w,expected,rtol=1e-12)
        self.assertAlmostEqual(fitted['naive_se']**2,w@C@w)
        self.assertGreater(np.max(abs(w-(1/np.diag(C))/(1/np.diag(C)).sum())),.01)

    def test_paired_contrast_keeps_shared_noise(self):
        rng=np.random.default_rng(9)
        shared=rng.normal(size=71);extra=rng.normal(size=71)*.01
        a=dict(row=dict(estimate=.2,low_signal=False),fluct=shared+extra,loo=None)
        b=dict(row=dict(estimate=.19,low_signal=False),fluct=shared,loo=None)
        row=analysis.contrast(a,b)
        self.assertAlmostEqual(row['paired_naive_se'],standard_error(extra))
        self.assertLess(row['paired_naive_se'],standard_error(shared)/20)

    def test_singular_covariance_is_rejected(self):
        with self.assertRaises((ValueError,np.linalg.LinAlgError)):
            whiten(np.ones((4,4)))

    def test_nonfinite_never_passes_replay(self):
        with self.assertRaises(ValueError):assert_close(float('nan'),0,'nan')
        with self.assertRaises(ValueError):assert_close(float('inf'),0,'inf')
        with self.assertRaises(ValueError):reproduce.scalar('nan',1.,'field')
        with self.assertRaises(ValueError):reproduce.scalar('0',None,'undefined')

    def test_reported_stromboli_curvature_condition_roundoff(self):
        reproduce.scalar('61138888376.26336',61139006697.25453,
                         'reported row 32',field='curvature_condition')

    def test_curvature_roundoff_does_not_relax_parameter_or_error_checks(self):
        for field in ('A','B','c','A_conditional_se','B_conditional_se','c_conditional_se'):
            with self.assertRaises(ValueError):
                reproduce.scalar('61138888376.26336',61139006697.25453,
                                 field,field=field)

    def test_material_curvature_condition_changes_still_fail(self):
        with self.assertRaises(ValueError):
            reproduce.curvature_condition(6e10,6.6e10,'material change')
        with self.assertRaises(ValueError):
            reproduce.curvature_condition(1000.,1000.01,'well conditioned')

    def test_curvature_classification_and_invalid_values_still_fail(self):
        with self.assertRaises(ValueError):
            reproduce.curvature_condition(1e12*(1.-1e-8),1e12*(1.+1e-8),'threshold')
        for bad in ('nan','inf','0','-1','0.9'):
            with self.assertRaises(ValueError):
                reproduce.curvature_condition(bad,6e10,'invalid')
        with self.assertRaises(ValueError):
            reproduce.scalar('',6e10,'unexpected undefined',field='curvature_condition')
        with self.assertRaises(ValueError):
            reproduce.scalar('6e10',None,'unexpected defined',field='curvature_condition')

    def test_curvature_roundoff_matches_small_matrix_perturbation(self):
        q,_=np.linalg.qr(np.array([[1.,2.,3.],[-2.,1.,4.],[4.,-3.,2.]]))
        h=q@np.diag([1.,1e-5,1e-11])@q.T
        perturb=np.array([[0.,1.,-1.],[1.,0.,1.],[-1.,1.,0.]])*np.finfo(float).eps
        first=float(np.linalg.cond(h))
        second=float(np.linalg.cond(h+perturb))
        reproduce.curvature_condition(first,second,'perturbed 3x3 Hessian')

    def test_approach_recovers_exact_synthetic_curve(self):
        rng=np.random.default_rng(16)
        d=rng.normal(size=(71,24))*.001
        t=np.arange(24,dtype=float)
        y=approach_model.model(t,[.2,.04,.25])
        candidates=approach_model.search(t,y,covariance(d),[.005,2])
        row,_=approach_model.assess(t,y,d,20,candidates,[.005,2])
        np.testing.assert_allclose([row['A'],row['B'],row['c']],[.2,.04,.25],rtol=3e-7)
        self.assertTrue(row['diagnostic_screen'])
        bad=copy.deepcopy(candidates);bad[0]['params'][0]+=.001
        with self.assertRaises(ValueError):approach_model.assess(t,y,d,20,bad,[.005,2])

    def test_flat_curve_does_not_resolve_decay_rate(self):
        rng=np.random.default_rng(17);d=rng.normal(size=(71,16))*.001
        t=np.arange(16,dtype=float);y=np.full(16,.2)
        c=approach_model.search(t,y,covariance(d),[.005,2])
        row,_=approach_model.assess(t,y,d,20,c,[.005,2])
        self.assertFalse(row['diagnostic_screen'])
        self.assertTrue(row['resolution_flags'])

    def test_python36_syntax_and_no_control_characters(self):
        for p in data.HERE.glob('*.py'):
            text=p.read_text()
            if sys.version_info >= (3,8):
                ast.parse(text,filename=str(p),feature_version=(3,6))
            else:
                ast.parse(text,filename=str(p))
            self.assertFalse(any(ord(c)<32 and c not in '\n' for c in text),p.name)


class PublishedDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=data.load()
        cls.tables,cls.summary=analysis.compute(cls.data)

    def test_full71_coverage_and_overlap(self):
        self.assertEqual(self.data['audit']['defined_energy_points'],1814)
        self.assertEqual(self.data['audit']['raw_overlap_records'],32589)
        self.assertEqual(self.data['audit']['covariance_entries_verified'],204800)
        self.assertEqual(len(self.tables['data_coverage']),68)
        for row in self.tables['data_coverage']:
            extended=row['R'] in data.LONG_RADII
            self.assertEqual(row['max_correlator_tau'],32 if extended else 16)
            self.assertEqual(row['max_effective_lower_tau'],31 if extended else 15)

    def test_all_undefined_and_weak_points_are_retained(self):
        energies=self.tables['effective_energies']
        self.assertEqual(len(energies),1856)
        bad=[r for r in energies if not r['valid_point']]
        self.assertEqual(len(bad),42)
        self.assertTrue(all(r['aE_eff']=='' and r['naive_se']=='' for r in bad))

    def test_published_fit_regression_inventory(self):
        self.assertEqual(self.summary['fit_regression'],dict(short_fits_verified=544,long_fits_verified=576))
        self.assertEqual(len(self.tables['constant_fit_catalogue']),2068)

    def test_missing_time_coverage_is_not_zero_or_plateau(self):
        row=next(r for r in self.tables['spatial_window_estimates']
                 if r['Nv']==100 and r['R']==3 and r['tau_min']==22)
        self.assertEqual(row['status'],'outside_measured_time')
        self.assertIsNone(row['estimate'])
        self.assertFalse(row['diagnostic_screen'])

    def test_R24_weak_Nv100_fit_remains_flagged(self):
        row=next(r for r in self.tables['constant_fit_catalogue']
                 if (r['Nv'],r['R'],r['tau_min'],r['tau_max'])==(100,24,14,19))
        self.assertEqual(row['status'],'defined')
        self.assertTrue(row['low_signal'])
        self.assertFalse(row['diagnostic_screen'])
        self.assertIsNone(row['jackknife_se'])

    def test_known_paired_cutoff_difference(self):
        row=next(r for r in self.tables['paired_cutoff_contrasts']
                 if (r['Nv'],r['R'],r['tau_min'],r['tau_max'])==(40,12,14,19))
        self.assertAlmostEqual(row['difference_over_se'],-3.39987911049,places=7)
        self.assertFalse(row['low_signal'])

    def test_each_raw_geometry_has_every_configuration(self):
        self.assertTrue(all(len(g)==71 for g in self.data['raw'].values()))
        self.assertEqual(data.CONFIGS[54],640)


if __name__=='__main__':unittest.main()
