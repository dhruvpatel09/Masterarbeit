#!/usr/bin/env python3
"""Independent statistical controls and the captured native cutoff kernels."""
import ast
import importlib.util
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from decimal import Decimal as D,localcontext
from pathlib import Path

import numpy as np
import estimators as s
import reproduce as a


class Tests(unittest.TestCase):
    def test_matched_normalization_and_log_cancel_common_noise(self):
        x=list(map(D,['1','2','4','7']))
        q,e=a.normalized(x,x);self.assertEqual((q,e),(D(1),D(0)))
        result,fluct,loo=s.energy([2*v for v in x],x)
        self.assertEqual(result['naive_se'],0)
        # Averaging identical 50-digit logarithms can leave last-digit residue.
        self.assertLess(result['jackknife_se'],D('1e-48'))
        self.assertEqual(len(loo),4)
        xs=np.array([2.,1.,3.,4.]);ys=np.array([1.,2.,1.,3.])
        ratio,error=a.normalized(list(map(D,map(str,xs))),list(map(D,map(str,ys))))
        gradient=np.array([1/ys.mean(),-xs.mean()/ys.mean()**2])
        covariance=np.cov(xs,ys,ddof=1)/len(xs)
        self.assertAlmostEqual(float(error),float(np.sqrt(gradient@covariance@gradient)),places=14)

    def test_signed_samples_and_invalid_log_jackknife_are_retained(self):
        x=list(map(D,['-0.9','0.1','1.1']));y=[D(1)]*3
        result,fluct,loo=s.energy(x,y)
        self.assertTrue(result['valid_point']);self.assertFalse(result['valid_jackknife'])
        self.assertEqual(loo,[]);self.assertEqual(a.count_positive_loo(x,y),2)
        self.assertEqual(len(fluct),3)
        result,_,_=s.energy([D(-2),D(-1),D(1)],y)
        self.assertFalse(result['valid_point'])

    def test_correlated_fit_matches_closed_form_weights(self):
        fluctuations=np.array([[1.,2.],[-1.,-2.],[2.,-1.],[-2.,1.],[.5,1.],[-.5,-1.]])
        C=np.cov(fluctuations,rowvar=False,ddof=1)/len(fluctuations)
        denominator=C[0,0]+C[1,1]-2*C[0,1]
        weights=np.array([C[1,1]-C[0,1],C[0,0]-C[0,1]])/denominator
        fit=s.gls([2.,3.],fluctuations,np.ones((2,1)))
        np.testing.assert_allclose(fit['weights'][0],weights,rtol=1e-12,atol=1e-14)
        self.assertAlmostEqual(fit['params'][0],float(weights@np.array([2.,3.])))

    def test_correlated_fit_rejects_singular_covariance(self):
        with self.assertRaises((ValueError,np.linalg.LinAlgError)):
            s.gls([1.,2.],np.array([[1.,1.],[-1.,-1.],[0.,0.]]),np.ones((2,1)))

    def test_covariance_matches_decimal_reference(self):
        x=[D('0.01'),D('0.07'),D('-0.02'),D('0.04')]
        y=[D('0.2'),D('-0.1'),D('0.5'),D('0.3')]
        expected=float(s.covariance(x,y))
        self.assertAlmostEqual(np.cov(np.array(x,float),np.array(y,float),ddof=1)[0,1]/4,expected,places=16)

    def test_python36_core_syntax(self):
        for path in a.HERE.glob('*.py'):
            options=dict(feature_version=(3,6)) if sys.version_info>=(3,8) else {}
            ast.parse(path.read_text(),filename=str(path),**options)

    @unittest.skipUnless(importlib.util.find_spec('scipy'),'optional SciPy chi-square tail control')
    def test_chi_square_reference_tail(self):
        from scipy.stats import chi2
        for df in range(1,10):
            for q in (0.,.1,2.,10.,100.):
                self.assertAlmostEqual(s.chi2_sf(q,df),chi2.sf(q,df),places=13)

    @unittest.skipUnless(shutil.which('gcc'),'GCC is needed for the captured native kernel control')
    def test_captured_native_contractions_and_periodic_payloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'kernels'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(a.SOURCE/'tests'),
                '-I'+str(a.SOURCE/'build_provenance'),str(a.SOURCE/'tests/kernels_and_shifts.c'),
                '-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)


if __name__=='__main__':unittest.main()
