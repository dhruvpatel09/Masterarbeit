#!/usr/bin/env python3
"""Controls for nonlinear errors, correlated fits and spatial covariance."""
import math
import unittest
from decimal import Decimal as D,localcontext
import numpy as np
import estimators as e
import reproduce as a


class Tests(unittest.TestCase):
    def test_sample_mean_error(self):
        m,s=e.mean_se(list(map(D,(1,2,3))))
        self.assertEqual(m,2)
        with localcontext() as c:
            c.prec=50;self.assertLess(abs(s*s-D(1)/3),D('1e-49'))

    def test_correlated_log_amplitude_cancels(self):
        r,d,j=e.energy(list(map(D,(2,4,8))),list(map(D,(1,2,4))))
        self.assertTrue(r['valid_jackknife']);self.assertLess(r['naive_se'],D('1e-48'));self.assertLess(r['jackknife_se'],D('1e-48'))

    def test_log_ratio_of_means(self):
        r,_,_=e.energy(list(map(D,(1,2,9))),list(map(D,(1,2,3))))
        with localcontext() as c:c.prec=50;self.assertEqual(r['aE_eff'],D(2).ln())

    def test_nonpositive_means_and_leave_one_out(self):
        r,_,_=e.energy(list(map(D,(-4,1,1))),[D(1)]*3)
        self.assertFalse(r['valid_point']);self.assertEqual(r['aE_eff'],'')
        r,_,_=e.energy(list(map(D,(-4,3,3))),[D(1)]*3)
        self.assertTrue(r['valid_point']);self.assertFalse(r['valid_jackknife'])

    def test_mismatched_samples(self):
        with self.assertRaises(ValueError):e.energy([D(1)]*3,[D(1)]*2)

    def test_chi_square_known_tails(self):
        for q in (.01,1,10,100):
            self.assertAlmostEqual(e.chi2_sf(q,2),math.exp(-q/2),places=14)
            self.assertAlmostEqual(e.chi2_sf(q,1),math.erfc(math.sqrt(q/2)),places=14)
            self.assertGreaterEqual(e.chi2_sf(q,6),e.chi2_sf(q,4))

    def test_chi_square_against_scipy(self):
        try:from scipy.stats import chi2
        except ImportError:self.skipTest('Optional SciPy reference')
        for df in range(1,20):
            for q in (.01,1,10,100,500,1500):
                value=e.chi2_sf(q,df);reference=chi2.sf(q,df)
                if reference==0:self.assertLess(value,1e-300)  # SciPy can underflow before subnormal tails.
                else:self.assertTrue(math.isclose(value,reference,rel_tol=3e-12,abs_tol=1e-320))

    def test_correlated_cornell_fit_and_units(self):
        R=np.arange(4,11,dtype=float);A=np.column_stack((np.ones(len(R)),R,-1/R));true=np.array((.15,.005,.24))
        i=np.arange(1,72)[:,None];j=np.arange(1,8)[None,:]
        d=.01*np.sin(i*j)+.02*np.cos(i/3)
        first=e.gls(A@true,d,A);second=e.gls(1000*A@true,1000*d,A)
        np.testing.assert_allclose(first['params'],true,rtol=1e-11,atol=1e-12)
        np.testing.assert_allclose(second['params'],1000*true,rtol=1e-11,atol=1e-9)
        self.assertLess(first['chi2'],1e-20)
        np.testing.assert_allclose(second['covariance'],1e6*first['covariance'],rtol=1e-10)

    def test_singular_covariance_rejected(self):
        with self.assertRaises(ValueError):e.gls(np.ones(4),np.ones((71,4)),np.ones((4,1)))

    def test_common_spatial_fluctuations_are_not_independent(self):
        blocks=[]
        for axis in (0,1,2):
            for index,cfg in enumerate(a.CONFIGS):
                for x in range(4):
                    for y in range(4):
                        for z in range(4):blocks.append(dict(tau=4,R=1,axis=axis,cfg=cfg,px=x,py=y,pz=z,L_RE=str(1+index*.01+x*.001)))
        rows=a.block_diagnostics(blocks)
        for r in rows:self.assertAlmostEqual(r['correct_over_ignore_regional_covariance'],math.sqrt(r['regional_channels']),places=10)


if __name__=='__main__':unittest.main()
