#!/usr/bin/env python3
"""Independent numerical checks for the bounded, correlated fit machinery."""
import unittest
import numpy as np
import fit_models as fm


def synthetic(states=1):
    t=np.arange(15,dtype=float)
    p=np.array([1.1,.18]) if states==1 else np.array([.65,.45,.13,.34])
    rng=np.random.RandomState(91)
    z=rng.normal(size=(71,len(t)))
    z-=z.mean(0)
    corr=.65**np.abs(t[:,None]-t[None,:])
    noise=z@np.linalg.cholesky(corr).T*.001
    return fm.Window(fm.model(p,t,states)+noise,t,scale=1.),p


class TestCorrelatedFits(unittest.TestCase):
    def test_model_first_and_second_derivatives(self):
        for states,p in ((1,np.array([.8,.17])),(2,np.array([.6,.4,.12,.29]))):
            t=np.linspace(0,13,16);J,H=fm.derivatives(p,t,states)
            for j in range(len(p)):
                step=np.zeros(len(p));step[j]=1e-5
                numerical=(fm.model(p+step,t,states)-fm.model(p-step,t,states))/(2e-5)
                np.testing.assert_allclose(J[:,j],numerical,rtol=2e-7,atol=1e-9)
                numerical_H=(fm.derivatives(p+step,t,states)[0]-fm.derivatives(p-step,t,states)[0])/(2e-5)
                np.testing.assert_allclose(H[:,:,j],numerical_H,rtol=2e-7,atol=1e-8)

    def test_whitening_equals_full_quadratic_form(self):
        w,p=synthetic();p=p+np.array([.003,.001]);r=fm.model(p,w.t,1)-w.mean
        self.assertAlmostEqual(float(w.residual(p,1)@w.residual(p,1)),float(r@np.linalg.solve(w.covariance,r)),places=6)

    def test_single_and_two_state_recovery(self):
        for states in (1,2):
            w,p=synthetic(states)
            single,_=fm.solve(w,1)
            best,_=fm.solve(w,states,single=single['p'])
            self.assertIsNotNone(best)
            np.testing.assert_allclose(best['p'],p,rtol=1e-5,atol=1e-7)

    def test_curvature_propagation_matches_refitting(self):
        # A curved model with nonzero residual tests the residual-Hessian term.
        w,p=synthetic();mean=w.mean+.001*np.sin(w.t*.6)
        best,_=fm.solve(w,1,mean=mean);p=best['p']
        local=w.local_errors(p,1,mean)
        direction=np.cos(w.t*.4);step=1e-5
        plus,_=fm.solve(w,1,starts=[p],mean=mean+step*direction)
        minus,_=fm.solve(w,1,starts=[p],mean=mean-step*direction)
        numerical=(plus['p']-minus['p'])/(2*step)
        np.testing.assert_allclose(local['gradient']@direction,numerical,rtol=3e-4,atol=2e-5)
        np.testing.assert_allclose(np.cov(local['fluctuations'],rowvar=False,ddof=1)/w.n,local['covariance'],rtol=1e-9,atol=1e-15)

    def test_rescaling_leaves_energy_and_objective_invariant(self):
        w,p=synthetic();raw=w.samples.copy()
        a=fm.Window(raw,w.times);b=fm.Window(raw*1e-20,w.times)
        pa,_=fm.solve(a,1);pb,_=fm.solve(b,1)
        np.testing.assert_allclose(pa['p'],pb['p'],rtol=1e-8,atol=1e-10)
        self.assertAlmostEqual(pa['q2'],pb['q2'],places=6)

    def test_signed_measurements_are_retained(self):
        rng=np.random.RandomState(123);t=np.arange(6,dtype=float)
        raw=np.exp(-.3*t)+rng.normal(size=(71,6))*1.5
        self.assertTrue(np.any(raw<0));w=fm.Window(raw,t,scale=1.)
        np.testing.assert_array_equal(raw,w.samples)
        best,_=fm.solve(w,1);self.assertIsNotNone(best)

    def test_singular_covariance_is_rejected(self):
        raw=np.repeat(np.arange(71,dtype=float)[:,None],5,axis=1)
        with self.assertRaises((ValueError,np.linalg.LinAlgError)):fm.Window(raw,np.arange(5))

    def test_two_state_boundary_is_explicit(self):
        w,p=synthetic();q=np.array([p[0],0.,p[1],.3])
        d=fm.diagnostics(w,2,q,float(w.residual(q,2)@w.residual(q,2)))
        self.assertFalse(d['identified']);self.assertIn('parameter_boundary',d['flags'])
        self.assertIn('vanishing_state_amplitude',d['flags'])

    def test_recentered_bootstrap_null_and_covariance(self):
        w,p=synthetic();p=p+np.array([.01,.003]);pred=fm.model(p,w.t,1)
        null=w.samples-w.mean+pred
        np.testing.assert_allclose(null.mean(0),pred,rtol=1e-13,atol=1e-15)
        chosen=np.random.RandomState(7).randint(0,71,size=71)
        draw=fm.Window(null[chosen],w.times,scale=1.)
        np.testing.assert_allclose(draw.covariance,np.cov(null[chosen],rowvar=False,ddof=1)/71)
        self.assertFalse(np.allclose(draw.covariance,w.covariance,rtol=1e-6,atol=1e-14))


if __name__=='__main__':unittest.main()
