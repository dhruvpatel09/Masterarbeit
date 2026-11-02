#!/usr/bin/env python3
import unittest
import numpy as np
import joint_models as jm


def synthetic(kind='shared_one',m=2):
    t=np.arange(0,17,2,dtype=float);nt=len(t)
    p=np.r_[np.linspace(.8,1.2,m),.18] if kind=='shared_one' else np.r_[np.linspace(.8,1.2,m),np.linspace(.15,.22,m)] if kind=='independent_one' else np.r_[np.linspace(.5,.7,m),np.linspace(.3,.6,m),.13,.33]
    rng=np.random.RandomState(901);z=rng.normal(size=(71,m*nt));z-=z.mean(0)
    c=.45**np.abs(np.arange(m*nt)[:,None]-np.arange(m*nt)[None,:])
    noise=z@np.linalg.cholesky(c).T*.001
    raw=(jm.model(p,t,m,kind)+noise).reshape(71,m,nt)
    return jm.JointWindow(raw,t,list(range(m)),scales=np.ones(m)),p


class TestJointModels(unittest.TestCase):
    def test_analytic_derivatives(self):
        for kind in jm.KINDS:
            w,p=synthetic(kind,3);J,H=jm.derivatives(p,w.dt,w.m,kind)
            for j in range(len(p)):
                delta=np.zeros(len(p));delta[j]=1e-5
                fd=(jm.model(p+delta,w.dt,w.m,kind)-jm.model(p-delta,w.dt,w.m,kind))/(2e-5)
                np.testing.assert_allclose(J[:,j],fd,rtol=3e-7,atol=1e-9)
                fdH=(jm.derivatives(p+delta,w.dt,w.m,kind)[0]-jm.derivatives(p-delta,w.dt,w.m,kind)[0])/(2e-5)
                np.testing.assert_allclose(H[:,:,j],fdH,rtol=3e-7,atol=2e-8)

    def test_cross_source_covariance_is_retained(self):
        w,p=synthetic();v=jm.model(p+.001,w.dt,w.m,'shared_one')-w.mean
        self.assertGreater(np.max(np.abs(w.covariance[:w.nt,w.nt:])),0)
        self.assertAlmostEqual(float(v@np.linalg.solve(w.covariance,v)),float((w.W@v)@(w.W@v)),places=6)

    def test_exact_shared_energy_recovery(self):
        for m in (2,4):
            w,p=synthetic(m=m);best,_=jm.solve(w,'shared_one');np.testing.assert_allclose(best['p'],p,rtol=1e-6,atol=1e-8)

    def test_exact_independent_and_two_state_recovery(self):
        for kind in ('independent_one','shared_two'):
            w,p=synthetic(kind);one,_=jm.solve(w,'shared_one');ind,_=jm.solve(w,'independent_one',shared=one)
            best=ind if kind=='independent_one' else jm.solve(w,'shared_two',shared=one,independent=ind)[0]
            np.testing.assert_allclose(best['p'],p,rtol=2e-5,atol=2e-7)

    def test_implicit_gradient_against_refitting(self):
        w,p=synthetic();mean=w.mean+.001*np.sin(w.t*.5);best,_=jm.solve(w,'shared_one',mean=mean);p=np.array(best['p'])
        local=w.local_errors(p,'shared_one',mean);direction=np.cos(w.t*.2);eps=1e-5
        plus,_=jm.solve(w,'shared_one',starts=[p],mean=mean+eps*direction)
        minus,_=jm.solve(w,'shared_one',starts=[p],mean=mean-eps*direction)
        fd=(np.array(plus['p'])-np.array(minus['p']))/(2*eps)
        np.testing.assert_allclose(local['gradient']@direction,fd,rtol=5e-4,atol=3e-5)

    def test_source_normalization_and_order_invariance(self):
        w,p=synthetic();raw=w.samples.reshape(71,w.m,w.nt)
        other=jm.JointWindow(raw[:,::-1]*np.array([1e-12,3.])[None,:,None],w.source_times,w.cutoffs[::-1])
        first,_=jm.solve(w,'shared_one');second,_=jm.solve(other,'shared_one')
        self.assertAlmostEqual(first['p'][-1],second['p'][-1],places=7)
        self.assertAlmostEqual(first['q2'],second['q2'],places=6)

    def test_nested_model_point_and_objective(self):
        w,p=synthetic('independent_one');one,_=jm.solve(w,'shared_one');ind,_=jm.solve(w,'independent_one',shared=one)
        self.assertLessEqual(ind['q2'],one['q2']+1e-6)
        nested=np.r_[one['p'][:w.m],np.repeat(one['p'][w.m],w.m)]
        np.testing.assert_allclose(jm.model(nested,w.dt,w.m,'independent_one'),jm.model(np.array(one['p']),w.dt,w.m,'shared_one'))

    def test_signed_values_retained_and_singular_sources_rejected(self):
        rng=np.random.RandomState(6);raw=rng.normal(size=(71,2,5))+np.exp(-.3*np.arange(5))
        w=jm.JointWindow(raw,np.arange(5),[40,100],scales=np.ones(2));np.testing.assert_array_equal(w.samples,raw.reshape(71,-1))
        raw[:,1]=raw[:,0]
        with self.assertRaises((ValueError,np.linalg.LinAlgError)):jm.JointWindow(raw,np.arange(5),[40,100])

    def test_null_recentering_preserves_matching(self):
        w,p=synthetic();null=w.samples-w.mean+jm.model(p,w.dt,w.m,'shared_one')
        np.testing.assert_allclose(null.mean(0),jm.model(p,w.dt,w.m,'shared_one'))
        idx=np.random.RandomState(2).randint(0,71,71);sample=jm.JointWindow(null[idx].reshape(71,w.m,w.nt),w.source_times,w.cutoffs,scales=np.ones(w.m))
        np.testing.assert_allclose(sample.covariance,np.cov(null[idx],rowvar=False,ddof=1)/71)

    def test_configuration_axis_preserved_in_time_source_selection(self):
        import analysis
        a=np.arange(71*4*33*12).reshape(71,4,33,12)
        for case in analysis.cases():
            actual=analysis.selected_samples(case,a)
            expected=np.array([[[a[i,analysis.CUTOFFS.index(n),t,analysis.RADII.index(case['R'])]
                for t in case['times']] for n in case['cutoffs']] for i in range(71)])
            np.testing.assert_array_equal(actual,expected)

    def test_zero_exceedances_retains_finite_draw_uncertainty(self):
        from analysis import tail_summary
        r=tail_summary([0.]*399,1.)
        self.assertEqual(r['exceedances'],0);self.assertEqual(r['p_add_one'],1/400)
        self.assertEqual(r['conditional_tail_low95'],0.)
        self.assertGreater(r['conditional_tail_high95'],.009)

    def test_paired_ratio_cancels_common_fluctuations_and_keeps_sign(self):
        from analysis import ratio_estimate
        reference=np.linspace(1,2,71);value,se=ratio_estimate(2*reference,reference)
        self.assertEqual(value,2.);self.assertEqual(se,0.)
        value,se=ratio_estimate(-reference,reference);self.assertEqual(value,-1.);self.assertEqual(se,0.)
        with self.assertRaises(ValueError):ratio_estimate(reference,np.zeros(71))


class TestSummaryReplay(unittest.TestCase):
    def test_summary_accepts_float_roundoff(self):
        import copy
        from reproduce import check_json
        expected={'cohort':[100,110], 'rows':[{'draws':399,'p':.0325,'optional':None,'success':True}]}
        actual=copy.deepcopy(expected)
        actual['rows'][0]['p']=float(np.nextafter(.0325,np.inf))
        self.assertNotEqual(actual,expected)
        check_json(actual,expected)

    def test_summary_rejects_changed_counts_structure_and_values(self):
        import copy
        from reproduce import check_json
        expected=[{'case_id':'pair_dense_R12_t12_32','draws':399,'p':.0125}]
        changes=[('draws',398),('draws',399.),('case_id','another_case'),('p',.02)]
        for field,value in changes:
            actual=copy.deepcopy(expected);actual[0][field]=value
            with self.assertRaises(ValueError):check_json(actual,expected)
        for actual in ([],expected+expected,[{'draws':399,'p':.0125}]):
            with self.assertRaises(ValueError):check_json(actual,expected)
        with self.assertRaises(ValueError):check_json([110,100],[100,110])

    def test_summary_rejects_nonfinite_and_boolean_numbers(self):
        from reproduce import check_json
        for value in (float('nan'),float('inf'),-float('inf'),True):
            with self.assertRaises(ValueError):check_json({'p':value},{'p':.5})


if __name__=='__main__':unittest.main()
