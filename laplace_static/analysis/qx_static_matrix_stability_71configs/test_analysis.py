#!/usr/bin/env python3
"""Focused tests of complex null tests, holdout identity and nonlinear fits."""
import ast,unittest
import numpy as np
from common import HERE,read,load_data,jk_cov,jk_se,previous_projections,close,check_equal_time_reality
from source_splits import make_assignments
from complex_checks import antihermitian_components,max_t_bootstrap
from correlator_fits import ExponentialFit,whitening,SingularCovariance
from analysis import constant_fit,curve_tables


class ComplexTests(unittest.TestCase):
    def test_complex_hermitian_matrix_has_zero_antihermitian_part(self):
        a=np.arange(16).reshape(4,4)+1j*np.arange(16)[::-1].reshape(4,4)
        a=a+a.conj().T
        names,x=antihermitian_components(a[None,:,:])
        self.assertEqual(len(names),16);np.testing.assert_array_equal(x,0)
        self.assertGreater(np.max(abs(a.imag)),0)

    def test_antisymmetric_real_matrix_is_not_hermitian(self):
        a=np.zeros((1,4,4),complex);a[0,0,1]=2;a[0,1,0]=-2
        names,x=antihermitian_components(a)
        self.assertEqual(x[0,names.index(('offdiag_re',0,1))],2)

    def test_shared_vector_bootstrap_preserves_duplicate_components(self):
        rng=np.random.RandomState(4);x=rng.normal(size=(71,2));x=np.column_stack([x,x[:,0]])
        draws=rng.randint(0,71,size=(199,71))
        result=max_t_bootstrap(x,draws)
        self.assertEqual(result[3][0],result[3][2])
        self.assertTrue(np.all(result[3]>=1/200.))

    def test_large_nonzero_mean_is_detected_by_joint_null(self):
        rng=np.random.RandomState(8);x=rng.normal(size=(71,3))+10
        draws=rng.randint(0,71,size=(199,71));result=max_t_bootstrap(x,draws)
        self.assertEqual(result[5],1/200.)

    def test_zero_variance_is_not_silently_regularized(self):
        with self.assertRaises(ValueError):max_t_bootstrap(np.ones((71,2)),np.tile(np.arange(71),(10,1)))

    def test_degenerate_draw_is_not_discarded(self):
        x=np.arange(71.)[:,None]
        with self.assertRaises(ValueError):max_t_bootstrap(x,np.zeros((1,71),int))


class FitTests(unittest.TestCase):
    def setUp(self):
        self.t=np.arange(7);self.cov=.001*.5**abs(self.t[:,None]-self.t[None,:])

    def test_exact_exponential_recovery(self):
        fit=ExponentialFit(self.cov,self.t).fit_many(2*np.exp(-.22*self.t))[0]
        self.assertAlmostEqual(fit['energy'],.22,places=11)
        self.assertAlmostEqual(fit['amplitude'],2,places=11)
        self.assertLess(fit['q'],1e-20)

    def test_signed_amplitude_and_negative_energy_are_allowed(self):
        fit=ExponentialFit(self.cov,self.t).fit_many(-2*np.exp(.13*self.t))[0]
        self.assertAlmostEqual(fit['energy'],-.13,places=11)
        self.assertAlmostEqual(fit['amplitude'],-2,places=11)

    def test_global_rescaling_preserves_energy_and_chi_square(self):
        y=2*np.exp(-.22*self.t)+np.array([.01,-.02,.01,0,.02,-.01,0])
        a=ExponentialFit(self.cov,self.t).fit_many(y)[0]
        b=ExponentialFit(self.cov*1e-14,self.t).fit_many(y*1e-7)[0]
        self.assertAlmostEqual(a['energy'],b['energy'],places=11)
        self.assertAlmostEqual(a['q'],b['q'],places=9)
        self.assertAlmostEqual(a['amplitude']*1e-7,b['amplitude'],places=15)

    def test_cholesky_whitening(self):
        w,c=whitening(self.cov);np.testing.assert_allclose(w.dot(self.cov).dot(w.T),np.eye(7),atol=1e-13)

    def test_singular_covariance_is_rejected(self):
        for cov in (np.ones((7,7)),np.zeros((7,7))):
            with self.assertRaises(SingularCovariance):ExponentialFit(cov,self.t)

    def test_boundary_solution_is_flagged(self):
        fit=ExponentialFit(self.cov,self.t).fit_many(np.exp(-3*self.t))[0]
        self.assertTrue(fit['boundary']);self.assertEqual(fit['energy'],1.5)

    def test_profile_derivatives_match_finite_difference(self):
        obj=ExponentialFit(self.cov,self.t);y=obj.W.dot(2*np.exp(-.22*self.t));E=.3;h=1e-5
        q,A,g,H=obj.evaluate(E,y);left=obj.evaluate(E-h,y);right=obj.evaluate(E+h,y)
        self.assertAlmostEqual(g,(right[0]-left[0])/(2*h),delta=1e-4)
        self.assertAlmostEqual(H,(right[2]-left[2])/(2*h),delta=.01)

    def test_finer_grid_finds_the_same_minimum(self):
        y=2*np.exp(-.22*self.t)+np.sin(self.t)*.01
        a=ExponentialFit(self.cov,self.t,grid_points=161).fit_many(y)[0]
        b=ExponentialFit(self.cov,self.t,grid_points=321).fit_many(y)[0]
        self.assertAlmostEqual(a['energy'],b['energy'],places=11)

    def test_nonfinite_fit_data_are_rejected(self):
        y=np.ones(7);y[3]=np.nan
        with self.assertRaises(ValueError):ExponentialFit(self.cov,self.t).fit_many(y)


class EqualTimeRealityTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.RandomState(73)
        self.means=np.exp(-.1*np.arange(33))[None,:].astype(complex)
        self.deletions=np.tile(self.means[:,None,:],(1,71,1))
        self.deletions+=rng.normal(size=(1,71,33))*.0001
        self.means[0,0]+=5.253035877957937e-19j
        self.deletions[0,:,0]+=rng.normal(size=71)*1e-18j
        self.means[0,1]+=2e-5j
        self.deletions[0,:,1]+=(2e-5+rng.normal(size=71)*1e-6)*1j

    def table(self,variant='alternating'):
        return curve_tables([dict(R=4,basis='pair',anchor=4,variant=variant)],self.means,self.deletions)[0]

    def test_every_source_variant_excludes_equal_time_imaginary_studentization(self):
        for variant in ('same_sample','alternating','seed1701','seed1702','seed1703','seed1704'):
            a=self.table(variant)
            self.assertEqual(a[0]['imaginary_pointwise_t'],'')
            if variant!='same_sample':
                self.assertEqual(a[0]['diagnostic'],'roundoff_imaginary_equal_time')
                self.assertNotEqual(a[0]['real_mean_over_se'],'')

    def test_raw_equal_time_imaginary_values_are_retained(self):
        a=self.table()
        self.assertEqual(a[0]['mean_im'],self.means[0,0].imag)
        np.testing.assert_allclose(a[0]['jackknife_se_im'],float(jk_se(self.deletions[0,:,0].imag)),rtol=8*np.finfo(float).eps,atol=0)

    def test_nonzero_time_imaginary_significance_is_still_computed(self):
        a=self.table()
        self.assertEqual(a[1]['imaginary_pointwise_t'],float(self.means[0,1].imag/jk_se(self.deletions[0,:,1].imag)))
        self.assertEqual(a[1]['diagnostic'],'physical_time_diagnostic')

    def test_substantive_central_equal_time_imaginary_part_is_rejected(self):
        self.means[0,0]+=1e-10j
        with self.assertRaisesRegex(ValueError,'roundoff bound'):self.table()

    def test_substantive_deleted_equal_time_imaginary_part_is_rejected(self):
        self.deletions[0,30,0]+=1e-10j
        with self.assertRaisesRegex(ValueError,'roundoff bound'):self.table()

    def test_roundoff_ratio_changes_do_not_relax_physical_comparisons(self):
        a=self.table()
        self.means[0,0]=self.means[0,0].real-3.542745126994888e-19j
        b=self.table()
        self.assertEqual(a[0]['imaginary_pointwise_t'],b[0]['imaginary_pointwise_t'])
        self.assertNotEqual(a[0]['mean_im'],b[0]['mean_im'])
        with self.assertRaisesRegex(ValueError,'archived=.*replayed='):
            close([.2],[.20001],'physical energy still strict')


class SamplingTests(unittest.TestCase):
    def test_balanced_folds_are_reproducible_and_complete(self):
        plan=read(HERE/'ANALYSIS_PLAN.json');a=make_assignments(plan)
        np.testing.assert_array_equal(a[0],np.arange(71)%2)
        np.testing.assert_array_equal(a,make_assignments(plan))
        for fold in a:self.assertEqual(sorted(np.bincount(fold).tolist()),[35,36])
        self.assertEqual(len({tuple(row) for row in a}),5)

    def test_deleted_configuration_does_not_reindex_fold_membership(self):
        a=make_assignments(read(HERE/'ANALYSIS_PLAN.json'))[0];keep=np.arange(71)!=4
        self.assertFalse(np.array_equal(a[keep],np.arange(70)%2))

    def test_jackknife_linear_mean_and_paired_cancellation(self):
        x=np.arange(71.)**2;r=(x.sum()-x)/70
        self.assertAlmostEqual(float(jk_se(r)),float(x.std(ddof=1)/np.sqrt(71)),places=10)
        np.testing.assert_array_equal(jk_se(np.column_stack([r,r])-np.column_stack([r,r])),[0,0])

    def test_normalization_point_is_not_reported_as_physical_signal_precision(self):
        means=np.exp(-.1*np.arange(33))[None,:].astype(complex)
        rng=np.random.RandomState(9);d=np.tile(means[:,None,:],(1,71,1));d[:,:,1:]+=rng.normal(size=(1,71,32))*.001
        tables,_,_,_=curve_tables([dict(R=4,basis='pair',anchor=4,variant='same_sample')],means,d)
        self.assertEqual(tables[0]['real_mean_over_se'],'')

    def test_published_input_shape(self):
        self.assertEqual(load_data().shape,(71,4,33,4,4))

    def test_undefined_masks_and_infinities_are_checked(self):
        close([1,np.nan],[1,np.nan],'valid mask')
        for other in ([1,2],[1,np.inf]):
            with self.assertRaises(ValueError):close([1,np.nan],other,'invalid mask')

    def test_python36_syntax(self):
        for p in HERE.glob('*.py'):
            try:ast.parse(p.read_text(),filename=p.name,feature_version=(3,6))
            except TypeError:ast.parse(p.read_text(),filename=p.name)


if __name__=='__main__':unittest.main()
