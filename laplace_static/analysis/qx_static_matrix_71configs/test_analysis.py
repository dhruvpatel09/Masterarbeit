#!/usr/bin/env python3
"""Focused checks of orientation, source training and paired uncertainty."""
import ast
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from import_review import HERE,check_file,import_archive
from replay import safe_name,compare_case_summary
from statistics_tools import metric_transform,jk_cov,jk_se,energies
from projections import train,quadratic,estimate
from analysis import baseline
from reproduce import normalization_signal


def synthetic(n=10):
    rng=np.random.RandomState(1904)
    data=[]
    for i in range(n):
        u=np.array([[1+.05*rng.normal(),.2+.02*rng.normal()],[.1j,1+.05*rng.normal()]])
        C=[]
        for t in range(10):
            # Two positive states plus configuration-dependent, nonsymmetric
            # imaginary fluctuations at nonzero time.
            a=u.dot(np.diag(np.exp(-t*np.array([.18+.01*rng.normal(),.5])))).dot(u.conj().T)
            if t:a[0,1]+=.0001j*rng.normal()
            C.append(a)
        data.append(C)
    return np.array(data)


class SourceTests(unittest.TestCase):
    def test_whitening_and_rescaling(self):
        C=synthetic().mean(axis=0);w,e,error=metric_transform(C[0])
        np.testing.assert_allclose(w.conj().T.dot(C[0]).dot(w),np.eye(2),atol=1e-12)
        v,_=train(C,4);scale=np.diag([4.,.03]);other=np.einsum('ai,tab,bj->tij',scale,C,scale)
        u,_=train(other,4)
        np.testing.assert_allclose(quadratic(C,v),quadratic(other,u),atol=1e-11)

    def test_unit_source_recovers_diagonal(self):
        C=synthetic();v=np.array([0.,1.])
        np.testing.assert_array_equal(quadratic(C,v),C[...,1,1])

    def test_projection_global_phase_invariance(self):
        C=synthetic();v,_=train(C.mean(axis=0),4)
        np.testing.assert_allclose(quadratic(C,v),quadratic(C,v*np.exp(1.39j)),atol=1e-13)

    def test_preserves_complex_orientation(self):
        C=synthetic();v=np.array([1,1j])/np.sqrt(2)
        self.assertGreater(np.max(abs(quadratic(C,v)-quadratic(C.swapaxes(-1,-2),v))),.01)

    def test_nonsymmetric_is_not_silently_symmetrized(self):
        C=synthetic();C[:,:,0,1]+=.1j
        v=np.array([1,2])/np.sqrt(5)
        original=quadratic(C,v);h=quadratic((C+C.swapaxes(-1,-2).conj())/2,v)
        self.assertGreater(np.max(abs(original.imag-h.imag)),.01)

    def test_holdout_does_not_train_on_itself(self):
        C=synthetic();ids=np.arange(10)
        _,before=estimate(C,ids,4,'cross_fit',True)
        changed=C.copy();changed[0]*=10
        _,after=estimate(changed,ids,4,'cross_fit',True)
        # Original index 0 is in test fold 0. Its source uses odd indices only.
        np.testing.assert_array_equal(before[0][3],after[0][3])
        self.assertGreater(np.linalg.norm(before[1][3]-after[1][3]),1e-6)

    def test_crossfit_mean_uses_all_heldout_configurations(self):
        C=synthetic(11);ids=np.arange(11)
        mean,records=estimate(C,ids,4,'cross_fit',True)
        expected=sum(quadratic(C[ids%2==fold],rec[3]).sum(axis=0) for fold,rec in enumerate(records))/11
        np.testing.assert_allclose(mean,expected,atol=1e-15)

    def test_complete_retraining_changes_uncertainty(self):
        C=synthetic();ids=np.arange(10);v,_=train(C.mean(axis=0),4)
        fixed=[];refit=[]
        for k in ids:
            keep=ids!=k;fixed.append(quadratic(C[keep],v).mean(axis=0))
            refit.append(estimate(C[keep],ids[keep],4,'same_sample'))
        self.assertGreater(np.max(abs(jk_se(np.real(fixed))-jk_se(np.real(refit)))),1e-6)

    def test_metric_rejects_singular_and_nonhermitian(self):
        for M in (np.ones((2,2)),np.array([[1,1j],[0,1]])):
            with self.assertRaises(ValueError):metric_transform(M)


class StatisticsTests(unittest.TestCase):
    def test_normalization_signal_roundoff(self):
        a=dict(procedure='same_sample',tau=0,mean_re=1.,mean_im=0.,
               jackknife_se_re=2e-15,jackknife_se_im=0.,real_mean_over_se=5e14)
        b=dict(a,jackknife_se_re=4e-15,real_mean_over_se=2.5e14)
        normalization_signal(a,b,'normalization fixture')
        normalization_signal(a,dict(a,jackknife_se_re=0.,real_mean_over_se=float('inf')),'zero-error fixture')

    def test_normalization_signal_does_not_accept_physical_change(self):
        a=dict(procedure='same_sample',tau=0,mean_re=1.,mean_im=0.,
               jackknife_se_re=2e-15,jackknife_se_im=0.,real_mean_over_se=5e14)
        for b in (dict(a,mean_re=1.001),dict(a,jackknife_se_re=1e-8,real_mean_over_se=1e8),
                  dict(a,real_mean_over_se=1e8),dict(a,mean_im=1e-5),
                  dict(a,jackknife_se_re=float('nan'))):
            with self.subTest(changed=b):
                with self.assertRaises(ValueError):normalization_signal(a,b,'changed fixture')

    def test_normalization_signal_is_restricted_to_training_normalization(self):
        a=dict(procedure='same_sample',tau=0,mean_re=1.,mean_im=0.,
               jackknife_se_re=2e-15,jackknife_se_im=0.,real_mean_over_se=5e14)
        for b in (dict(a,tau=1),dict(a,procedure='cross_fit')):
            with self.assertRaises(ValueError):normalization_signal(a,b,'wrong point')

    def test_jackknife_mean_equals_independent_se(self):
        x=np.arange(10.)**2;loo=(x.sum()-x)/9
        self.assertAlmostEqual(float(jk_se(loo)),float(x.std(ddof=1)/np.sqrt(10)),places=12)

    def test_paired_covariance_retains_cancellation(self):
        x=np.arange(10.)**2; y=x+.01*np.arange(10.)
        paired=float(jk_se(x-y));uncorrelated=float(np.hypot(jk_se(x),jk_se(y)))
        self.assertLess(paired,uncorrelated/100)
        np.testing.assert_allclose(np.diag(jk_cov(np.column_stack([x,y]))),[jk_se(x)**2,jk_se(y)**2])

    def test_log_ratio_after_ensemble_average(self):
        x=np.array([[2.,1.],[1.,2.]])
        self.assertEqual(float(energies(x.mean(axis=0))[0]),0)

    def test_nonpositive_resample_is_not_dropped(self):
        E=energies(np.array([[2,1,.1],[2,0,-.1]],float))
        self.assertTrue(np.isfinite(E[0]).all());self.assertFalse(np.isfinite(E[1]).any())

    def test_signed_individual_values_are_retained(self):
        x=np.array([[3.,2.],[-1.,-.5],[2.,1.]])
        self.assertTrue(np.isfinite(energies(x.mean(axis=0))).all())


class ReceiptTests(unittest.TestCase):
    def test_native_pivot_roundoff_is_portable(self):
        recorded=dict(configuration=110,scans=[dict(piter=10,rows=8448,
            minimum_equal_time_relative_pivot=0.006686613774462862)])
        fresh=json.loads(json.dumps(recorded))
        fresh['scans'][0]['minimum_equal_time_relative_pivot']=float(np.nextafter(
            recorded['scans'][0]['minimum_equal_time_relative_pivot'],np.inf))
        self.assertNotEqual(recorded,fresh)  # Former whole-dict equality failed.
        compare_case_summary(recorded,fresh,110)

    def test_native_pivot_material_change_is_rejected(self):
        for value in (.0067,-.006686613774462862,0.):
            with self.subTest(value=value):
                old=dict(scans=[dict(minimum_equal_time_relative_pivot=.006686613774462862)])
                new=dict(scans=[dict(minimum_equal_time_relative_pivot=value)])
                with self.assertRaisesRegex(ValueError,r'n110 summary.scans\[0\].minimum_equal_time_relative_pivot'):
                    compare_case_summary(old,new,110)

    def test_native_counts_types_and_inventory_remain_exact(self):
        for old,new in [({'rows':8448},{'rows':8447}),({'rows':1},{'rows':True}),
                        ({'rows':8448},{'rows':8448.}),({'scans':[10]},{'scans':[]}),
                        ({'configuration':110},{'configuration':120}),({'piter':10},{'piter':40}),
                        ({'piter':10},{'piter':10,'extra':1})]:
            with self.subTest(old=old,new=new):
                with self.assertRaises(ValueError):compare_case_summary(old,new,110)

    def test_native_other_floats_remain_exact(self):
        for field,value in [('diagonal_reference_maximum_absolute_difference',0.),
                            ('maximum_relative_residual',8.586715922464113e-12),
                            ('unit_max',7.232611195272963e-16),
                            ('maximum_equal_time_relative_hermitian_defect',0.)]:
            old={field:value};new={field:float(np.nextafter(value,np.inf))}
            with self.subTest(field=field):
                with self.assertRaises(ValueError):compare_case_summary(old,new,110)

    def test_native_nonfinite_diagnostics_are_rejected(self):
        for value in (float('nan'),float('inf'),float('-inf')):
            old=dict(scans=[dict(minimum_equal_time_relative_pivot=.006686613774462862)])
            new=dict(scans=[dict(minimum_equal_time_relative_pivot=value)])
            with self.assertRaises(ValueError):compare_case_summary(old,new,110)

    def test_pivot_tolerance_cannot_apply_to_another_field_path(self):
        old={'minimum_equal_time_relative_pivot':.006686613774462862}
        new={'minimum_equal_time_relative_pivot':float(np.nextafter(old['minimum_equal_time_relative_pivot'],np.inf))}
        with self.assertRaises(ValueError):compare_case_summary(old,new,110)

    def test_archive_paths(self):
        for name in ('../a','/tmp/a','a/../b','a\\b','a//b'):
            with self.assertRaises(ValueError):safe_name(name)
        self.assertEqual(str(safe_name('a/b')),'a/b')

    def test_import_idempotence_and_tamper_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'inputs').mkdir();data=b'abcdef';source=root/'raw.zip';source.write_bytes(data)
            chunks=[dict(path='inputs/review_chunks/part_'+str(i)+'.bin',bytes=3,sha256=hashlib.sha256(part).hexdigest()) for i,part in enumerate((b'abc',b'def'))]
            (root/'inputs/SOURCES.json').write_text(json.dumps(dict(archive_bytes=6,archive_sha256=hashlib.sha256(data).hexdigest(),chunks=chunks)))
            self.assertEqual(import_archive(source,root),2);self.assertEqual(import_archive(source,root),2)
            first=root/chunks[0]['path'];first.write_bytes(b'bad')
            with self.assertRaises(ValueError):import_archive(source,root)
            self.assertEqual(first.read_bytes(),b'bad')
            source.write_bytes(b'truncated')
            with self.assertRaises(ValueError):import_archive(source,root)

    def test_primary_geometry_and_exact_published_diagonals(self):
        a=np.load(str(HERE/'inputs/primary_matrices.npz'),allow_pickle=False)
        self.assertEqual(a['C'].shape,(71,4,33,4,4));self.assertTrue(np.isfinite(a['C']).all())
        b,e,_,_=baseline(a['C'],a['diagonal_text'])
        self.assertEqual((len(b),len(e)),(528,512))

    def test_python36_syntax(self):
        for p in HERE.glob('*.py'):
            try:ast.parse(p.read_text(),filename=p.name,feature_version=(3,6))
            except TypeError:ast.parse(p.read_text(),filename=p.name)


if __name__=='__main__':unittest.main()
