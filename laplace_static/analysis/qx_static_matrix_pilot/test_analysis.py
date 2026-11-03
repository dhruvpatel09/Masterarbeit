#!/usr/bin/env python3
"""Exact spectral controls, complex-mode semantics and real-pilot replay tests."""
import csv
import json
import math
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock
import numpy as np
from scipy import linalg
import analysis as a
from replay import compare,extract_checked,replay
from reproduce import read_mode_labels,verify_exports,verify_table


class SpectralTests(unittest.TestCase):
    def series(self):
        z=np.array([[1,.2+.1j,.1],[.15j,.9,.25],[.2,.1j,1.2]],complex)
        energies=np.array([.2,.4,.7])+1j*np.array([.01,-.02,.03])
        return np.array([z@np.diag(np.exp(-energies*t))@z.conj().T for t in range(33)]),energies

    def test_exact_complex_spectrum_and_reference_times(self):
        matrices,energies=self.series()
        for tref in (0,1,2,4):
            result=a.solve_series(matrices,tref)
            for t in range(tref+1,33):
                np.testing.assert_allclose(result['values'][t],np.exp(-energies*(t-tref)),rtol=2e-8,atol=1e-12)
                self.assertLess(max(result['residuals'][t]),1e-13)
            effective,_=a.energy(result['values'][tref+1,0],result['values'][tref+2,0],.05)
            self.assertAlmostEqual(effective.real,.2,places=12);self.assertAlmostEqual(effective.imag,.01,places=12)

    def test_complex_nonsymmetric_entries_are_preserved(self):
        matrices,_=self.series();saved=matrices.copy();result=a.solve_series(matrices,0)
        np.testing.assert_array_equal(matrices,saved)
        self.assertGreater(np.max(np.abs(result['values'][4].imag)),.001)
        self.assertGreater(linalg.norm(matrices[4]-matrices[4].conj().T),.001)

    def test_basis_rescaling_preserves_generalized_eigenvalues(self):
        matrices,_=self.series();scale=np.diag([.1,7.,2.]);changed=np.array([scale@c@scale for c in matrices])
        first=a.solve_series(matrices,0);second=a.solve_series(changed,0)
        np.testing.assert_allclose(first['values'][1:],second['values'][1:],rtol=1e-8,atol=1e-12)

    def test_equal_time_whitening_and_rank_failure(self):
        matrices,_=self.series();transform,_,eigenvalues=a.metric_transform(matrices[0])
        np.testing.assert_allclose(transform.conj().T@matrices[0]@transform,np.eye(3),atol=1e-12)
        self.assertGreater(eigenvalues[0],0)
        for matrix in (np.ones((2,2),complex),np.diag([1.,-1.]),np.array([[1,1j],[0,1]],complex)):
            with self.assertRaises(ValueError):a.metric_transform(matrix)

    def test_negative_and_large_phase_values_do_not_become_positive_energies(self):
        for first,second in ((-1+0j,-.5+0j),(1+.3j,.5+.1j),(0j,1+0j)):
            _,valid=a.energy(first,second,.05);self.assertFalse(valid)
        effective,valid=a.energy(1+0j,.5+0j,.05);self.assertTrue(valid);self.assertAlmostEqual(effective.real,np.log(2))

    def test_vector_matching_is_invariant_under_phase_and_permutation(self):
        vectors=np.eye(3,dtype=complex);changed=vectors[:,[2,0,1]]*np.exp(1j*np.array([.2,.9,-.5]))
        order,score=a.match(vectors,changed,np.eye(3));self.assertEqual(list(order),[1,2,0]);np.testing.assert_allclose(score,1)

    def test_numeric_replay_accepts_roundoff_and_rejects_material_changes(self):
        compare({'value':1.+1e-13},{'value':1.})
        with self.assertRaises(ValueError):compare({'value':1.01},{'value':1.})
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'table.csv';path.write_text('value,label\n1.000000000001,a\n')
            verify_table(path,[dict(value=1.,label='a')])
            with self.assertRaises(ValueError):verify_table(path,[dict(value=1.1,label='a')])

    def test_principal_phase_comparison_wraps_pi_and_rejects_changed_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'angles.csv'
            path.write_text('phase,effective_energy_im,label\n%.17g,%.17g,a\n'%(-math.pi+1e-13,math.pi-1e-13))
            verify_table(path,[dict(phase=math.pi-1e-13,effective_energy_im=-math.pi+1e-13,label='a')])
            with self.assertRaises(ValueError):verify_table(path,[dict(phase=math.pi,effective_energy_im=.7182888468862978,label='a')])
            path.write_text('effective_energy_im\n-0.7182888468862978\n')
            with self.assertRaises(ValueError):verify_table(path,[dict(effective_energy_im=.7182888468862978)])
            path.write_text('phase\n%.17g\n'%(3*math.pi))
            with self.assertRaises(ValueError):verify_table(path,[dict(phase=math.pi)])
            for value in ('nan','inf','-inf'):
                path.write_text('phase\n'+value+'\n')
                with self.assertRaises(ValueError):verify_table(path,[dict(phase=0.)])

    def test_spectrum_replay_requires_every_independently_computed_eigenvalue(self):
        labels=np.array([.4+.2j,.4-.2j,.8+0j]);fresh=labels[[2,1,0]]
        order=a.match_spectrum(labels,fresh);np.testing.assert_allclose(fresh[order],labels)
        changed=labels.copy();changed[0]+=.01
        with self.assertRaises(ValueError):a.match_spectrum(changed,fresh)
        with self.assertRaises(ValueError):a.match_spectrum(labels[:2],fresh)
        with self.assertRaises(ValueError):a.match_spectrum(np.array([np.nan,1,2]),fresh)

    def test_series_replay_labels_survive_eigenvector_phases_and_permutation(self):
        matrices,_=self.series();reference=a.solve_series(matrices,0);original=a.linalg.eig
        def permuted(*args,**kwargs):
            values,vectors=original(*args,**kwargs);order=np.arange(len(values))[::-1]
            return values[order],vectors[:,order]*np.exp(1j*np.linspace(.2,.9,len(values)))
        with mock.patch.object(a.linalg,'eig',side_effect=permuted):
            reproduced=a.solve_series(matrices,0,labels=reference['values'])
        np.testing.assert_allclose(reproduced['values'][1:],reference['values'][1:],rtol=1e-10,atol=1e-12)
        np.testing.assert_allclose(reproduced['vectors'][1:],reference['vectors'][1:],rtol=1e-8,atol=1e-10)
        labels=reference['values'].copy();labels[4,0]+=.1
        with self.assertRaises(ValueError):a.solve_series(matrices,0,labels=labels)

    def test_review_extraction_rejects_unexpected_or_escaping_entries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for i,name in enumerate(('../escape','different/source.py','laplace_static/validation/qx_static_matrix/../../escape')):
                path=root/('bad%d.zip'%i)
                with zipfile.ZipFile(path,'w') as z:z.writestr(name,'x')
                with self.assertRaises(ValueError):extract_checked(path,root/'output')


class PilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows,cls.blocks,cls.validation=replay();cls.tables,cls.summary=a.compute(cls.rows,cls.validation)

    def test_raw_mpi_replay_and_exact_diagonal_controls(self):
        v=self.validation;self.assertTrue(v['success']);self.assertEqual((len(self.rows),len(self.blocks)),(33792,49152))
        self.assertEqual(v['p10_diagonal_reference_rows'],6336);self.assertEqual(v['p10_diagonal_maximum_absolute_difference'],0.)
        self.assertLess(v['p10_p40_full_matrix_maximum_component_difference'],1e-20)
        self.assertEqual((v['direct_tests'],v['legacy_tests'],v['audited_modes']),(6144,1536,9900))

    def test_explicit_source_column_and_sink_row_orientation(self):
        plan=json.loads((a.HERE/'ANALYSIS_PLAN.json').read_text());data,_=a.assemble(self.rows,plan)
        row=next(r for r in self.rows if r['cfg']=='170' and r['piter']=='10' and r['tau']=='4' and r['R']=='12' and r['axis']=='-1' and r['Nv_src']=='40' and r['Nv_sink']=='100')
        self.assertEqual(data[0,1,4,3,0],complex(float(row['L_RE']),float(row['L_IM'])))
        self.assertNotEqual(data[0,1,4,3,0],data[0,1,4,0,3])

    def test_full_rank_metrics_and_smaller_basis_conditioning(self):
        for radius,conditions in self.summary['equal_time_diagonally_normalized_conditions'].items():
            self.assertGreater(conditions['full4'],1);self.assertLess(conditions['wide2'],conditions['full4'])
        self.assertEqual(self.summary['empirical_covariance_rank'],2)
        self.assertEqual(self.summary['statistical_sample_size'],3)
        self.assertFalse(self.summary['production_submitted'])

    def test_all_numerical_residuals_and_heldout_training_separation(self):
        self.assertLess(self.summary['maximum_generalized_eigenpair_relative_residual'],1e-10)
        selected=[r for r in self.tables['fixed_projections'] if r['evaluation']=='heldout']
        self.assertTrue(selected)
        for row in selected:self.assertFalse(set(row['training_configurations'].split(';'))&set(row['evaluation_configurations'].split(';')))

    def test_portable_real_pilot_replay_under_eigensolver_order_and_roundoff(self):
        plan=json.loads((a.HERE/'ANALYSIS_PLAN.json').read_text());labels=read_mode_labels(a.HERE/'results/complex_modes.csv',plan)
        original=a.linalg.eig;calls=[0]
        def varied(*args,**kwargs):
            values,vectors=original(*args,**kwargs);order=np.roll(np.arange(len(values)),calls[0]%len(values));calls[0]+=1
            values=values[order].copy();vectors=vectors[:,order]*np.exp(1j*np.linspace(.17,.83,len(values)))
            tiny=abs(values.imag)<=32*np.finfo(float).eps*max(abs(values))
            values[tiny]=values[tiny].conjugate()
            return values,vectors
        with mock.patch.object(a.linalg,'eig',side_effect=varied):tables,summary=a.compute(self.rows,self.validation,mode_labels=labels)
        verify_exports(tables,summary,self.validation)

    def test_saved_mode_label_inventory_rejects_duplicate_missing_and_wrong_geometry(self):
        plan=json.loads((a.HERE/'ANALYSIS_PLAN.json').read_text());source=a.HERE/'results/complex_modes.csv'
        with source.open(newline='') as f:rows=list(csv.DictReader(f))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'labels.csv'
            for values in (rows[:-1],rows+[rows[0]],rows+[dict(rows[0],tau='0')]):
                with path.open('w',newline='') as f:
                    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(values)
                with self.assertRaises(ValueError):read_mode_labels(path,plan)


if __name__=='__main__':unittest.main()
