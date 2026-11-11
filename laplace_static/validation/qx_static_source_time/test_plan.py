#!/usr/bin/env python3
"""Check geometry, source isolation and failure handling without MPI or Slurm."""
import ast,csv,hashlib,json,math,tempfile,unittest
from pathlib import Path
from plan import HERE,ROOT,CONFIGS,PILOT,SOURCES,TAUS,RADII,CUTOFFS,BASES,ANCHORS,interval,scans,audit_times,matrix_keys,build_plan
from preflight import inspect_inputs,save_report


class GeometryTests(unittest.TestCase):
    def test_all_intervals_stay_inside_open_time(self):
        for t0 in SOURCES:
            for tau in TAUS:
                r=interval(t0,tau)
                self.assertEqual(r['t1']-r['t0'],tau)
                self.assertGreater(r['left_boundary_distance'],0)
                self.assertGreaterEqual(r['right_boundary_distance'],37)

    def test_unplanned_sources_and_separations_are_rejected(self):
        for t0,tau in ((0,4),(191,0),(90,33),(90,-1)):
            with self.assertRaises(ValueError):interval(t0,tau)

    def test_mode_audit_inventory_includes_every_endpoint(self):
        for cfg in PILOT:
            times=set(audit_times(cfg))
            for scan in scans(cfg):
                for tau in TAUS:self.assertTrue({scan['t0'],scan['t0']+tau}<=times)
        self.assertEqual([len(audit_times(c)) for c in PILOT],[66,66,97])

    def test_control_configuration_audits_the_reference_gap(self):
        self.assertNotIn(110,audit_times(170));self.assertIn(110,audit_times(650))
        self.assertEqual(sum(len(audit_times(c)) for c in PILOT),229)

    def test_projection_and_fresh_reference_controls_are_explicit(self):
        all_scans=[(c,s['t0'],s['piter'],s['purpose']) for c in PILOT for s in scans(c)]
        self.assertEqual(len(all_scans),8)
        self.assertEqual([s for s in all_scans if s[2]==40],[(650,122,40,'projection_control')])
        self.assertEqual([s for s in all_scans if s[1]==90],[(650,90,10,'fresh_reference_regression')])

    def test_full_ordered_matrix_and_axis_inventory(self):
        keys=matrix_keys();self.assertEqual(len(keys),8448)
        self.assertIn((40,100,32,24,-1),keys);self.assertIn((100,40,32,24,-1),keys)
        self.assertEqual({k[-1] for k in keys},{-1,0,1,2})

    def test_unknown_configuration_is_rejected(self):
        with self.assertRaises(ValueError):scans(180)

    def test_plan_is_a_design_not_a_submission(self):
        plan=build_plan();self.assertFalse(plan['native_extension_implemented']);self.assertFalse(plan['submission_implemented'])
        self.assertEqual(plan['planned_total_mpi_launches'],11)
        self.assertEqual(json.loads((HERE/'PLAN.json').read_text()),plan)


class FrozenReferenceTests(unittest.TestCase):
    def test_training_excludes_every_pilot_configuration(self):
        data=json.loads((HERE/'inputs/FIXED_SOURCES.json').read_text())
        self.assertEqual(data['training_configuration_labels'],[c for c in CONFIGS if c not in PILOT])
        self.assertEqual(len(data['training_configuration_labels']),68)
        self.assertFalse(set(data['training_configuration_labels'])&set(PILOT))

    def test_all_fixed_vectors_have_the_declared_basis_and_normalization(self):
        data=json.loads((HERE/'inputs/FIXED_SOURCES.json').read_text());vectors=data['source_candidates']
        self.assertEqual({(v['R'],v['basis'],v['anchor']) for v in vectors},{(r,b,a) for r in RADII for b in BASES for a in ANCHORS})
        self.assertEqual(len(vectors),36)
        for v in vectors:
            self.assertEqual(v['Nv'],list(BASES[v['basis']]))
            self.assertEqual(len(v['coefficients_re']),len(v['Nv']));self.assertEqual(len(v['coefficients_im']),len(v['Nv']))
            self.assertTrue(all(math.isfinite(x) for x in v['coefficients_re']+v['coefficients_im']))
            self.assertAlmostEqual(v['reference_metric_norm'],1,delta=1e-9)
            self.assertLessEqual(v['diagnostics']['rayleigh_residual'],1e-10)

    def test_baseline_has_every_ordered_pair_with_correct_physical_labels(self):
        with (HERE/'inputs/reference_t90_axis_average.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        keys={(int(d['configuration']),int(d['R']),int(d['tau']),int(d['Nv_src']),int(d['Nv_sink'])) for d in rows}
        self.assertEqual(len(rows),6336)
        self.assertEqual(keys,{(c,r,t,s,k) for c in PILOT for r in RADII for t in TAUS for s in CUTOFFS for k in CUTOFFS})
        for d in rows:
            self.assertEqual(int(d['t1']),90+int(d['tau']));self.assertEqual(d['t0'],'90');self.assertEqual(d['axis'],'-1')
            self.assertEqual(d['Nsrc'],str(3*64**3));self.assertEqual(d['piter'],'10')
            self.assertTrue(all(math.isfinite(float(d[k])) for k in ('L_RE','L_IM')))

    def test_review_preserves_the_known_late_time_warning(self):
        rows=json.loads((HERE/'inputs/DIAGNOSTIC_REVIEW.json').read_text())['radii']
        r=next(x for x in rows if x['R']==32)
        self.assertLess(r['late_minus_early'],-3*r['paired_window_se'])
        self.assertGreater(r['late_gaussian_tail_approx'],.5)

    def test_reference_equal_time_matrices_preserve_complex_gram_symmetry(self):
        with (HERE/'inputs/reference_t90_axis_average.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        matrix={(int(d['configuration']),int(d['R']),int(d['Nv_src']),int(d['Nv_sink'])):complex(float(d['L_RE']),float(d['L_IM'])) for d in rows if d['tau']=='0'}
        for (cfg,r,s,k),value in matrix.items():
            other=matrix[cfg,r,k,s].conjugate()
            self.assertLessEqual(abs(value-other),1e-20+5e-12*max(abs(value),abs(other)))


class PhysicalInventoryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.gauges=self.root/'gauges';self.gauges.mkdir()
        self.runs=[dict(configuration=c,run=str(self.root/('unused%d'%c))) for c in CONFIGS]
        for cfg in PILOT:
            run=self.root/('n%d'%cfg);ev=run/'eigenvectors';values=run/'eigenvalues';ev.mkdir(parents=True);values.mkdir()
            table=values/'values.h5';table.write_bytes(b'synthetic-file-inventory-only')
            marker=run/'COMPLETED.json';marker.write_text('{"synthetic":true}\n')
            (self.gauges/('qX_p2n%d'%cfg)).write_bytes(b'synthetic-gauge-inventory-only')
            for t in audit_times(cfg):(ev/('qX_p2n%d_evec_t%d.h5'%(cfg,t))).write_bytes(b'synthetic-mode-inventory-only')
            entry=self.runs[CONFIGS.index(cfg)]
            if cfg==640:entry.update(eigenvectors=str(ev),eigenvalues=str(table),eigenvalue_sha256=hashlib.sha256(table.read_bytes()).hexdigest())
            else:entry.update(run=str(run),completion_marker_sha256=hashlib.sha256(marker.read_bytes()).hexdigest())

    def tearDown(self):self.tmp.cleanup()

    def test_complete_inventory_is_not_claimed_as_native_validation(self):
        r=inspect_inputs(self.runs,self.gauges)
        self.assertTrue(r['success']);self.assertEqual(r['inspected_mode_files'],229)
        self.assertFalse(r['native_mode_audit_completed']);self.assertFalse(r['native_extension_built']);self.assertFalse(r['jobs_submitted'])

    def test_missing_required_shifted_time_blocks_that_case(self):
        (self.root/'n170/eigenvectors/qX_p2n170_evec_t154.h5').unlink()
        r=inspect_inputs(self.runs,self.gauges)
        self.assertFalse(r['success']);self.assertEqual(r['available_configurations'],[640,650])
        self.assertIn('t154.h5',r['cases'][0]['errors'][0])

    def test_reference_only_time_is_required_on_n650(self):
        (self.root/'n650/eigenvectors/qX_p2n650_evec_t110.h5').unlink()
        r=inspect_inputs(self.runs,self.gauges)
        self.assertEqual(r['available_configurations'],[170,640])

    def test_empty_physical_file_is_rejected(self):
        (self.gauges/'qX_p2n170').write_bytes(b'')
        self.assertFalse(inspect_inputs(self.runs,self.gauges)['cases'][0]['available'])

    def test_modified_completion_marker_is_rejected(self):
        (self.root/'n650/COMPLETED.json').write_text('{}')
        r=inspect_inputs(self.runs,self.gauges)
        self.assertFalse(r['cases'][2]['available']);self.assertIn('marker changed',r['cases'][2]['errors'][0])

    def test_recovered_n640_admission_checksum_is_enforced(self):
        (self.root/'n640/eigenvalues/values.h5').write_bytes(b'changed')
        r=inspect_inputs(self.runs,self.gauges)
        self.assertFalse(r['cases'][1]['available']);self.assertIn('admission receipt changed',r['cases'][1]['errors'][0])

    def test_multiple_eigenvalue_tables_are_rejected(self):
        (self.root/'n170/eigenvalues/extra.h5').write_bytes(b'extra')
        r=inspect_inputs(self.runs,self.gauges)
        self.assertFalse(r['cases'][0]['available'])

    def test_duplicate_or_incomplete_configuration_registry_is_rejected(self):
        for bad in (self.runs[:-1],[self.runs[0]]+self.runs[:-1]):
            with self.assertRaises(ValueError):inspect_inputs(bad,self.gauges)

    def test_output_is_kept_inside_the_repository(self):
        with self.assertRaises(ValueError):save_report(self.root/'outside.json',{})

    def test_report_cannot_overwrite_a_tracked_plan_file(self):
        with self.assertRaises(ValueError):save_report(HERE/'PLAN.json',{})


class SyntaxTests(unittest.TestCase):
    def test_python36_syntax(self):
        for p in HERE.glob('*.py'):
            try:ast.parse(p.read_text(),filename=p.name,feature_version=(3,6))
            except TypeError:ast.parse(p.read_text(),filename=p.name)


if __name__=='__main__':unittest.main()
