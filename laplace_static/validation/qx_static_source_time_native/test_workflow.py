#!/usr/bin/env python3
"""Native temporal paths, strict analytic replay and complete pilot lifecycle."""
import argparse,ast,contextlib,io,json,os,shutil,subprocess,sys,tempfile,unittest,zipfile
from pathlib import Path
from unittest import mock
import campaign as c
import check_outputs as check
import generate_sources as helper
from parents import HERE,ROOT,MATRIX,DESIGN,load,generator,design
from receipts import digest,read,write
from pilot_plan import PILOT,CUTOFFS,RADII,TAUS,execution_scans,expected_keys,direct_keys,block_keys
from fixed_projections import quadratic
synthetic=load(HERE/'tests/synthetic_outputs.py','qx_source_time_synthetic_fixture')


class OutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.folder=Path(cls.tmp.name)
        synthetic.write_references(cls.folder/'refs.csv');cls.reference=check.load_references(cls.folder/'refs.csv')
        cls.spec=execution_scans(650)[2];cls.text=synthetic.scan_fixture(650,cls.spec)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def replay(self,text,spec=None,reference=None):
        path=self.folder/'scan.out';path.write_text(text)
        return check.scan(path,650,self.spec if spec is None else spec,self.reference if reference is None else reference)

    def test_all_source_times_and_projection_control_replay(self):
        for spec in execution_scans(650):
            values,blocks,summary=self.replay(synthetic.scan_fixture(650,spec),spec)
            self.assertEqual((len(values),len(blocks)),(8448,18432))
            self.assertEqual(summary['fresh_reference_matches'],2112 if spec['t0']==90 else 0)
            self.assertEqual(summary['endpoint_loads'],33)
            self.assertTrue(any(float(r['L_RE'])<0 for r in values));self.assertTrue(any(float(r['L_IM'])!=0 for r in values))

    def test_wrong_source_sink_and_hdf5_load_times_are_rejected(self):
        for prefix,change in [('MATRIX_RT',{'t0':'90'}),('MATRIX_RT',{'t1':'90'}),('MATRIX_DIRECT',{'t0':'90'}),
                ('MATRIX_BLOCK',{'t1':'90'}),('MODE_LOAD',{'loaded_t':'90'}),('SU3_WINDOW',{'t0':'90'}),('SU3_TIME',{'t':'90'})]:
            with self.subTest(prefix=prefix,change=change):
                with self.assertRaises(ValueError):self.replay(synthetic.replace_record(self.text,prefix,change))

    def test_missing_duplicate_and_wrong_matrix_orientation_are_rejected(self):
        line=next(r for r in self.text.splitlines(True) if r.startswith('MATRIX_RT '))
        for changed in (self.text.replace(line,'',1),self.text+line,
                self.text.replace('source_cutoffs=40,60,80,100','source_cutoffs=40,80,60,100',1)):
            with self.assertRaises(ValueError):self.replay(changed)

    def test_missing_endpoint_load_or_duplicate_field_is_rejected(self):
        line=next(r for r in self.text.splitlines(True) if r.startswith('MODE_LOAD '))
        for changed in (self.text.replace(line,'',1),self.text+line,line.replace('status=PASS','status=PASS status=PASS')+self.text):
            with self.assertRaises(ValueError):self.replay(changed)

    def test_counts_nan_axis_reductions_and_half_box_cancellation(self):
        for change,predicate in [({'Nsrc':'1'},lambda r:True),({'L_RE':'nan'},lambda r:True),
                ({'L_RE':'1'},lambda r:r['axis']=='-1'),({'L_IM':'1e-5'},lambda r:r['R']=='32')]:
            with self.assertRaises(ValueError):self.replay(synthetic.replace_record(self.text,'MATRIX_RT',change,predicate))

    def test_fresh_reference_compares_all_16_ordered_pairs(self):
        spec=execution_scans(650)[0];text=synthetic.scan_fixture(650,spec)
        reference=dict(self.reference);key=(650,40,100,32,24);reference[key]=dict(reference[key],L_RE='1')
        with self.assertRaises(ValueError):self.replay(text,spec,reference)

    def test_direct_checks_cannot_inflate_or_hide_errors(self):
        for changed in [synthetic.replace_record(self.text,'MATRIX_DIRECT',{'tolerance':'1'}),
                synthetic.replace_record(self.text,'MATRIX_DIRECT',{'direct_re':'1','abs_error':'0'}),
                synthetic.replace_record(self.text,'MATRIX_LEGACY',{'rank':'1'})]:
            with self.assertRaises(ValueError):self.replay(changed)

    def test_rank_blocks_include_r24_and_correct_rank_coordinates(self):
        line=next(r for r in self.text.splitlines(True) if r.startswith('MATRIX_BLOCK ') and 'R=24 ' in r)
        for changed in (self.text.replace(line,'',1),synthetic.replace_record(self.text,'MATRIX_BLOCK',{'px':'7'}),
                synthetic.replace_record(self.text,'MATRIX_BLOCK',{'rank':'1'}),synthetic.replace_record(self.text,'MATRIX_BLOCK',{'L_RE':'1'})):
            with self.assertRaises(ValueError):self.replay(changed)

    def test_projection_limits_and_periodic_roundtrip_remain_strict(self):
        for prefix,change in [('SU3_TIME',{'unit_max':'2e-12'}),('SU3_WINDOW',{'count':'1'}),
                ('MATRIX_BEGIN',{'preserved_change':'1'}),('PERIODIC_ROUNDTRIP',{'status':'FAIL'})]:
            with self.assertRaises(ValueError):self.replay(synthetic.replace_record(self.text,prefix,change))

    def test_p10_p40_requires_identical_physical_endpoints(self):
        low,_,_=self.replay(self.text)
        high,_,_=self.replay(synthetic.scan_fixture(650,execution_scans(650)[3]),execution_scans(650)[3])
        self.assertEqual(check.projection_pair(low,high),0)
        changed=[dict(r) for r in high];changed[0]['t0']='90'
        with self.assertRaises(ValueError):check.projection_pair(low,changed)
        changed=[dict(r) for r in high];changed[0]['L_IM']='1'
        with self.assertRaises(ValueError):check.projection_pair(low,changed)

    def test_audits_cover_exact_66_and_97_times_and_all100_modes(self):
        for cfg,size in ((170,66),(640,66),(650,97)):
            path=self.folder/'audit.out';path.write_text(synthetic.audit_fixture(cfg))
            summary=check.audit(path,cfg);self.assertEqual((summary['slices'],summary['modes']),(size,size*100))

    def test_audit_rejects_missing_modes_unsorted_values_and_inconsistent_residuals(self):
        text=synthetic.audit_fixture(650);path=self.folder/'bad_audit.out'
        first=next(r for r in text.splitlines(True) if r.startswith('MODE '))
        variants=[text.replace(first,'',1),text+first,synthetic.replace_record(text,'MODE',{'residual_relative':'1.1e-9'}),
            synthetic.replace_record(text,'MODE',{'residual_l2':'1e-10'}),
            synthetic.replace_record(text,'MODE',{'lambda':'.5','rayleigh_re':'.5','residual_l2':str(.5*5e-12)}),
            synthetic.replace_record(text,'SLICE',{'gram_offdiag':'1.1e-9'}),synthetic.replace_record(text,'SLICE',{'max_relative_residual':'1e-13'})]
        for changed in variants:
            path.write_text(changed)
            with self.assertRaises(ValueError):check.audit(path,650)

    def test_sink_row_source_column_quadratic_form_is_not_transposed(self):
        matrix={(40,40):1,(40,100):3+.2j,(100,40):-2+.7j,(100,100):.4}
        value,bound=quadratic(matrix,{40:1+1j,100:2-.5j})
        self.assertAlmostEqual(value.real,3.95);self.assertAlmostEqual(value.imag,-11.15)


class NativeSourceTests(unittest.TestCase):
    def test_python36_and_shell_syntax(self):
        for path in HERE.rglob('*.py'):
            if any(s in ('build','campaigns','__pycache__') for s in path.relative_to(HERE).parts):continue
            try:ast.parse(path.read_text(),filename=str(path),feature_version=(3,6))
            except TypeError:ast.parse(path.read_text(),filename=str(path))
        for name in ('build.sh','submit_micro.sh'):subprocess.run(['bash','-n',str(HERE/name)],check=True)

    def test_generated_io_kernels_and_audit_algebra_are_unchanged(self):
        with tempfile.TemporaryDirectory() as new,tempfile.TemporaryDirectory() as old:
            helper.generate(new);generator.generate(old)
            for name in ('kernels.h','io_helpers.h','mpi_helpers.h','metrics.h'):
                self.assertEqual((Path(new)/name).read_bytes(),(Path(old)/name).read_bytes())
            base=(Path(old)/'audit_selected.c').read_text()
            oldlist='static const int selected_times[]={'+','.join(map(str,range(90,123)))+'};'
            for cfg,size in ((170,66),(650,97)):
                newlist='static const int selected_times[]={'+','.join(map(str,design.audit_times(cfg)))+'};'
                self.assertEqual((Path(new)/('audit%d.c'%size)).read_text().replace(newlist,oldlist),base)
            self.assertNotIn('#define MATRIX_T0',(Path(new)/'matrix_plan.h').read_text())

    @unittest.skipUnless(shutil.which('gcc'),'GCC required for native synthetic validation')
    def test_native_paths_mixed_contractions_and_periodic_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);helper.generate(out);parent=generator.load_parent().parent
            functions='\n'.join(parent.function((out/'kernels.h').read_text(),name) for name in ('mat3_zero','mat3_mul'))
            functions+='\n'+parent.function((out/'io_helpers.h').read_text(),'temporal_product')
            (out/'time_functions.h').write_text(functions)
            exe=out/'time_paths'
            subprocess.run(['gcc','-std=c11','-O2','-I'+tmp,str(HERE/'tests/temporal_paths.c'),'-lm','-o',str(exe)],check=True);subprocess.run([str(exe)],check=True)
            tests=HERE.parent/'qx_static_extension/tests'
            for modes in (4,100):
                exe=out/('mixed%d'%modes)
                subprocess.run(['gcc','-std=c11','-O2','-DNV='+str(modes),'-I'+str(tests),'-I'+tmp,'-I'+str(MATRIX),str(MATRIX/'tests/mixed_contractions.c'),'-lm','-o',str(exe)],check=True)
                subprocess.run([str(exe)],check=True)
            exe=out/'periodic'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(tests),'-I'+tmp,str(tests/'kernels_and_shifts.c'),'-lm','-o',str(exe)],check=True);subprocess.run([str(exe)],check=True)
            include=[tests,out,ROOT/'solver/include',ROOT/'laplace_static/validation/qx_n640_diagnostic/inputs/qcd_include',
                ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230',MATRIX]
            subprocess.run(['gcc','-std=c11','-fsyntax-only','-Werror=implicit-function-declaration','-Werror=format']+['-I'+str(p) for p in include]+[str(HERE/'scan.c')],check=True)


class CampaignTests(unittest.TestCase):
    def isolated(self):
        stack=contextlib.ExitStack();self.addCleanup(stack.close)
        root=Path(stack.enter_context(tempfile.TemporaryDirectory()));here=root/'laplace_static/validation/qx_static_source_time_native';here.mkdir(parents=True)
        (here/'README.md').write_text('Synthetic campaign lifecycle only.\n');write(here/'FILE_MANIFEST.json',dict(files={'README.md':{}},dependencies={}))
        build=here/'build';build.mkdir();write(build/'MANIFEST.json',{})
        design_dir=root/'laplace_static/validation/qx_static_source_time';(design_dir/'inputs').mkdir(parents=True)
        write(design_dir/'PLAN.json',design.build_plan());shutil.copyfile(DESIGN/'inputs/FIXED_SOURCES.json',design_dir/'inputs/FIXED_SOURCES.json')
        reference=design_dir/'inputs/reference_t90_axis_average.csv';synthetic.write_references(reference)
        write(design_dir/'FILE_MANIFEST.json',dict(files={p.relative_to(design_dir).as_posix():{} for p in design_dir.rglob('*') if p.is_file()},dependencies={}))
        physical=root/'physical';physical.mkdir();inputs={}
        for cfg in PILOT:
            folder=physical/str(cfg);folder.mkdir();ev=folder/'vectors';ev.mkdir()
            gauge=folder/('qX_p2n%d'%cfg);gauge.write_bytes(b'synthetic inventory only')
            values=folder/'values.h5';values.write_bytes(b'synthetic table only')
            marker=folder/'COMPLETED.json';write(marker,dict(synthetic=True))
            vectors=[]
            for t in design.audit_times(cfg):
                p=ev/('qX_p2n%d_evec_t%d.h5'%(cfg,t));p.write_bytes(b'synthetic vector inventory only');vectors.append(dict(c.inventory.fingerprint(p),time=t))
            admission=dict(kind='recovered_n640_eigenvalue_checksum',expected_sha256=digest(values)) if cfg==640 else dict(kind='saved_mode_completion_marker',marker=str(marker),sha256=digest(marker))
            inputs[cfg]=dict(configuration=cfg,available=True,required_times=design.audit_times(cfg),errors=[],vectors=vectors,
                gauge=c.inventory.fingerprint(gauge),eigenvalues=dict(c.inventory.fingerprint(values),sha256=digest(values)),admission=admission)
        stack.enter_context(mock.patch.multiple(c,ROOT=root,HERE=here,BUILD=build,DESIGN=design_dir,REFERENCE=reference,SOURCES=design_dir/'inputs/FIXED_SOURCES.json'))
        stack.enter_context(mock.patch.object(c,'verify_sources'));stack.enter_context(mock.patch.object(c,'build_manifest',return_value='synthetic-build-hash'))
        stack.enter_context(mock.patch.object(c,'physical_cases',return_value=inputs))
        stack.enter_context(mock.patch.dict(os.environ,dict(USER='fixture',SLURM_NTASKS='64',SLURM_JOB_NUM_NODES='8',SLURM_JOB_ID='123456',SLURM_ARRAY_TASK_ID='0')))
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        args=argparse.Namespace(name='pilot',concurrent=3,index=0);c.prepare(args)
        return root,here,args,inputs

    def fake_mpi(self,command,**kwargs):
        cfg=int(kwargs['env']['QX_CFG']);piter=int(kwargs['env']['HYP_PITER']);t0=int(kwargs['env']['QX_T0'])
        if Path(command[3]).name.startswith('audit'):text=synthetic.audit_fixture(cfg)
        else:
            spec=next(s for s in execution_scans(cfg) if s['t0']==t0 and s['piter']==piter);text=synthetic.scan_fixture(cfg,spec)
        kwargs['stdout'].write(text);return subprocess.CompletedProcess(command,0)

    def test_full_lifecycle_failed_attempt_retry_reuse_and_archive(self):
        root,here,args,inputs=self.isolated();target=here/'campaigns/pilot'
        with self.assertRaises(SystemExit):c.collect(args)
        self.assertFalse(read(target/'SUMMARY.json')['success'])
        def bad(command,**kwargs):
            kwargs['stdout'].write(synthetic.audit_fixture(170).replace('slices=66','slices=65'))
            return subprocess.CompletedProcess(command,0)
        with mock.patch.object(c.subprocess,'run',side_effect=bad):
            with self.assertRaises(ValueError):c.run(args)
        self.assertTrue((target/'case_000/attempt_0001/FAILED.json').exists());self.assertFalse((target/'case_000/RUNNING.json').exists())
        with mock.patch.object(c.subprocess,'run',side_effect=self.fake_mpi):
            for index in range(3):os.environ['SLURM_ARRAY_TASK_ID']=str(index);c.run(args)
        with mock.patch.object(c.subprocess,'run',side_effect=AssertionError('Validated reuse must not launch MPI')):c.run(args)
        c.collect(args);summary=read(target/'SUMMARY.json')
        self.assertTrue(summary['success']);self.assertEqual((summary['matrix_rows'],summary['rank_block_rows']),(67584,147456))
        self.assertEqual((summary['fixed_source_rows'],summary['paired_fixed_source_rows']),(13068,9504))
        self.assertEqual(summary['statistical_configurations'],3);self.assertEqual(summary['source_training_configurations'],68)
        launches=read(target/'case_002/attempt_0001/LAUNCHES.json')
        self.assertEqual([r['stem'] for r in launches],['audit97','reference90_p10','shift58_p10','shift122_p10','shift122_p40'])
        with zipfile.ZipFile(str(root/'qx_static_source_time_pilot_review.zip')) as z:
            self.assertIsNone(z.testzip());self.assertTrue(any(n.endswith('FAILED.json') for n in z.namelist()));self.assertTrue(any(n.endswith('shift122_p40.err') for n in z.namelist()))
        exported=target/'case_000/attempt_0002/matrix_correlators.csv';exported.write_text('changed')
        with self.assertRaises(ValueError):c.status(args)

    def test_changed_physical_input_and_prepared_scan_order_are_rejected(self):
        root,here,args,inputs=self.isolated()
        Path(inputs[170]['vectors'][0]['path']).write_bytes(b'changed inventory')
        with self.assertRaises(ValueError):c.run(args)
        target=here/'campaigns/pilot';data=read(target/'PREPARED.json');data['cases'][2]['scans'].reverse()
        write(target/'PREPARED.json',data);(target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
        with self.assertRaises(ValueError):c.prepared(args.name)

    def test_allocation_and_campaign_name_guards(self):
        root,here,args,inputs=self.isolated()
        with mock.patch.dict(os.environ,dict(SLURM_NTASKS='1')):
            with self.assertRaises(ValueError):c.run(args)
        for name in ('../other','Pilot','a/b'):
            with self.assertRaises(ValueError):c.folder(name)

    def test_failed_fresh_regression_stops_before_shifted_launches(self):
        root,here,args,inputs=self.isolated();os.environ['SLURM_ARRAY_TASK_ID']='2'
        launched=[]
        def wrong_reference(command,**kwargs):
            launched.append(Path(command[3]).name)
            if launched[-1].startswith('audit'):return self.fake_mpi(command,**kwargs)
            spec=execution_scans(650)[0]
            kwargs['stdout'].write(synthetic.replace_record(synthetic.scan_fixture(650,spec),'MATRIX_RT',{'L_RE':'1'},lambda r:r['axis']=='-1'))
            return subprocess.CompletedProcess(command,0)
        with mock.patch.object(c.subprocess,'run',side_effect=wrong_reference):
            with self.assertRaises(ValueError):c.run(args)
        self.assertEqual(launched,['audit97','scan'])
        case=here/'campaigns/pilot/case_002'
        self.assertFalse((case/'LATEST.json').exists());self.assertTrue((case/'attempt_0001/FAILED.json').exists())

    def test_nonzero_mpi_exit_retains_logs_without_a_validated_marker(self):
        root,here,args,inputs=self.isolated()
        def failed(command,**kwargs):
            kwargs['stderr'].write('Synthetic MPI failure\n');return subprocess.CompletedProcess(command,17)
        with mock.patch.object(c.subprocess,'run',side_effect=failed):
            with self.assertRaises(ValueError):c.run(args)
        case=here/'campaigns/pilot/case_000';attempt=case/'attempt_0001'
        self.assertFalse((case/'LATEST.json').exists())
        self.assertEqual(read(attempt/'LAUNCHES.json')[0]['exit_code'],17)
        self.assertEqual((attempt/'audit66.err').read_text(),'Synthetic MPI failure\n')

    def test_collection_refuses_an_active_lock_before_writing(self):
        root,here,args,inputs=self.isolated();target=here/'campaigns/pilot'
        case=target/'case_000';case.mkdir();write(case/'RUNNING.json',dict(job_id='123456',array_task_id=0))
        with self.assertRaises(ValueError):c.collect(args)
        self.assertFalse((target/'SUMMARY.json').exists());self.assertFalse((root/'qx_static_source_time_pilot_review.zip').exists())

    def test_submission_receipts_block_duplicates_and_ambiguous_retries(self):
        root,here,args,inputs=self.isolated();target=here/'campaigns/pilot'
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout='777\n',stderr='')):c.submit(args)
        self.assertTrue((target/'SUBMISSION_001.json').exists())
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout='777 777\n',stderr='')):
            with self.assertRaises(ValueError):c.submit(args)
        write(target/'SUBMISSION_001.json',dict(state='submitting'))
        with mock.patch.object(c.subprocess,'run',side_effect=AssertionError('Ambiguous receipt must not resubmit')):
            with self.assertRaises(ValueError):c.submit(args)

    def test_lock_recovery_requires_a_confirmed_stopped_job(self):
        root,here,args,inputs=self.isolated();case=here/'campaigns/pilot/case_000';case.mkdir()
        write(case/'RUNNING.json',dict(job_id='123456',array_task_id=0))
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout='123456\n',stderr='')):
            with self.assertRaises(ValueError):c.recover_lock(args)
        self.assertTrue((case/'RUNNING.json').exists())
        with mock.patch.object(c.subprocess,'run',return_value=subprocess.CompletedProcess([],0,stdout='',stderr='')):c.recover_lock(args)
        self.assertFalse((case/'RUNNING.json').exists());self.assertEqual(len(list(case.glob('RECOVERED_RUNNING_*.json'))),1)


if __name__=='__main__':unittest.main()
