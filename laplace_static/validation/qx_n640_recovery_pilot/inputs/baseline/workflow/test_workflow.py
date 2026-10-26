#!/usr/bin/env python3
"""Test trace semantics, control gates and observation-only numerical equivalence."""
import ctypes.util
import contextlib
import importlib.util
import io
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import unittest
import zipfile
from unittest import mock
from pathlib import Path
import campaign
import instrument_lanczos

HERE=Path(__file__).resolve().parent
SPEC=json.loads((HERE/'INPUTS.json').read_text())
PARAM=SPEC['parameters']


def fixture(code=0,variant='traced',bad_audit=False):
    cfg=640 if code else 650
    case={'configuration':cfg,'time':1,'variant':variant}
    lines=['QX_CASE configuration=%d time=1 variant=%s Nv=100 ranks=64 spatial_ranks=64 P=1,4,4,4 res=1e-9 redTol=1e-12 deg=10 upper=.15 tolReortho=1e-10 maxRestarts=30 nLanczosFactor=7 redFreqFactor=10 nAPE=20 alphaAPE=.5 mpi_thread_min=1'%(cfg,variant),
           'QX_GAUGE plaquette=.66 load_seconds=1',
           'QX_LINKS stage=raw count=786432 nonfinite_links=0 unitarity_max=1e-15 determinant_max=1e-15',
           'QX_APE steps=20 alpha=.5 seconds=1',
           'QX_LINKS stage=ape count=786432 nonfinite_links=0 unitarity_max=1e-15 determinant_max=1e-15']
    cycles=30 if code==2 else 2 if code==3 else 0 if code==1 else 1
    final_res=2e-8 if code==2 else 1e-12
    if variant=='traced':
        for cycle in range(1,cycles+1):
            lines.append('QX_CYCLE_BEGIN cycle=%d m=700'%cycle)
            failed=(code in (3,4) and cycle==cycles)
            method='DSTEMR' if cycle==1 else 'DSYEVR'
            lines.append('QX_REDUCED cycle=%d method=%s info=%d m=700'%(cycle,method,-71 if failed else 0))
            if failed:continue
            checked=100 if cycle==cycles else 1
            residual=final_res if cycle==cycles else 1e-3
            for index in range(99,99-checked,-1):
                eigenvalue=.2 if code==5 else .01+index*.0001
                lines.append('QX_RITZ cycle=%d index=%d lambda=%.17e residual_norm=1e-12 relative_residual=%.17e'%(cycle,index,eigenvalue,residual))
            lines.append('QX_CYCLE_END cycle=%d checked=%d max_checked_res=%.17e lambda_Nv=%.17e restarts=%d nRed=%d maxRed=700 reorthogonalizations=0 seconds=1'%(cycle,checked,residual,.2 if code==5 else .0199,cycle-1,cycle))
    normal=code in (0,2,5)
    lines.append('QX_RETURN code=%d code_min=%d code_max=%d restarts=%d nRed=%d maxRed=%d maxRes=%.17e usedLapack=1 seconds=1'%(code,code,code,cycles-1 if normal else 0,cycles if normal else 0,700 if normal else 0,final_res if normal else 100.))
    if normal:
        for index in range(100):
            eigenvalue=.2 if code==5 else .01+index*.0001
            lines.append('QX_AUDIT_MODE index=%d lambda=%.17e norm=1 relative_residual=%s finite_positive=%d'%(index,eigenvalue,'nan' if bad_audit else '%.17e'%final_res,0 if bad_audit else 1))
        lines.append('QX_AUDIT modes=100 invalid=%d max_relative_residual=%.17e max_norm_defect=0 max_gram_offdiag=0 valid_selected_slice=%d'%(100 if bad_audit else 0,0. if bad_audit else final_res,int(code==0 and not bad_audit)))
    else:lines.append('QX_AUDIT skipped=1 reason=solver_return_code')
    lines.append('QX_COMPLETE configuration=%d time=1 solver_code=%d'%(cfg,code))
    return '\n'.join(lines)+'\n',case


class ProtocolTests(unittest.TestCase):
    def test_completed_nonconvergence_is_not_configuration_admission(self):
        text,case=fixture(2)
        report=campaign.parse_log(text,case,PARAM)
        self.assertTrue(report['diagnostic_completed'])
        self.assertFalse(report['selected_slice_valid'])
        self.assertFalse(report['production_configuration_admitted'])
        self.assertEqual(report['solver_result'],'residual_nonconvergence')
        self.assertEqual(report['cycles'],30)

    def test_return_codes_are_preserved(self):
        for code in (0,1,3,4,5):
            text,case=fixture(code)
            report=campaign.parse_log(text,case,PARAM)
            self.assertEqual(report['solver_code'],code)
            self.assertEqual(report['selected_slice_valid'],code==0)

    def test_nonfinite_audit_blocks_false_success(self):
        text,case=fixture(0,bad_audit=True)
        report=campaign.parse_log(text,case,PARAM)
        self.assertFalse(report['selected_slice_valid'])

    def test_truncated_or_wrong_geometry_and_parameters_rejected(self):
        text,case=fixture(0)
        changes=[text.replace('QX_COMPLETE','INCOMPLETE'),text.replace('res=1e-9','res=1e-8'),
                 text.replace('ranks=64 spatial','ranks=32 spatial'),
                 '\n'.join(line for line in text.splitlines() if not line.startswith('QX_AUDIT_MODE index=99 '))]
        for changed in changes:
            with self.assertRaises(ValueError):campaign.parse_log(changed,case,PARAM)

    def test_partial_checks_not_reported_as_full_maximum(self):
        text,case=fixture(2)
        report=campaign.parse_log(text,case,PARAM)
        self.assertEqual(report['cycle_history'][0]['checked'],'1')
        self.assertEqual(report['cycle_history'][-1]['checked'],'100')
        with self.assertRaises(ValueError):campaign.parse_log(text.replace('cycle=1 checked=1 ','cycle=1 checked=100 '),case,PARAM)

    def test_original_library_control_and_gate(self):
        a,ca=fixture(0,'library');b,cb=fixture(0,'traced')
        reports=[campaign.parse_log(a,ca,PARAM),campaign.parse_log(b,cb,PARAM)]
        self.assertEqual(campaign.compare_controls(reports,1e-10)['maximum_absolute_eigenvalue_difference'],0)
        reports[1]['audited_modes'][0]['lambda']='.1'
        with self.assertRaises(ValueError):campaign.compare_controls(reports,1e-10)

    def test_preserved_inputs_and_instrumentation_anchor(self):
        for name,digest in SPEC['preserved_inputs'].items():self.assertEqual(campaign.sha(HERE/name),digest)
        source=(HERE/'inputs/qcd_lanczos.c').read_text()
        traced=instrument_lanczos.instrument(source)
        self.assertEqual(source.split('int qcd_lanczosCheby(')[0],traced.split('int qcd_lanczosCheby(')[0])
        self.assertEqual(source.split('int qcd_lanczosReortho(')[1],traced.split('int qcd_lanczosReortho(')[1])
        with self.assertRaises(ValueError):instrument_lanczos.instrument(traced)

    def test_runtime_output_not_in_package_pins(self):
        before=campaign.package_pins()
        with tempfile.NamedTemporaryFile(dir=str(HERE),prefix='diagnostic_fixture_',suffix='.out') as f:
            self.assertEqual(campaign.package_pins(),before)

    def test_campaign_names_are_local(self):
        for name in ('../outside','a/b',''):
            with self.assertRaises(ValueError):campaign.directory(name)


class NativeTests(unittest.TestCase):
    def test_driver_syntax_with_captured_qcd_headers(self):
        compiler=shutil.which('gcc')
        if compiler is None:self.skipTest('GCC not available')
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);(folder/'diagnostic_config.h').write_text(campaign.config_header(SPEC))
            subprocess.run([compiler,'-std=c11','-Werror=implicit-function-declaration','-fsyntax-only',
                '-I'+str(folder),'-I'+str(HERE/'inputs/qcd_include'),'-I'+str(HERE/'tests/shim'),str(HERE/'test_n640.c')],check=True)

    def test_observation_only_source_equivalence(self):
        compiler=shutil.which('gcc');lapack=ctypes.util.find_library('lapack');blas=ctypes.util.find_library('blas')
        if not compiler or not lapack or not blas:self.skipTest('GCC and system LAPACK/BLAS needed for dense fixture')
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary);instrument_lanczos.generate(folder/'traced.c')
            for variant,source in [('original',HERE/'inputs/qcd_lanczos.c'),('traced',folder/'traced.c')]:
                subprocess.run([compiler,'-std=c11','-O2','-I'+str(HERE/'tests/shim'),'-I'+str(HERE/'inputs/qcd_include'),
                    str(source),str(HERE/'tests/dense_backend.c'),str(HERE/'tests/synthetic_main.c'),
                    '-l:'+lapack,'-l:'+blas,'-lm','-o',str(folder/variant)],check=True)
            for name,restarts,upper,factor,failure,expected in [
                    ('success',30,'.15',6,None,0),('nonconvergence',1,'.15',4,None,2),
                    ('cutoff',30,'.025',6,None,5),('dsyevr',30,'.15',4,'DSYEVR',3),('dstemr',30,'.15',4,'DSTEMR',4)]:
                environment=os.environ.copy()
                environment.pop('QX_TEST_REDUCED_FAILURE',None)
                if failure:environment['QX_TEST_REDUCED_FAILURE']=failure
                outputs=[]
                for variant in ('original','traced'):
                    output=folder/(name+'_'+variant+'.bin')
                    result=subprocess.run([str(folder/variant),str(output),str(restarts),upper,str(factor)],env=environment,
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,check=True)
                    outputs.append(output.read_bytes())
                    self.assertEqual(struct.unpack('=i',outputs[-1][:4])[0],expected)
                    if name=='success':
                        self.assertAlmostEqual(float(re.search(r'lambda0=(\S+)',result.stdout).group(1)),.01,places=12)
                        self.assertAlmostEqual(float(re.search(r'lambda2=(\S+)',result.stdout).group(1)),.03,places=12)
                self.assertEqual(outputs[0],outputs[1],name+' changed numerical output under observation-only patch')


class LifecycleTests(unittest.TestCase):
    def test_preparation_control_gate_run_collection_and_source_guard(self):
        # Fake launchers exercise workflow bookkeeping, not distributed solver numerics.
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)/'repo';base=root/'laplace_static/validation/qx_n640_diagnostic'
            shutil.copytree(str(HERE),str(base),ignore=shutil.ignore_patterns('build','campaigns','__pycache__','*.out','*.err'))
            spec=importlib.util.spec_from_file_location('staged_qx_diagnostic',base/'campaign.py')
            staged=importlib.util.module_from_spec(spec);spec.loader.exec_module(staged)
            data=json.loads((base/'INPUTS.json').read_text())
            reference=root/'fixture_reference';reference.write_text('fixture\n')
            data['repository_references']={'fixture_reference':staged.sha(reference)}
            gauges=root/'gauges';gauges.mkdir()
            for configuration in (640,650):(gauges/('qX_p2n%d'%configuration)).write_bytes(struct.pack('<4i',192,64,64,64))
            data['gauge_directory']=str(gauges);staged.write_json(base/'INPUTS.json',data)
            build=base/'build';build.mkdir()
            for name in ('test_n640','test_n640_reference','test_n640_reference.o','test_n640.o','qcd_lanczos_trace.o',
                         'qcd_lanczos_trace.c','diagnostic_config.h','compiler.txt','modules.txt','linked_libraries.txt','link.map'):
                (build/name).write_text('fixture build artifact\n')
            executables=root/'bin';executables.mkdir()
            launcher=executables/'mpiexec'
            launcher.write_text('#!/usr/bin/env python3\nimport os,sys\nsys.path.insert(0,'+repr(str(HERE))+')\n'
                'from test_workflow import fixture\nconfiguration=int(sys.argv[-2])\n'
                'variant="library" if "reference" in sys.argv[-4] else "traced"\n'
                'bad=configuration==650 and os.environ.get("QX_FIXTURE_BAD_CONTROL")=="1"\n'
                'text,_=fixture(2 if configuration==640 else 0,variant,bad_audit=bad)\nsys.stdout.write(text)\n')
            launcher.chmod(0o700)
            args=argparse_namespace('baseline')
            with contextlib.redirect_stdout(io.StringIO()),mock.patch.dict(os.environ,{'PATH':str(executables)+os.pathsep+os.environ['PATH']}):
                staged.record_build();staged.prepare(args)
                (base/'diagnostic_fixture.out').write_text('scheduler output created after preparation\n')
                staged.run(args);staged.collect(args)
                summary=staged.read_json(base/'campaigns/baseline/SUMMARY.json')
                self.assertTrue(summary['control_comparison']['success'])
                self.assertEqual(summary['cases'][-1]['solver_result'],'residual_nonconvergence')
                self.assertFalse(summary['production_configuration_admitted'])
                with zipfile.ZipFile(root/'qx_n640_diagnostic_baseline_review.zip') as archive:
                    self.assertIsNone(archive.testzip())
                    self.assertIn('campaigns/baseline/outputs/n640_t1_traced.err',archive.namelist())
                failed=argparse_namespace('control_failure');staged.prepare(failed)
                with mock.patch.dict(os.environ,{'QX_FIXTURE_BAD_CONTROL':'1'}):
                    with self.assertRaisesRegex(ValueError,'Control failed'):staged.run(failed)
                launches=staged.read_json(base/'campaigns/control_failure/LAUNCHES.json')
                self.assertEqual(len(launches),1)
                self.assertEqual(launches[0]['configuration'],650)
                staged.collect(failed)
                self.assertFalse(staged.read_json(base/'campaigns/control_failure/SUMMARY.json')['diagnostic_completed'])
                self.assertTrue((root/'qx_n640_diagnostic_control_failure_review.zip').is_file())
                reference.write_text('changed source fixture\n')
                with self.assertRaisesRegex(ValueError,'Production source'):staged.prepared('baseline')


def argparse_namespace(name):
    import argparse
    return argparse.Namespace(name=name)


if __name__=='__main__':unittest.main()
