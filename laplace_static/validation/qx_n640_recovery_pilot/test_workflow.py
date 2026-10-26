#!/usr/bin/env python3
"""Protocol, provenance and source-transformation checks for n640 recovery."""
import ast
import copy
import importlib.util
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

HERE = Path(__file__).resolve().parent


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


c = module('qx_recovery_controller',HERE/'campaign.py')
p = module('qx_recovery_patcher',HERE/'patch_sources.py')
SPEC=json.loads((HERE/'INPUTS.json').read_text())
TIMES=SPEC['selected_times']


def logs():
    values=[float(m['lambda']) for m in c.baseline()['audited_modes']]
    gen=['RECOVERY_START cfg=640 ranks=64 P=1,4,4,4 selected=14 res=1.00000000000000006e-09']
    audit=['QUALITY_START cfg=640 grid=1x4x4x4 Nv=100 APE=20,0.5 FUNNELED=1 smear_seconds=1']
    for t in TIMES:
        gen.append('RECOVERY_SOLVE cfg=640 t=%d code=0 restarts=1 nRed=2 maxRes=1e-11 lambda_Nv=%.17e'%(t,values[-1]))
        gen.append('RECOVERY_WRITTEN cfg=640 t=%d Nv=100'%t)
        for i,value in enumerate(values):
            audit.append('MODE cfg=640 t=%d i=%d lambda=%.17e norm=1 residual_l2=%.17e residual_relative=1e-11 rayleigh_re=%.17e rayleigh_im=0'%(t,i,value,value*1e-11,value))
        audit.append('SLICE cfg=640 t=%d Nv=100 max_relative_residual=1e-11 max_norm_defect=0 gram_diag_defect=0 gram_offdiag=1e-12 max_rayleigh_im=0 lambda_min=%.17e lambda_max=%.17e status=PASS'%(t,values[0],values[-1]))
    audit.append('QUALITY_END cfg=640 slices=14 status=PASS')
    return '\n'.join(gen)+'\n','\n'.join(audit)+'\n'


class WorkflowTests(unittest.TestCase):
    def test_real_baseline_replay(self):
        report=c.baseline()
        self.assertTrue(report['selected_slice_valid'])
        self.assertEqual(report['cycles'],4)
        self.assertLess(float(report['independent_audit']['max_relative_residual']),1e-9)
        self.assertFalse(report['production_configuration_admitted'])

    def test_all_preserved_input_hashes(self):
        for name,digest in SPEC['preserved_inputs'].items():
            self.assertEqual(c.sha(HERE/name),digest,name)

    def test_writer_patch_does_not_change_numerical_call(self):
        source=(HERE/'inputs/mental_funneled.c').read_text()
        result=p.writer(source)
        call='lret=qcd_lanczosCheby(lambda, v, &lstat, &u3d, params.Nv, &lpar);'
        self.assertEqual(source.count(call),1);self.assertEqual(result.count(call),1)
        self.assertIn('if(!qx_selected(lt))continue;',result)
        self.assertIn('continue; /* This pilot writes eigenpairs only. */',result)
        for writer in ('DistEigvecsHdf5Writer_CreateFile','DistEigvecsHdf5Writer_WriteDistEigvecs','DistEigvalsHdf5Writer_WriteDistEigvals'):
            self.assertEqual(source.count(writer),result.count(writer))
        self.assertIn('params.P[0]!=1',result)
        self.assertNotIn('geo.P[0]',result)
        with self.assertRaises(ValueError):p.writer(result)

    def test_checker_patch_strict_and_selected_only(self):
        result=p.checker((HERE/'inputs/test_saved_modes.c').read_text())
        self.assertNotIn('1.1e-9',result)
        self.assertIn('if(!qx_selected(t))continue;',result)
        self.assertIn('rel>1e-9',result)
        self.assertNotIn('solver_mpiInit',result)
        self.assertIn('DistEigvecsHdf5Reader_ReadDistEigvecsRange',result)
        self.assertIn('values[t][i]',result)
        self.assertEqual(len(TIMES),14)
        self.assertTrue(set(range(95,102)).issubset(TIMES))
        self.assertIn(1,TIMES)
        self.assertIn('192',c.input_text(Path('/tmp/new_recovery')))

    def test_valid_protocol_remains_partial(self):
        gen,audit=logs();result=c.replay_logs(gen,audit,TIMES)
        self.assertEqual(result['audited_modes'],1400)
        self.assertTrue(result['selected_saved_modes_valid'])
        self.assertFalse(result['complete_configuration'])
        self.assertFalse(result['production_configuration_admitted'])
        self.assertEqual(result['t1_baseline_max_absolute_eigenvalue_difference'],0.)

    def test_incomplete_and_duplicate_mode_rejected(self):
        gen,audit=logs()
        lines=audit.splitlines()
        index=next(i for i,line in enumerate(lines) if line.startswith('MODE '))
        for changed in (lines[:index]+lines[index+1:],lines[:index]+[lines[index],lines[index]]+lines[index+1:]):
            with self.assertRaises(ValueError):c.replay_logs(gen,'\n'.join(changed),TIMES)

    def test_residual_numerical_and_status_failures_rejected(self):
        gen,audit=logs()
        for badgen,badaudit in [
            (gen.replace('code=0','code=2',1),audit),
            (gen.replace('maxRes=1e-11','maxRes=1e-8',1),audit),
            (gen,audit.replace('residual_relative=1e-11','residual_relative=2e-9',1)),
            (gen,audit.replace('residual_relative=1e-11','residual_relative=nan',1)),
            (gen,audit.replace('norm=1 ','norm=1.1 ',1)),
            (gen,audit.replace('gram_offdiag=1e-12','gram_offdiag=1e-6',1)),
            (gen,audit.replace('status=PASS','status=FAIL',1)),
            (gen.replace('selected=14','selected=13'),audit)]:
            with self.subTest(change=(badgen[:100],badaudit[:100])):
                with self.assertRaises(ValueError):c.replay_logs(badgen,badaudit,TIMES)

    def test_baseline_eigenvalue_comparison(self):
        gen,audit=logs()
        rows=audit.splitlines()
        for i,line in enumerate(rows):
            if line.startswith('MODE cfg=640 t=1 i=0 '):
                row=c.records(line,'MODE')[0];value=float(row['lambda'])+1e-6
                rows[i]=line.replace('lambda='+row['lambda'],'lambda=%.17e'%value).replace('residual_l2='+row['residual_l2'],'residual_l2=%.17e'%(value*1e-11))
            if line.startswith('SLICE cfg=640 t=1 '):
                row=c.records(line,'SLICE')[0]
                rows[i]=line.replace('lambda_min='+row['lambda_min'],'lambda_min=%.17e'%(float(row['lambda_min'])+1e-6))
        with self.assertRaisesRegex(ValueError,'baseline'):c.replay_logs(gen,'\n'.join(rows),TIMES)

    def test_name_safety_and_runtime_exclusions(self):
        for name in ('../x','/tmp/x','x/y','x\n',''):
            with self.assertRaises(ValueError):c.directory(name)
        self.assertEqual(c.directory('selected14').name,'selected14')
        before=c.package_pins()
        with tempfile.NamedTemporaryFile(dir=str(HERE),suffix='.out') as tmp:
            self.assertEqual(before,c.package_pins())
        self.assertIn('inputs/baseline/campaigns/baseline/outputs/n640_t1_traced.out',before)

    def test_generator_sources_and_python_syntax(self):
        with tempfile.TemporaryDirectory() as tmp:
            p.generate(Path(tmp))
            self.assertIn('0,1,32,64,95,96,97,98,99,100,101,128,160,191',(Path(tmp)/'selection.h').read_text())
            self.assertTrue((Path(tmp)/'mental_recovery.c').stat().st_size>10000)
        for name in ('campaign.py','patch_sources.py'):
            ast.parse((HERE/name).read_text())
        for name in ('build.sh','submit_micro.sh'):
            subprocess.check_call(['bash','-n',str(HERE/name)])

    def test_collection_exports_failed_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/pilot';here.mkdir(parents=True)
            folder=here/'campaigns/selected14';folder.mkdir(parents=True)
            (folder/'outputs').mkdir();(folder/'outputs/generate.err').write_text('diagnostic failure\n')
            build=here/'build';build.mkdir()
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'BUILD',build),patch.object(c,'prepared',side_effect=ValueError('generator failed')):
                c.collect(type('Args',(),{'name':'selected14'})())
            summary=json.loads((folder/'SUMMARY.json').read_text())
            self.assertFalse(summary['selected_saved_modes_valid'])
            import zipfile
            with zipfile.ZipFile(root/'qx_n640_recovery_pilot_selected14_review.zip') as z:
                self.assertIn('campaign/outputs/generate.err',z.namelist())

    def test_prepare_run_collect_protocol_lifecycle(self):
        gen,audit=logs();baseline_report=c.baseline()
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/pilot';here.mkdir(parents=True)
            build=here/'build';build.mkdir();(build/'MANIFEST.json').write_text('{}\n')
            spec=copy.deepcopy(SPEC);spec['storage_root']=str(root/'saved')
            gauge={'path':str(root/'gauge'),'bytes':28991029272,'mtime_ns':123,'header_hex':'c0000000400000004000000040000000'}
            args=SimpleNamespace(name='selected14')
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'BUILD',build), \
                 patch.object(c,'check_inputs',return_value=spec),patch.object(c,'verify_build',return_value={}), \
                 patch.object(c,'gauge_record',return_value=gauge), \
                 patch.object(c,'baseline',return_value=baseline_report), \
                 patch.object(c.shutil,'disk_usage',return_value=SimpleNamespace(free=40*1024**3)):
                c.prepare(args)
                folder,plan=c.prepared(args.name)
                def fake_mpi(command,**kwargs):
                    self.assertEqual(command[:5],['mpiexec','-n','64','-ppn','8'])
                    if command[5]==str(build/'mental_recovery'):
                        kwargs['stdout'].write(gen)
                        saved=Path(plan['run_directory'])/'eigenvectors'
                        for t in TIMES:(saved/('qX_p2n640_evec_t%d.h5'%t)).write_bytes(b'x'*2048)
                        Path(plan['eigenvalue_file']).write_bytes(b'x'*2048)
                    else:
                        self.assertEqual(command[5],str(build/'check_saved'))
                        self.assertEqual(command[-1],plan['eigenvalue_file'])
                        kwargs['stdout'].write(audit)
                    return SimpleNamespace(returncode=0)
                with patch.object(c.subprocess,'run',side_effect=fake_mpi):c.run(args)
                c.collect(args)
                summary=json.loads((folder/'SUMMARY.json').read_text())
                self.assertTrue(summary['selected_saved_modes_valid'])
                self.assertFalse(summary['complete_configuration'])
                self.assertEqual(len(json.loads((folder/'SAVED_FILES.json').read_text())),15)
                with self.assertRaises(ValueError):c.prepare(args)
                (folder/'mental.in').write_text('changed\n')
                with self.assertRaisesRegex(ValueError,'Input changed'):c.prepared(args.name)


if __name__=='__main__':unittest.main()
