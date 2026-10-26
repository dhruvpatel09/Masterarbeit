#!/usr/bin/env python3
"""Coverage, numerical-protocol, native selection and assembly-parser checks."""
import ast
import contextlib
import copy
import importlib.util
import io
import json
import os
import shutil
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


c=load('qx_full_controller',HERE/'campaign.py')
p=load('qx_full_sources',HERE/'patch_sources.py')
SPEC=json.loads((HERE/'INPUTS.json').read_text())
PILOT=HERE/'inputs/pilot/campaign'


def fixtures():
    values,_=c.validate_audit((PILOT/'outputs/audit.out').read_text(),SPEC['reused_times'])
    for t in range(192):
        if t not in SPEC['reused_times']:
            for i in range(100):values[(t,i)]=values[(98,i)]
    return values


VALUES=fixtures()


def fake_audit(times):
    rows=['QUALITY_START cfg=640 grid=1x4x4x4 Nv=100 APE=20,0.5 FUNNELED=1 smear_seconds=1']
    for t in times:
        for i in range(100):
            value=VALUES[(t,i)]
            rows.append('MODE cfg=640 t=%d i=%d lambda=%.17e norm=1 residual_l2=%.17e residual_relative=1e-11 rayleigh_re=%.17e rayleigh_im=0'%(t,i,value,value*1e-11,value))
        rows.append('SLICE cfg=640 t=%d Nv=100 max_relative_residual=1e-11 max_norm_defect=0 gram_diag_defect=0 gram_offdiag=1e-12 max_rayleigh_im=0 lambda_min=%.17e lambda_max=%.17e status=PASS'%(t,VALUES[(t,0)],VALUES[(t,99)]))
    rows.append('QUALITY_END cfg=640 slices=%d status=PASS'%len(times))
    return '\n'.join(rows)+'\n'


def fake_generation(times):
    rows=['RECOVERY_START cfg=640 ranks=64 P=1,4,4,4 selected=%d res=1.00000000000000006e-09'%len(times)]
    for t in times:
        rows.append('RECOVERY_SOLVE cfg=640 t=%d code=0 restarts=0 nRed=1 maxRes=1e-11 lambda_Nv=%.17e'%(t,VALUES[(t,99)]))
        rows.append('RECOVERY_WRITTEN cfg=640 t=%d Nv=100'%t)
    return '\n'.join(rows)+'\n'


class Tests(unittest.TestCase):
    def setUp(self):
        output=contextlib.redirect_stdout(io.StringIO());output.__enter__()
        self.addCleanup(output.__exit__,None,None,None)

    def test_actual_pilot_replay(self):
        plan,summary=c.pilot_replay()
        self.assertTrue(summary['selected_saved_modes_valid'])
        self.assertEqual(summary['audited_modes'],1400)
        self.assertEqual(summary['t1_baseline_max_absolute_eigenvalue_difference'],0)
        self.assertEqual(plan,SPEC['pilot_plan'])

    def test_exact_disjoint_coverage(self):
        c.check_spec(SPEC)
        todo=sum(SPEC['slots'],[])
        self.assertEqual(len(todo),178);self.assertEqual(len(set(todo)),178)
        self.assertFalse(set(todo)&set(SPEC['reused_times']))
        self.assertEqual(sorted(todo+SPEC['reused_times']),list(range(192)))
        self.assertEqual(len(SPEC['slots']),16)
        self.assertLessEqual(max(map(len,SPEC['slots'])),12)
        changed=copy.deepcopy(SPEC);changed['slots'][0].append(1)
        with self.assertRaises(ValueError):c.check_spec(changed)

    def test_preserved_input_hashes(self):
        for name,digest in SPEC['preserved_inputs'].items():self.assertEqual(c.sha(HERE/name),digest,name)

    def test_source_transform_is_selection_only(self):
        source=(HERE/'inputs/mental_recovery.c').read_text();result=p.transform_writer(source)
        call='lret=qcd_lanczosCheby(lambda, v, &lstat, &u3d, params.Nv, &lpar);'
        self.assertEqual(source.count(call),1);self.assertEqual(result.count(call),1)
        self.assertEqual(source.count('DistEigvecsHdf5Writer_WriteDistEigvecs'),result.count('DistEigvecsHdf5Writer_WriteDistEigvecs'))
        self.assertIn('qx_slot()<0',result)
        with self.assertRaises(ValueError):p.transform_writer(result)
        checker=p.transform_checker((HERE/'inputs/check_saved.c').read_text())
        self.assertIn('static int selected(int t) {return qx_selected(t);}',checker)
        self.assertIn('require(slices!=qx_selected_count()',checker)
        self.assertNotIn('1.1e-9',checker)

    @unittest.skipUnless(shutil.which('gcc'),'native gcc unavailable')
    def test_native_selection_all_slots_and_final(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'selection.h').write_text(p.header(SPEC['reused_times']))
            (d/'probe.c').write_text('#include "selection.h"\nint main(void){printf("%d\\n",qx_selected_count());for(int t=0;t<192;++t)if(qx_selected(t))printf("%d ",t);return 0;}\n')
            subprocess.run(['gcc','-std=c11','-Wall','-Wextra',str(d/'probe.c'),'-o',str(d/'probe')],check=True)
            for index,times in list(enumerate(SPEC['slots']))+[('all',list(range(192)))]:
                result=subprocess.check_output([str(d/'probe')],env=dict(os.environ,QX_RECOVERY_SLOT=str(index))).decode().split()
                self.assertEqual(int(result[0]),len(times));self.assertEqual(list(map(int,result[1:])),times)
            for bad in ('','-1','16','1x'):
                r=subprocess.run([str(d/'probe')],env=dict(os.environ,QX_RECOVERY_SLOT=bad),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                self.assertNotEqual(r.returncode,0)

    @unittest.skipUnless(shutil.which('gcc'),'native gcc unavailable')
    def test_native_assembly_values_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);shutil.copy2(HERE/'audited_values.h',d/'audited_values.h')
            (d/'probe.c').write_text('#include "audited_values.h"\nstatic double v[192][100];int main(int n,char **a){return n!=2 || read_audited_values(a[1],v);}\n')
            subprocess.run(['gcc','-std=c11','-Wall','-Wextra',str(d/'probe.c'),'-lm','-o',str(d/'probe')],check=True)
            lines=['t,i,lambda']+['%d,%d,%.17e'%(t,i,VALUES[(t,i)]) for t in range(192) for i in range(100)]
            good='\n'.join(lines)+'\n';data=d/'values.csv';data.write_text(good)
            self.assertEqual(subprocess.run([str(d/'probe'),str(data)]).returncode,0)
            for bad in ('\n'.join(lines[:-1])+'\n',good+'0,0,0.03\n',good.replace('0,0,','0,1,',1),
                        good.replace(lines[1],'0,0,nan',1),good.replace(lines[1],'0,0,0.15',1)):
                data.write_text(bad)
                self.assertNotEqual(subprocess.run([str(d/'probe'),str(data)]).returncode,0)

    def test_full_audit_and_eigenvalue_comparison(self):
        text=fake_audit(list(range(192)))
        values,slices=c.validate_audit(text,list(range(192)),VALUES)
        self.assertEqual(len(values),19200);self.assertEqual(len(slices),192)
        changed=dict(VALUES);changed[(37,3)]+=1e-8
        with self.assertRaisesRegex(ValueError,'eigenvalue mismatch'):c.validate_audit(text,list(range(192)),changed)

    def test_audit_rejects_missing_bad_nonfinite_and_gram(self):
        times=SPEC['slots'][0];text=fake_audit(times)
        for changed in (text.replace('residual_relative=1e-11','residual_relative=1.01e-9',1),
                        text.replace('residual_relative=1e-11','residual_relative=nan',1),
                        text.replace('norm=1 ','norm=1.01 ',1),text.replace('gram_offdiag=1e-12','gram_offdiag=1e-6',1),
                        text.replace('status=PASS','status=FAIL',1),'\n'.join(text.splitlines()[2:])):
            with self.assertRaises(ValueError):c.validate_audit(changed,times)

    def test_generation_requires_all_times_and_original_target(self):
        times=SPEC['slots'][0];text=fake_generation(times)
        self.assertEqual(len(c.validate_generation(text,times)),len(times))
        for changed in (text.replace('code=0','code=2',1),text.replace('maxRes=1e-11','maxRes=1.1e-9',1),
                        text.replace('res=1.00000000000000006e-09','res=1e-8'),text.replace('RECOVERY_WRITTEN','BAD',1)):
            with self.assertRaises(ValueError):c.validate_generation(changed,times)

    def test_submission_dependency_and_duplicate_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'PREPARED.json').write_text('{}\n')
            args=SimpleNamespace(name='full192',concurrent=4)
            with patch.object(c,'prepared',return_value=(folder,{})),patch.object(c,'sbatch',side_effect=['100','101']) as calls:
                c.submit(args)
                self.assertEqual(calls.call_args_list[0][0][0][2],'--array=0-15%4')
                self.assertEqual(calls.call_args_list[1][0][0][2],'--dependency=afterok:100')
                self.assertEqual(c.read(folder/'SUBMISSION.json')['final_job'],'101')
                with self.assertRaises(ValueError):c.submit(args)

    def test_saved_inventory_change_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            pth=Path(tmp)/'mode';pth.write_bytes(b'initial')
            record=c.stat_record(pth);c.check_saved_inventory([record]);pth.write_bytes(b'changed length')
            with self.assertRaises(ValueError):c.check_saved_inventory([record])

    def test_prepare_and_input_immutability(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/full';here.mkdir(parents=True)
            shutil.copytree(HERE/'inputs',here/'inputs')
            build=here/'build';build.mkdir();(build/'MANIFEST.json').write_text('{}\n')
            gauge=root/'gauge'
            with gauge.open('wb') as f:
                f.write(struct.pack('<4i',192,64,64,64));f.truncate(28991029272)
            old=c.pilot_module();spec=copy.deepcopy(SPEC)
            spec['storage_root']=str(root/'saved')
            spec['pilot_plan']['gauge']=old.gauge_record(gauge)
            saved=root/'pilot';(saved/'eigenvectors').mkdir(parents=True)
            inventory=[]
            for t in spec['reused_times']:
                file=saved/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t);file.write_bytes(b'x'*2048)
                inventory.append(c.stat_record(file))
            file=saved/'values.h5';file.write_bytes(b'e'*2048);inventory.append(c.stat_record(file))
            spec['pilot_plan']['run_directory']=str(saved);spec['pilot_plan']['eigenvalue_file']=str(file)
            old_read=c.read
            def read_fixture(path):
                if Path(path)==here/'inputs/pilot/campaign/SAVED_FILES.json':return inventory
                return old_read(path)
            args=SimpleNamespace(name='full192')
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'BUILD',build), \
                 patch.object(c,'check_inputs',return_value=spec),patch.object(c,'verify_build'),patch.object(c,'read',side_effect=read_fixture), \
                 patch.object(c.shutil,'disk_usage',return_value=SimpleNamespace(free=300*1024**3)):
                c.prepare(args);folder,plan=c.prepared('full192')
                self.assertEqual(sum(len(x['times']) for x in plan['slots']),178)
                self.assertEqual(len(plan['slots']),16)
                with self.assertRaises(ValueError):c.prepare(args)
                (folder/'slot_00/mental.in').write_text('changed\n')
                with self.assertRaisesRegex(ValueError,'Slot input changed'):c.prepared('full192')

    def test_protocol_lifecycle_all_slots_assembly_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/full';here.mkdir(parents=True)
            shutil.copytree(HERE/'inputs',here/'inputs')
            build=here/'build';build.mkdir();(build/'MANIFEST.json').write_text('{}\n')
            folder=here/'campaigns/full192';folder.mkdir(parents=True)
            run=root/'saved';run.mkdir();pilot=run/'pilot';(pilot/'eigenvectors').mkdir(parents=True)
            for t in SPEC['reused_times']:(pilot/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)).write_bytes(b'p'*2048)
            pvalues=pilot/'pilot_eval.h5';pvalues.write_bytes(b'e'*2048)
            plan={'name':'full192','run_directory':str(run),'gauge':{'path':str(root/'gauge')},'reused_times':SPEC['reused_times'],
                  'pilot_plan':{'run_directory':str(pilot),'eigenvalue_file':str(pvalues)},'slots':[]}
            for index,times in enumerate(SPEC['slots']):
                saved=run/('slot_%02d'%index);(saved/'eigenvectors').mkdir(parents=True);(saved/'eigenvalues').mkdir()
                slot=folder/('slot_%02d'%index);slot.mkdir();(slot/'mental.in').write_text('fixture input\n')
                plan['slots'].append({'index':index,'times':times,'run_directory':str(saved),'eigenvalue_file':str(saved/'eigenvalues/values.h5')})
            def fake_mpi(command,**kwargs):
                exe=Path(command[5]).name;slot=kwargs['env']['QX_RECOVERY_SLOT']
                if exe=='mental_recovery':
                    row=plan['slots'][int(slot)];kwargs['stdout'].write(fake_generation(row['times']))
                    for t in row['times']:(Path(row['run_directory'])/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)).write_bytes(b'x'*2048)
                    Path(row['eigenvalue_file']).write_bytes(b'e'*2048)
                elif exe=='check_saved':
                    times=list(range(192)) if slot=='all' else plan['slots'][int(slot)]['times']
                    kwargs['stdout'].write(fake_audit(times))
                else:
                    self.assertEqual(exe,'assemble_values');self.assertEqual(command[1:5],['-n','1','-ppn','1'])
                    kwargs['stdout'].write('ASSEMBLED_VALUES times=192 Nv=100 rows=19200\n')
                return SimpleNamespace(returncode=0)
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'BUILD',build), \
                 patch.object(c,'prepared',return_value=(folder,plan)),patch.object(c.subprocess,'run',side_effect=fake_mpi):
                for index in range(16):c.run_slot(SimpleNamespace(name='full192',slot=index))
                c.run_final(SimpleNamespace(name='full192'))
                summary=c.final_replay(folder,plan)
                self.assertTrue(summary['complete_eigenpair_set_valid'])
                self.assertFalse(summary['production_configuration_admitted'])
                self.assertEqual(summary['audited_modes'],19200)
                self.assertEqual(len(list((run/'assembled/eigenvectors').glob('*.h5'))),192)
                c.collect(SimpleNamespace(name='full192'))
                self.assertTrue((root/'qx_n640_recovery_full192_review.zip').is_file())

    def test_failures_remain_exportable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);here=root/'laplace_static/validation/full';here.mkdir(parents=True)
            folder=here/'campaigns/full192';folder.mkdir(parents=True)
            build=here/'build';build.mkdir();(folder/'failure.err').write_text('failed MPI attempt\n')
            with patch.object(c,'HERE',here),patch.object(c,'ROOT',root),patch.object(c,'BUILD',build), \
                 patch.object(c,'prepared',side_effect=ValueError('slot failed')):
                c.collect(SimpleNamespace(name='full192'))
            self.assertFalse(c.read(folder/'SUMMARY.json')['complete_eigenpair_set_valid'])
            import zipfile
            with zipfile.ZipFile(root/'qx_n640_recovery_full192_review.zip') as z:self.assertIn('campaign/failure.err',z.namelist())

    def test_python_shell_sources_and_name_safety(self):
        for name in ('campaign.py','patch_sources.py'):ast.parse((HERE/name).read_text())
        for name in ('build.sh','submit_micro.sh','submit_final.sh'):subprocess.run(['bash','-n',str(HERE/name)],check=True)
        for name in ('../x','x/y','x\n',''):
            with self.assertRaises(ValueError):c.directory(name)
        with tempfile.TemporaryDirectory() as tmp:p.generate(Path(tmp))


if __name__=='__main__':unittest.main()
