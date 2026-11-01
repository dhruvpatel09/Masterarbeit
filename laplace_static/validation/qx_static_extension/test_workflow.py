#!/usr/bin/env python3
"""Coverage, contraction and provenance checks for the temporal/cutoff scan."""
import ast
import json
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import campaign as c
import check_extension as check
import generate_sources as helper
from measurement_plan import (AUDIT_TIMES, BLOCK_RADII, BLOCK_TAUS, CUTOFFS,
                              RADII, cutoffs, direct_keys, expected_keys)


def line(prefix,**values):
    return prefix+' '+' '.join('%s=%s'%(k,v) for k,v in values.items())+'\n'


def fixture(cfg,piter,refs):
    output=[line('EXTENSION_BEGIN',cfg=cfg,piter=piter,saved_Nv=100,
        cutoffs=','.join(map(str,cutoffs(piter))),ranks=64,grid='1x4x4x4',t0=90,
        tau_min=0,tau_max=32,radii=','.join(map(str,RADII)),temporal_bc='open',
        hyp_steps=1,temporal_change=.1,preserved_change=0)]
    for scheme,ns,taus,rs in [('regression',[100],(2,4,6),(1,2,4)),
                               ('fixed',cutoffs(piter),range(33),RADII)]:
        for n in ns:
            for tau in taus:
                for r in rs:
                    real=[];imag=[]
                    for axis in (0,1,2):
                        if n==100 and (scheme=='regression' or tau<=16):
                            ref=refs[cfg,piter,scheme,n,tau,r,axis]
                            x,y=float(ref['L_RE']),float(ref['L_IM'])
                        else:
                            x=1e-5*(n/100.)**2*math.exp(-.2*tau)/(1+.01*r)*(1+.001*axis)
                            y=0 if r==32 else 1e-8*(axis-1)
                        real.append(x);imag.append(y)
                        output.append(line('RT',cfg=cfg,piter=piter,scheme=scheme,
                            t0=90 if scheme=='fixed' else 98-tau//2,tau=tau,R=r,
                            axis=axis,Nv=n,Nsrc=64**3,L_RE=x,L_IM=y))
                        if scheme=='fixed' and n==100 and tau in BLOCK_TAUS and r in BLOCK_RADII:
                            for rank in range(64):
                                output.append(line('RT_BLOCK',cfg=cfg,piter=piter,Nv=100,t0=90,
                                    tau=tau,R=r,axis=axis,rank=rank,px=rank//16,
                                    py=(rank//4)%4,pz=rank%4,Nsrc=4096,L_RE=x,L_IM=y))
                    output.append(line('RT',cfg=cfg,piter=piter,scheme=scheme,
                        t0=90 if scheme=='fixed' else 98-tau//2,tau=tau,R=r,axis=-1,
                        Nv=n,Nsrc=3*64**3,L_RE=math.fsum(real)/3,L_IM=math.fsum(imag)/3))
    for scheme,n,t,r,a,s in sorted(direct_keys(piter)):
        output.append(line('DIRECT_TEST',scheme=scheme,Nv=n,tau=t,R=r,axis=a,
            site=s,abs_error=0,tolerance=1e-16,status='PASS'))
    for axis in range(3):output.append(line('PERIODIC_ROUNDTRIP',axis=axis,cells=64,status='PASS'))
    for t in range(90,122):
        output.append(line('SU3_TIME',cfg=cfg,piter=piter,t=t,count=64**3,
            unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,
            unit_above_1em2=0,det_above_1em2=0))
    output.append(line('SU3_WINDOW',cfg=cfg,piter=piter,t0=90,t1=122,count=32*64**3,
        unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,
        unit_above_1em2=0,det_above_1em2=0))
    output.append(line('EXTENSION_END',cfg=cfg,piter=piter,records=len(expected_keys(piter)),status='PASS'))
    return ''.join(output)


def audit_fixture(cfg):
    output=[]
    for t in AUDIT_TIMES:
        for i in range(100):
            output.append(line('MODE',cfg=cfg,t=t,i=i,**{'lambda':.01+.001*i,
                'norm':1,'residual_l2':1e-13,'residual_relative':5e-12,
                'rayleigh_re':.01+.001*i,'rayleigh_im':0}))
        output.append(line('SLICE',cfg=cfg,t=t,Nv=100,max_relative_residual=5e-12,
            max_norm_defect=1e-15,gram_diag_defect=1e-15,gram_offdiag=1e-15,
            max_rayleigh_im=0,lambda_min=.01,lambda_max=.109,status='PASS'))
    return ''.join(output)+line('QUALITY_END',cfg=cfg,slices=33,status='PASS')


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs=c.references()
        cls.scans={p:fixture(650,p,cls.refs) for p in (10,40)}

    def test_all71_inventory_and_pilot_reuse(self):
        cases,slots=c.case_plan('full',c.PILOT,2)
        self.assertEqual(sorted([r['configuration'] for r in cases]+c.PILOT),list(range(100,801,10)))
        self.assertEqual((len(cases),len(slots)),(68,34))
        self.assertEqual([r['configuration'] for r in c.case_plan('pilot')[0]],[170,640,650])
        self.assertEqual(len(self.refs),63048)
        for cfg in c.TARGET:
            self.assertEqual(sum(k[0]==cfg and k[1]==10 for k in self.refs),852)
        self.assertEqual((len(expected_keys(10)),len(expected_keys(40))),(6372,1620))
        self.assertEqual((len(direct_keys(10)),len(direct_keys(40))),(396,108))

    def test_complete_scans_and_projection_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out';parsed={}
            for piter in (10,40):
                p.write_text(self.scans[piter]);values,blocks,summary=check.scan(p,650,piter,self.refs)
                parsed[piter]=values
                self.assertEqual((len(values),len(blocks)),(len(expected_keys(piter)),4800))
                self.assertEqual(summary['archived_reference_rows'],852)
                self.assertEqual(summary['cutoffs'],cutoffs(piter))
            self.assertLess(check.projection_pair(parsed[10],parsed[40]),1e-19)
            changed=[dict(r) for r in parsed[40]]
            row=next(r for r in changed if r['scheme']=='fixed' and r['tau']=='32')
            row['L_RE']=str(float(row['L_RE'])*1.01)
            with self.assertRaises(ValueError):check.projection_pair(parsed[10],changed)

    def test_missing_duplicate_and_cutoff_rows_rejected(self):
        original=self.scans[10]
        first=next(s for s in original.splitlines(True) if s.startswith('RT ') and 'Nv=40 ' in s)
        block=next(s for s in original.splitlines(True) if s.startswith('RT_BLOCK '))
        direct=next(s for s in original.splitlines(True) if s.startswith('DIRECT_TEST ') and 'Nv=80 ' in s)
        variants=[original.replace(first,'',1),original+first,original.replace('cutoffs=40,60,80,100','cutoffs=40,60,100',1),
            original.replace(block,'',1),original.replace(direct,'',1),original.replace('unit_max=3e-16','unit_max=3e-3',1)]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out'
            for variant in variants:
                p.write_text(variant)
                with self.assertRaises(ValueError):check.scan(p,650,10,self.refs)

    def test_overlapping_baseline_must_match(self):
        refs=dict(self.refs);key=(650,10,'fixed',100,16,24,0)
        refs[key]=dict(refs[key],L_RE='1')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out';p.write_text(self.scans[10])
            with self.assertRaises(ValueError):check.scan(p,650,10,refs)

    def test_strict_audit_including_order_and_new_times(self):
        text=audit_fixture(640)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'audit.out';p.write_text(text)
            self.assertEqual(check.audit(p,640)['modes'],3300)
            first=next(s for s in text.splitlines(True) if s.startswith('MODE ') and 't=122 ' in s)
            variants=[text.replace('residual_relative=5e-12','residual_relative=1.01e-9',1),
                text.replace(first,'',1),text.replace('lambda=0.01 norm=1','lambda=0.5 norm=1',1)]
            for variant in variants:
                p.write_text(variant)
                with self.assertRaises(ValueError):check.audit(p,640)

    def test_case_and_log_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp);c.write(target/'PREPARED.json',{'test':True})
            entry=dict(index=0,configuration=650,piters=[10,40])
            case=target/'case_000';case.mkdir();out=case/'attempt_0001';out.mkdir()
            c.write(out/'CASE.json',dict(entry,prepared_sha256=c.digest(target/'PREPARED.json')))
            launches=[]
            for kind,piter,text in [('audit',0,audit_fixture(650)),('scan',10,self.scans[10]),('scan',40,self.scans[40])]:
                stem='audit' if not piter else 'p%d'%piter
                (out/(stem+'.out')).write_text(text);(out/(stem+'.err')).write_text('')
                launches.append(dict(kind=kind,piter=piter,stem=stem,exit_code=0,seconds=1,
                    out_sha256=c.digest(out/(stem+'.out')),err_sha256=c.digest(out/(stem+'.err'))))
            c.write(out/'LAUNCHES.json',launches)
            values,blocks,summary=c.replay_case(out,entry,self.refs)
            self.assertEqual((len(values),len(blocks)),(7992,9600))
            c.write(out/'SUMMARY.json',dict(success=True,launches_sha256=c.digest(out/'LAUNCHES.json'),case_sha256=c.digest(out/'CASE.json')))
            c.write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=c.digest(out/'SUMMARY.json')))
            self.assertEqual(c.complete_case(target,entry),out)
            (out/'p40.err').write_text('modified')
            with self.assertRaises(ValueError):c.replay_case(out,entry,self.refs)
            (out/'LAUNCHES.json').write_text('[{}]')
            with self.assertRaises(ValueError):c.complete_case(target,entry)

    def test_python36_syntax_and_source_pins(self):
        for path in c.HERE.glob('*.py'):
            options=dict(feature_version=(3,6)) if sys.version_info>=(3,8) else {}
            ast.parse(path.read_text(),filename=str(path),**options)
        c.verify_sources()
        inputs=json.loads((c.HERE/'INPUTS.json').read_text())
        self.assertEqual(inputs['geometry']['temporal_separations'],list(range(33)))
        self.assertEqual(inputs['geometry']['mode_cutoffs'],CUTOFFS)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is needed for the native control')
    def test_native_cutoffs_and_periodic_shifts(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp);exe=Path(tmp)/'kernels'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(c.HERE/'tests'),'-I'+tmp,
                str(c.HERE/'tests/kernels_and_shifts.c'),'-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is needed for the syntax control')
    def test_scan_syntax_and_unchanged_payload_io(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp)
            include=[c.HERE/'tests',Path(tmp),c.ROOT/'solver/include',
                c.ROOT/'laplace_static/validation/qx_n640_diagnostic/inputs/qcd_include',
                c.ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230']
            subprocess.run(['gcc','-std=c11','-fsyntax-only','-Werror=implicit-function-declaration',
                '-Werror=format']+['-I'+str(p) for p in include]+[str(c.HERE/'scan.c')],check=True)
            audit=(Path(tmp)/'audit_selected.c').read_text()
            self.assertIn('static const int selected_times[]={'+','.join(map(str,AUDIT_TIMES))+'};',audit)
            self.assertNotIn('1.1e-9',audit)
            with tempfile.TemporaryDirectory() as old:
                helper.parent.generate(old)
                for name in ('io_helpers.h','mpi_helpers.h'):
                    self.assertEqual((Path(tmp)/name).read_bytes(),(Path(old)/name).read_bytes())


if __name__=='__main__':unittest.main()
