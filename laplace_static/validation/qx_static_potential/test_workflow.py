#!/usr/bin/env python3
"""Coverage, physics-kernel and distributed-shift controls for the static scan."""
import json
import math
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import campaign as c
import check_outputs as check
import generate_helpers as helper


def line(prefix,**values):return prefix+' '+' '.join('%s=%s'%(k,v) for k,v in values.items())+'\n'


def fixture(cfg,piter,refs):
    text=line('STATIC_BEGIN',cfg=cfg,piter=piter,Nv=100,ranks=64,grid='1x4x4x4',t0=90,tau_min=0,tau_max=16,R_min=1,R_max=32,temporal_bc='open',hyp_steps=1,temporal_change=.1,preserved_change=0)
    for scheme,taus,rs in [('regression',(2,4,6),(1,2,4)),('fixed',range(17),range(1,33))]:
        for tau in taus:
            for r in rs:
                real=[];imag=[]
                for axis in (0,1,2):
                    if scheme=='regression':
                        ref=refs[cfg,piter,tau,r,axis];x=float(ref['L_RE']);y=float(ref['L_IM'])
                    else:x=1e-5*math.exp(-.2*tau)/(1+.01*r)*(1+.001*axis);y=0 if r==32 else 1e-8*(axis-1)
                    real.append(x);imag.append(y)
                    text+=line('RT',cfg=cfg,piter=piter,scheme=scheme,t0=90 if scheme=='fixed' else 98-tau//2,tau=tau,R=r,axis=axis,Nv=100,Nsrc=64**3,L_RE=x,L_IM=y)
                    if scheme=='fixed' and tau in (0,4,8,16) and r in (1,4,8,16,32):
                        for rank in range(64):text+=line('RT_BLOCK',cfg=cfg,piter=piter,t0=90,tau=tau,R=r,axis=axis,rank=rank,px=rank//16,py=(rank//4)%4,pz=rank%4,Nsrc=4096,L_RE=x,L_IM=y)
                    if tau in (0,4,16) and r in (1,4,16,17,31,32):
                        for site in (0,4095):text+=line('DIRECT_TEST',scheme=scheme,tau=tau,R=r,axis=axis,site=site,abs_error=0,tolerance=1e-16,status='PASS')
                if scheme=='regression':
                    ref=refs[cfg,piter,tau,r,-1];x=float(ref['L_RE']);y=float(ref['L_IM'])
                else:x=math.fsum(real)/3;y=math.fsum(imag)/3
                text+=line('RT',cfg=cfg,piter=piter,scheme=scheme,t0=90 if scheme=='fixed' else 98-tau//2,tau=tau,R=r,axis=-1,Nv=100,Nsrc=3*64**3,L_RE=x,L_IM=y)
    for axis in range(3):text+=line('PERIODIC_ROUNDTRIP',axis=axis,cells=64,status='PASS')
    for t in range(90,106):text+=line('SU3_TIME',cfg=cfg,piter=piter,t=t,count=64**3,unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,unit_above_1em2=0,det_above_1em2=0)
    text+=line('SU3_WINDOW',cfg=cfg,piter=piter,t0=90,t1=106,count=16*64**3,unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,unit_above_1em2=0,det_above_1em2=0)
    text+=line('STATIC_END',cfg=cfg,piter=piter,records=2212,status='PASS')
    return text


def audit_fixture(cfg):
    text=''
    for t in range(90,107):
        for i in range(100):text+=line('MODE',cfg=cfg,t=t,i=i,**{'lambda':.01+.001*i,'norm':1,'residual_l2':1e-13,'residual_relative':5e-12,'rayleigh_re':.01+.001*i,'rayleigh_im':0})
        text+=line('SLICE',cfg=cfg,t=t,Nv=100,max_relative_residual=5e-12,max_norm_defect=1e-15,gram_diag_defect=1e-15,gram_offdiag=1e-15,max_rayleigh_im=0,lambda_min=.01,lambda_max=.109,status='PASS')
    return text+line('QUALITY_END',cfg=cfg,slices=17,status='PASS')


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.refs=c.references();cls.scan=fixture(650,10,cls.refs)

    def test_full71_plan_reuses_pilot(self):
        cases,slots=c.case_plan('full',c.PILOT,2)
        self.assertEqual(sorted([r['configuration'] for r in cases]+c.PILOT),list(range(100,801,10)))
        self.assertEqual(len(cases),68);self.assertEqual(len(slots),34)
        self.assertEqual([r['configuration'] for r in c.case_plan('pilot')[0]],[170,640,650])
        self.assertEqual(len(self.refs),5112)

    def test_source_parser_handles_comment_braces(self):
        source='static int f(void) { /* } */ const char *s="}"; // }\n return 0; }\nstatic int next(void) {return 1;}'
        result=helper.function(source,'f')
        self.assertIn('return 0',result);self.assertNotIn('next',result)

    def test_complete_scan_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out';p.write_text(self.scan)
            values,blocks,summary=check.scan(p,650,10,self.refs)
            self.assertEqual((len(values),len(blocks)),(2212,3840));self.assertEqual(summary['direct_tests'],120)

    def test_missing_duplicate_or_changed_geometry_rejected(self):
        first=next(x for x in self.scan.splitlines(True) if x.startswith('RT '))
        block=next(x for x in self.scan.splitlines(True) if x.startswith('RT_BLOCK '))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out'
            for corrupt in (self.scan.replace(first,'',1),self.scan+first,self.scan.replace('t0=90 tau_min=0','t0=89 tau_min=0',1),self.scan.replace(block,'',1),self.scan.replace('unit_max=3e-16','unit_max=3e-3',1)):
                p.write_text(corrupt)
                with self.assertRaises(ValueError):check.scan(p,650,10,self.refs)

    def test_regression_and_projection_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'scan.out';p.write_text(self.scan);values,_,_=check.scan(p,650,10,self.refs)
            changed=[dict(r) for r in values]
            changed[0]['L_RE']=str(float(changed[0]['L_RE'])*1.001)
            with self.assertRaises(ValueError):check.projection_pair(values,changed)
            self.assertEqual(check.projection_pair(values,values),0)
            ref=dict(self.refs);key=(650,10,2,1,0);ref[key]=dict(ref[key],L_RE='1')
            with self.assertRaises(ValueError):check.scan(p,650,10,ref)

    def test_saved_mode_coverage_and_strict_residual(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'audit.out';text=audit_fixture(640);p.write_text(text)
            self.assertEqual(check.audit(p,640)['modes'],1700)
            p.write_text(text.replace('residual_relative=5e-12','residual_relative=1.01e-9',1))
            with self.assertRaises(ValueError):check.audit(p,640)

    def test_completed_case_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp);c.write(target/'PREPARED.json',{'test':True});entry=dict(index=0,configuration=640,piters=[10])
            case=target/'case_000';case.mkdir();out=case/'attempt_0001';out.mkdir()
            c.write(out/'CASE.json',dict(entry,prepared_sha256=c.digest(target/'PREPARED.json')))
            c.write(out/'LAUNCHES.json',[])
            c.write(out/'SUMMARY.json',dict(success=True,launches_sha256=c.digest(out/'LAUNCHES.json'),case_sha256=c.digest(out/'CASE.json')))
            c.write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=c.digest(out/'SUMMARY.json')))
            self.assertEqual(c.complete_case(target,entry),out)
            (out/'LAUNCHES.json').write_text('[{}]')
            with self.assertRaises(ValueError):c.complete_case(target,entry)

    def test_replay_requires_complete_logs_and_projection_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);launches=[];entry=dict(configuration=650,piters=[10,40])
            for kind,piter,text in [('audit',0,audit_fixture(650)),('scan',10,self.scan),('scan',40,fixture(650,40,self.refs))]:
                stem='audit' if not piter else 'p%d'%piter
                (out/(stem+'.out')).write_text(text);(out/(stem+'.err')).write_text('')
                launches.append(dict(kind=kind,piter=piter,stem=stem,exit_code=0,seconds=1,
                    out_sha256=c.digest(out/(stem+'.out')),err_sha256=c.digest(out/(stem+'.err'))))
            c.write(out/'LAUNCHES.json',launches)
            values,blocks,summary=c.replay_case(out,entry,self.refs)
            self.assertEqual((len(values),len(blocks)),(4424,7680));self.assertEqual(summary['quality']['slices'],17)
            (out/'p40.out').write_text(self.scan)
            with self.assertRaises(ValueError):c.replay_case(out,entry,self.refs)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is needed for the native kernel control')
    def test_native_contractions_and_periodic_shifts(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp);exe=Path(tmp)/'kernels'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(c.HERE/'tests'),'-I'+tmp,str(c.HERE/'tests/kernels_and_shifts.c'),'-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is needed for the syntax control')
    def test_scan_syntax_against_preserved_qcd_interfaces(self):
        with tempfile.TemporaryDirectory() as tmp:
            helper.generate(tmp)
            include=[c.HERE/'tests',Path(tmp),c.ROOT/'solver/include',c.ROOT/'laplace_static/validation/qx_n640_diagnostic/inputs/qcd_include',c.ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230']
            subprocess.run(['gcc','-std=c11','-fsyntax-only','-Werror=implicit-function-declaration']+['-I'+str(p) for p in include]+[str(c.HERE/'scan.c')],check=True)
            audit=(Path(tmp)/'audit_selected.c').read_text()
            self.assertIn('static const int selected_times[]={'+','.join(map(str,range(90,107)))+'};',audit)
            self.assertNotIn('1.1e-9',audit)


if __name__=='__main__':unittest.main()
