"""Analytic parser fixtures only; these are not gauge-field measurements."""
import csv
import math
from matrix_plan import (PILOT,CUTOFFS,RADII,TAUS,AXES,AUDIT_TIMES,
                         expected_keys,direct_keys,legacy_keys,block_keys)


def line(prefix,**values):
    return prefix+' '+' '.join('%s=%s'%(key,value) for key,value in values.items())+'\n'


def value(cfg,source,sink,tau,radius,axis):
    if axis==-1:
        return sum(value(cfg,source,sink,tau,radius,a) for a in (0,1,2))/3
    scale=(1+.0001*cfg)*(1+.1*axis)/(1+.01*radius)
    if tau==0:
        # A strictly positive Gram matrix of nested real vectors.
        return complex(scale*(min(source,sink)*1e-8+source*sink*2e-10))
    real=scale*source*sink*1e-10*math.exp(-.2*tau)*math.cos(.08*tau)*(1+.003*tau*(source-sink)/100)
    imag=0 if radius==32 else scale*source*sink*1e-12*math.exp(-.2*tau)*(source-2*sink)/100
    return complex(real,imag)


def references(configs=PILOT):
    result={}
    for cfg in configs:
        for n in CUTOFFS:
            for tau in TAUS:
                for r in RADII:
                    for a in AXES:
                        z=value(cfg,n,n,tau,r,a)
                        result[cfg,n,tau,r,a]=dict(cfg=str(cfg),piter='10',scheme='fixed',t0='90',tau=str(tau),R=str(r),axis=str(a),Nv=str(n),Nsrc=str(64**3*(3 if a==-1 else 1)),L_RE=str(z.real),L_IM=str(z.imag))
    return result


def write_references(path,configs=PILOT):
    rows=list(references(configs).values())
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def scan_fixture(cfg,piter):
    output=[line('MATRIX_BEGIN',cfg=cfg,piter=piter,saved_Nv=100,
        source_cutoffs='40,60,80,100',sink_cutoffs='40,60,80,100',ranks=64,
        grid='1x4x4x4',t0=90,tau_min=0,tau_max=32,radii='4,12,24,32',
        temporal_bc='open',hyp_steps=1,temporal_change=.1,preserved_change=0)]
    for s,k,t,r,a in sorted(expected_keys()):
        z=value(cfg,s,k,t,r,a)
        output.append(line('MATRIX_RT',cfg=cfg,piter=piter,scheme='fixed',t0=90,
            tau=t,R=r,axis=a,Nv_src=s,Nv_sink=k,Nsrc=64**3*(3 if a==-1 else 1),L_RE=z.real,L_IM=z.imag))
    for prefix,keys in [('MATRIX_DIRECT',direct_keys()),('MATRIX_LEGACY',legacy_keys())]:
        for s,k,t,r,a,x in sorted(keys):
            z=value(cfg,s,k,t,r,a);bound=1e-20+5e-12*abs(z)
            output.append(line(prefix,cfg=cfg,piter=piter,tau=t,R=r,axis=a,rank=0,
                site=x,Nv_src=s,Nv_sink=k,fast_re=z.real,fast_im=z.imag,
                direct_re=z.real,direct_im=z.imag,abs_error=0.,tolerance=bound,status='PASS'))
    for a in range(3):output.append(line('PERIODIC_ROUNDTRIP',axis=a,cells=64,status='PASS'))
    for t in range(90,122):
        output.append(line('SU3_TIME',cfg=cfg,piter=piter,t=t,count=64**3,
            unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,
            unit_above_1em2=0,det_above_1em2=0))
    output.append(line('SU3_WINDOW',cfg=cfg,piter=piter,t0=90,t1=122,count=32*64**3,
        unit_max=3e-16,unit_rms=1e-16,det_max=4e-16,det_rms=2e-16,
        unit_above_1em2=0,det_above_1em2=0))
    for s,k,t,r,a in sorted(block_keys()):
        z=value(cfg,s,k,t,r,a)
        for rank in range(64):
            local=z+complex((rank-31.5)*1e-12,(rank-31.5)*1e-14)
            output.append(line('MATRIX_BLOCK',cfg=cfg,piter=piter,t0=90,tau=t,R=r,axis=a,
                Nv_src=s,Nv_sink=k,rank=rank,px=rank//16,py=(rank//4)%4,pz=rank%4,
                Nsrc=4096,L_RE=local.real,L_IM=local.imag))
    output.append(line('MATRIX_END',cfg=cfg,piter=piter,records=len(expected_keys()),status='PASS'))
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


def replace_record(text,prefix,changes,predicate=lambda row:True):
    """Change exactly one parsed diagnostic without editing similar records."""
    lines=text.splitlines(True)
    for i,original in enumerate(lines):
        if not original.startswith(prefix+' '):continue
        row=dict(token.split('=',1) for token in original.split()[1:])
        if predicate(row):
            row.update(changes);lines[i]=line(prefix,**row);return ''.join(lines)
    raise AssertionError('Synthetic record not found: '+prefix)
