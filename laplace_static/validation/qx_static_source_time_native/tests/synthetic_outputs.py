"""Analytic parser/lifecycle fixtures only; never physical gauge measurements."""
import csv,math
from parents import legacy,load,MATRIX
from pilot_plan import PILOT,CUTOFFS,RADII,TAUS,AXES,expected_keys,direct_keys,legacy_keys,block_keys
from parents import design
old=load(MATRIX/'tests/synthetic_outputs.py','qx_parent_synthetic_matrix_outputs')
line=old.line
replace_record=old.replace_record


def value(cfg,source,sink,tau,radius,axis,t0=90):
    return old.value(cfg,source,sink,tau,radius,axis)*(1+.001*(t0-90))*math.exp(-.0004*(t0-90)*tau)


def write_references(path):
    output=[]
    for cfg in PILOT:
        for s in CUTOFFS:
            for k in CUTOFFS:
                for t in TAUS:
                    for r in RADII:
                        z=value(cfg,s,k,t,r,-1)
                        output.append(dict(configuration=cfg,t0=90,t1=90+t,tau=t,R=r,axis=-1,piter=10,Nv_src=s,Nv_sink=k,Nsrc=3*64**3,L_RE=z.real,L_IM=z.imag))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(output[0]),lineterminator='\n');writer.writeheader();writer.writerows(output)


def audit_fixture(cfg):
    times=design.audit_times(cfg);output=[line('QUALITY_START',cfg=cfg,grid='1x4x4x4',Nv=100,APE='20,0.5',FUNNELED=1,smear_seconds=1)]
    for t in times:
        for i in range(100):
            lam=.01+.001*i
            output.append(line('MODE',cfg=cfg,t=t,i=i,**{'lambda':lam,'norm':1,'residual_l2':lam*5e-12,
                'residual_relative':5e-12,'rayleigh_re':lam,'rayleigh_im':0}))
        output.append(line('SLICE',cfg=cfg,t=t,Nv=100,max_relative_residual=5e-12,max_norm_defect=0,gram_diag_defect=1e-15,
            gram_offdiag=1e-15,max_rayleigh_im=0,lambda_min=.01,lambda_max=.109,status='PASS'))
    return ''.join(output)+line('QUALITY_END',cfg=cfg,slices=len(times),status='PASS')


def scan_fixture(cfg,spec):
    t0,piter=spec['t0'],spec['piter']
    output=[line('MATRIX_BEGIN',protocol='source_time_v1',cfg=cfg,piter=piter,saved_Nv=100,source_cutoffs='40,60,80,100',
        sink_cutoffs='40,60,80,100',ranks=64,grid='1x4x4x4',t0=t0,t1_max=t0+32,tau_min=0,tau_max=32,radii='4,12,24,32',
        temporal_bc='open',hyp_steps=1,temporal_change=.1,preserved_change=0)]
    output.append(line('MODE_LOAD',cfg=cfg,piter=piter,t0=t0,loaded_t=t0,role='source',status='PASS'))
    for t in range(1,33):output.append(line('MODE_LOAD',cfg=cfg,piter=piter,t0=t0,loaded_t=t0+t,role='sink',status='PASS'))
    for s,k,t,r,a in sorted(expected_keys()):
        z=value(cfg,s,k,t,r,a,t0)
        output.append(line('MATRIX_RT',cfg=cfg,piter=piter,scheme='fixed',t0=t0,t1=t0+t,tau=t,R=r,axis=a,Nv_src=s,Nv_sink=k,
            Nsrc=64**3*(3 if a==-1 else 1),L_RE=z.real,L_IM=z.imag))
    for prefix,keys in [('MATRIX_DIRECT',direct_keys()),('MATRIX_LEGACY',legacy_keys())]:
        for s,k,t,r,a,x in sorted(keys):
            z=value(cfg,s,k,t,r,a,t0)
            output.append(line(prefix,cfg=cfg,piter=piter,t0=t0,t1=t0+t,tau=t,R=r,axis=a,rank=0,site=x,Nv_src=s,Nv_sink=k,
                fast_re=z.real,fast_im=z.imag,direct_re=z.real,direct_im=z.imag,abs_error=0,tolerance=1e-20+5e-12*abs(z),status='PASS'))
    for a in range(3):output.append(line('PERIODIC_ROUNDTRIP',axis=a,cells=64,status='PASS'))
    for t in range(t0,t0+32):
        output.append(line('SU3_TIME',cfg=cfg,piter=piter,t=t,count=64**3,unit_max=3e-16,unit_rms=1e-16,
            det_max=4e-16,det_rms=2e-16,unit_above_1em2=0,det_above_1em2=0))
    output.append(line('SU3_WINDOW',cfg=cfg,piter=piter,t0=t0,t1=t0+32,count=32*64**3,unit_max=3e-16,unit_rms=1e-16,
        det_max=4e-16,det_rms=2e-16,unit_above_1em2=0,det_above_1em2=0))
    for s,k,t,r,a in sorted(block_keys()):
        z=value(cfg,s,k,t,r,a,t0)
        for rank in range(64):
            local=z+complex((rank-31.5)*1e-12,(rank-31.5)*1e-14)
            output.append(line('MATRIX_BLOCK',cfg=cfg,piter=piter,t0=t0,t1=t0+t,tau=t,R=r,axis=a,Nv_src=s,Nv_sink=k,rank=rank,
                px=rank//16,py=rank//4%4,pz=rank%4,Nsrc=4096,L_RE=local.real,L_IM=local.imag))
    return ''.join(output)+line('MATRIX_END',protocol='source_time_v1',cfg=cfg,piter=piter,t0=t0,t1_max=t0+32,records=len(expected_keys()),status='PASS')
