"""Conditional fixed-source pilot exports; no fits or extra statistical samples."""
import math,sys
from pilot_plan import CUTOFFS,TAUS,execution_scans


def quadratic(matrix,coefficients):
    terms=[v.conjugate()*matrix[k,s]*w for k,v in coefficients.items() for s,w in coefficients.items()]
    value=complex(math.fsum(t.real for t in terms),math.fsum(t.imag for t in terms))
    if not math.isfinite(value.real) or not math.isfinite(value.imag):raise ValueError('Nonfinite fixed-source projection')
    return value,64*sys.float_info.epsilon*math.fsum(abs(t) for t in terms)


def derive(values,reference,sources):
    configs=sorted({int(r['cfg']) for r in values});matrices={};specs={}
    for cfg in configs:
        for spec in execution_scans(cfg):specs[cfg,spec['t0'],spec['piter']]=spec
    for row in values:
        if row['axis']!='-1':continue
        cfg,t0,piter,r,t,s,k=(int(row[key]) for key in ('cfg','t0','piter','R','tau','Nv_src','Nv_sink'))
        spec=specs[cfg,t0,piter];group=(cfg,spec['label'],t0,piter,spec['purpose'])
        matrices.setdefault(group,{}).setdefault((r,t),{})[k,s]=complex(float(row['L_RE']),float(row['L_IM']))
    for (cfg,s,k,t,r),row in reference.items():
        if cfg not in configs:continue
        group=(cfg,'reused_t90',90,10,'published_reference')
        matrices.setdefault(group,{}).setdefault((r,t),{})[k,s]=complex(float(row['L_RE']),float(row['L_IM']))
    vectors=sources['source_candidates'];output=[];lookup={}
    for (cfg,label,t0,piter,purpose),matrix in sorted(matrices.items()):
        for vector in vectors:
            r=vector['R'];coeff={n:complex(re,im) for n,re,im in zip(vector['Nv'],vector['coefficients_re'],vector['coefficients_im'])}
            for tau in TAUS:
                value,bound=quadratic(matrix[r,tau],coeff)
                if tau==0 and abs(value.imag)>bound:raise ValueError('Equal-time fixed-source projection exceeds its roundoff bound')
                row=dict(configuration=cfg,label=label,purpose=purpose,t0=t0,piter=piter,R=r,basis=vector['basis'],anchor=vector['anchor'],
                    tau=tau,training_configurations=68,Q_RE=value.real,Q_IM=value.imag)
                output.append(row);lookup[cfg,label,r,vector['basis'],vector['anchor'],tau]=value
    contrasts=[]
    for row in output:
        if row['purpose']=='published_reference':continue
        cfg,label,r,b,a,t=(row[k] for k in ('configuration','label','R','basis','anchor','tau'))
        other='shift122_p10' if row['purpose']=='projection_control' else 'reused_t90'
        difference=lookup[cfg,label,r,b,a,t]-lookup[cfg,other,r,b,a,t]
        contrasts.append(dict(configuration=cfg,label=label,reference_label=other,purpose=row['purpose'],t0=row['t0'],piter=row['piter'],
            R=r,basis=b,anchor=a,tau=t,delta_Q_RE=difference.real,delta_Q_IM=difference.imag))
    return output,contrasts
