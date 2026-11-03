#!/usr/bin/env python3
"""Render the signed complex pilot diagnostics from the numerical CSV exports."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
RESULTS=HERE/'results'
RADII=[4,12,24,32]
COLORS=['#1d5a90','#db6719','#2a8767','#8b4e99']


def read(name):
    with (RESULTS/(name+'.csv')).open(newline='') as f:return list(csv.DictReader(f))


def select(rows,**where):return [r for r in rows if all(r[k]==str(v) for k,v in where.items())]


def save(fig,name):
    for extension in ('png','svg'):fig.savefig(RESULTS/(name+'.'+extension),dpi=170,bbox_inches='tight')
    plt.close(fig)


def energy_curve(rows):
    rows=[r for r in rows if int(r['tau'])<16]
    x=np.array([int(r['tau'])+.5 for r in rows]);y=np.array([float(r['effective_energy_re']) if r['energy_display_valid']=='1' else np.nan for r in rows])
    return x,y


def render():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'semibold','svg.hashsalt':'qx_matrix_pilot'})
    metric=read('metric_spectrum');diagnostics=read('matrix_diagnostics');modes=read('complex_modes');fixed=read('fixed_projections')
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    names=[('full4','40 / 60 / 80 / 100'),('spread3','40 / 60 / 100'),('wide2','40 / 100')]
    for i,(name,label) in enumerate(names):
        values=[float(select(metric,sample='mean',R=r,basis=name,mode=0)[0]['condition']) for r in RADII]
        axes[0,0].plot(RADII,values,'o-',color=COLORS[i],label=label)
    axes[0,0].set(xlabel=r'$R/a$',ylabel=r'$\kappa_2[D^{-1}C(0)D^{-1}]$',yscale='log',title='Equal-time basis conditioning');axes[0,0].legend(title='Mode cutoffs',fontsize=9)
    for j in range(4):
        values=[float(select(metric,sample='mean',R=r,basis='full4',mode=j)[0]['eigenvalue']) for r in RADII]
        axes[0,1].plot(RADII,values,'o-',color=COLORS[j],label='eigenvalue %d'%(j+1))
    axes[0,1].set(xlabel=r'$R/a$',ylabel='Normalized metric eigenvalue',yscale='log',title='Full four-source metric is positive definite');axes[0,1].legend(fontsize=8)
    for i,radius in enumerate(RADII):
        chosen=select(diagnostics,sample='mean',R=radius,basis='full4');x=[int(r['tau']) for r in chosen]
        axes[1,0].plot(x,[max(float(r['relative_antihermitian_norm']),1e-8) for r in chosen],color=COLORS[i],label='R/a=%d'%radius)
        axes[1,1].plot(x,[float(r['condition']) for r in chosen],color=COLORS[i],label='R/a=%d'%radius)
    axes[1,0].set(xlabel=r'$\tau/a$',ylabel=r'$\|A-A^\dagger\|_F/(2\|A\|_F)$',yscale='log',ylim=(1e-5,2),title='Finite-time non-Hermiticity (original matrices)')
    axes[1,1].set(xlabel=r'$\tau/a$',ylabel=r'$\kappa_2[A(\tau)]$',yscale='log',title='Finite-time singular-value condition number')
    for ax in axes[1]:ax.legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.18)
    fig.suptitle('qX mixed-source pilot: n170, n640, n650\nA = equal-time-diagonally-scaled matrix; three-configuration mean',fontsize=13)
    save(fig,'matrix_conditioning_pilot')

    fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
    for col,radius in enumerate(RADII):
        for j in range(4):
            chosen=select(modes,sample='mean',R=radius,basis='full4',t_ref=0,mode=j)
            x=[int(r['tau']) for r in chosen]
            axes[0,col].plot(x,[float(r['lambda_re']) for r in chosen],color=COLORS[j],label='tracked mode %d'%j)
            axes[1,col].plot(x,[float(r['lambda_im']) for r in chosen],color=COLORS[j])
        axes[0,col].set(title='R/a=%d'%radius,yscale='symlog',ylim=(-.02,1.1));axes[0,col].set_yscale('symlog',linthresh=1e-3)
        axes[1,col].set_yscale('symlog',linthresh=1e-4);axes[1,col].set(xlabel=r'$\tau/a$')
        for ax in axes[:,col]:ax.axhline(0,color='.5',lw=.6);ax.grid(alpha=.18)
    axes[0,0].set_ylabel(r'Re $\lambda(\tau,0)$; symlog');axes[1,0].set_ylabel(r'Im $\lambda(\tau,0)$; symlog');axes[0,0].legend(fontsize=8)
    fig.suptitle('Original complex generalized pencils: four cutoffs, reference separation 0\nAll signs and imaginary parts retained; mode labels follow vector overlaps',fontsize=13)
    save(fig,'complex_principal_correlators_pilot')

    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    models=[('Nv40','Nv=40',COLORS[0]),('Nv100','Nv=100','.45'),('fixed_full4','Fixed four-source projection',COLORS[1]),('fixed_wide2','Fixed 40/100 projection',COLORS[2])]
    for ax,radius in zip(axes.flat,RADII):
        for model,label,color in models:
            chosen=select(fixed,sample='mean',R=radius,model=model,evaluation='training');x,y=energy_curve(chosen)
            ax.plot(x,y,'.-',ms=3,color=color,label=label)
        ax.set(title='R/a=%d'%radius,xlabel=r'$\tau/a+1/2$',ylabel=r'Re $\log[q(\tau)/q(\tau+a)]$',xlim=(0,16));ax.grid(alpha=.18)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Fixed-vector effective-energy diagnostics; optimization at tau/a=4\nPositive real part and |phase| <= 0.05 at both endpoints; gaps fail this display gate',fontsize=12)
    save(fig,'fixed_source_effective_energies_pilot')

    fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
    for row,basis in enumerate(('full4','wide2')):
        for col,radius in enumerate(RADII):
            ax=axes[row,col]
            for color,tref in zip(COLORS,(0,1,2,4)):
                chosen=select(modes,sample='mean',R=radius,basis=basis,t_ref=tref,mode=0);x,y=energy_curve(chosen)
                ax.plot(x,y,'.-',ms=3,color=color,label=r'$\tau_{ref}/a=%d$'%tref)
            ax.set(title='%s; R/a=%d'%('Four sources' if basis=='full4' else '40/100 sources',radius),xlabel=r'$\tau/a+1/2$',xlim=(0,16));ax.grid(alpha=.18)
    for ax in axes[:,0]:ax.set_ylabel('Leading tracked-mode effective energy')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Reference-time and basis sensitivity: three-configuration mean\nComplex nonsymmetric pencils; gaps follow the same positive-real / phase display gate',fontsize=12)
    save(fig,'reference_time_sensitivity_pilot')

    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ax,radius in zip(axes.flat,RADII):
        for i,cfg in enumerate((170,640,650)):
            chosen=select(fixed,sample='drop%d'%cfg,R=radius,model='fixed_full4',evaluation='heldout');x,y=energy_curve(chosen)
            ax.plot(x,y,'.--',ms=3,color=COLORS[i],alpha=.85,label='Test n%d; train on other two'%cfg)
        chosen=select(fixed,sample='mean',R=radius,model='fixed_full4',evaluation='training');x,y=energy_curve(chosen)
        ax.plot(x,y,'k-',lw=1.6,label='Full three-configuration projection')
        ax.set(title='R/a=%d'%radius,xlabel=r'$\tau/a+1/2$',ylabel='Fixed-source effective energy',xlim=(0,16));ax.grid(alpha=.18)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Held-out-configuration checks of the fixed four-source projection\nSensitivity diagnostics; three test configurations; no confidence intervals',fontsize=12)
    save(fig,'fixed_source_holdout_pilot')


if __name__=='__main__':render()
