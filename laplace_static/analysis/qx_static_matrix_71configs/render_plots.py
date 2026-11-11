#!/usr/bin/env python3
"""Publication-file figures for the fully paired source and matrix diagnostics."""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from replay import HERE,RADII,CUTOFFS
from statistics_tools import energies,jk_se

P=HERE/'results'
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':10,'legend.fontsize':8,
    'figure.dpi':120,'savefig.dpi':180,'svg.hashsalt':'qx-matrix-full71',
    'axes.spines.top':False,'axes.spines.right':False})
# Combine colour with shape, fill and line pattern; nearly coincident curves
# must remain distinguishable in print as well as on a colour display.
COLORS={'Nv40':'#666666','Nv100':'#7B3294','pair':'#0072B2',
        'triple':'#AA4499','full':'#D55E00'}
MARKERS={'Nv40':'^','Nv100':'D','pair':'o','triple':'v','full':'s'}
LINES={'Nv40':':','Nv100':'-.','pair':'-','triple':':','full':'--'}


def style(label):
    return dict(color=COLORS[label],marker=MARKERS[label],linestyle=LINES[label],
                markerfacecolor=COLORS[label] if label=='full' else 'white',
                markeredgecolor=COLORS[label],markeredgewidth=.9)


def read(name):
    with (P/(name+'.csv')).open(newline='') as stream:return list(csv.DictReader(stream))


def select(rows,**where):
    return [r for r in rows if all(r[k]==str(v) for k,v in where.items())]


def val(rows,key):return np.array([float(r[key]) if r[key]!='' else np.nan for r in rows])


def save(fig,name,title,note):
    fig.suptitle(title,fontsize=15,y=.99)
    fig.text(.05,.014,note,fontsize=8,color='#444444',va='bottom')
    fig.tight_layout(rect=(.01,.06,.99,.94))
    for ext in ('png','svg'):fig.savefig(P/(name+'.'+ext),metadata={'Date':None} if ext=='svg' else None)
    plt.close(fig)


def base_axes(size=(12,8)):
    fig,axes=plt.subplots(2,2,figsize=size)
    for ax,r in zip(axes.flat,RADII):
        ax.set_title(r'$R/a={}$'.format(r));ax.grid(alpha=.16)
    return fig,list(axes.flat)


def metric():
    rows=read('metric_spectra');fig,axes=base_axes()
    for ax,r in zip(axes,RADII):
        for basis,offset in [('pair',-.1),('triple',0),('full',.1)]:
            rs=select(rows,R=r,basis=basis,sample='ensemble')
            x=val(rs,'mode')+1+offset;y=val(rs,'correlation_eigenvalue')
            ax.plot(x,y,ms=5,lw=1.1,label=basis+' basis',**style(basis))
            ranges=[]
            for d in rs:
                other=[float(a['correlation_eigenvalue']) for a in select(rows,R=r,basis=basis,mode=d['mode']) if a['sample']!='ensemble']
                ranges.append([min(other),max(other)])
            lo,hi=np.asarray(ranges).T
            ax.vlines(x,lo,hi,color=COLORS[basis],lw=4,alpha=.3)
        full=select(rows,R=r,basis='full',sample='ensemble')[0]
        ax.text(.97,.08,r'Full-basis $\kappa={:.1f}$'.format(float(full['correlation_condition'])),transform=ax.transAxes,ha='right')
        ax.set(yscale='log',xlabel='Equal-time correlation-metric mode',ylabel='Metric eigenvalue')
        ax.set_xticks([1,2,3,4])
    axes[0].legend(loc='upper left')
    save(fig,'matrix_conditioning_71configs','Equal-time metrics: all 71 configurations',
         'Markers: full ensemble. Thick bars: range across all 71 delete-one samples (not confidence intervals). No metric-mode deletion.')


def complexity():
    rows=read('matrix_defects');fig,axes=base_axes()
    for ax,r in zip(axes,RADII):
        rs=select(rows,R=r);x=val(rs,'tau')
        for field,label,series in [('relative_antihermitian_norm',r'$\|A-A^\dagger\|/(2\|A\|)$','pair'),
                ('relative_imaginary_norm',r'$\|\mathrm{Im}\,A\|/\|A\|$','full')]:
            y=val(rs,field);y[y<1e-14]=np.nan
            ax.plot(x,y,ms=3.5,lw=1.,label=label,**style(series))
        ax.set(yscale='log',xlabel=r'$\tau/a$',ylabel='Relative Frobenius norm',xlim=(0,32))
        if r==32:ax.text(.05,.84,'Imaginary part is at roundoff;\nreal matrix can still be nonsymmetric.',transform=ax.transAxes,fontsize=9)
    axes[0].legend()
    save(fig,'matrix_complexity_71configs','The full complex matrix is retained at nonzero time',
         'A = D(0)^(-1) C(t) D(0)^(-1). Central algebraic defects, not hypothesis tests. Values below 1e-14 are omitted as roundoff.')


def energy_plot(detail=False):
    base=read('diagonal_effective_energies');projected=read('projected_effective_energies')
    fig,axes=base_axes((13,8.5))
    for ax,r in zip(axes,RADII):
        candidates=[('Nv40',select(base,R=r,Nv=40),'jackknife_se',-.18),('Nv100',select(base,R=r,Nv=100),'jackknife_se',-.06),
                    ('pair',select(projected,R=r,basis='pair',anchor=4,procedure='cross_fit'),'jackknife_se',.06),
                    ('full',select(projected,R=r,basis='full',anchor=4,procedure='cross_fit'),'jackknife_se',.18)]
        if detail:candidates=[item for item in candidates if item[0] in ('pair','full')]
        for label,rs,err,shift in candidates:
            x=val(rs,'midpoint');y=val(rs,'aE_eff');se=val(rs,err)
            visible=(x>=6)&(x<=22) if detail else np.ones(len(x),bool)
            ax.errorbar(x[visible]+shift,y[visible],yerr=se[visible],ms=4,lw=.95,capsize=1.5,
                        label=label+' held out' if label in ('pair','full') else label,**style(label))
            if label in ('pair','full'):
                flagged=np.array([d['phase_within_005']=='False' or d['signal_at_least_3']=='False' for d in rs]) & visible
                ax.plot(x[flagged]+shift,y[flagged],'x',color='black',ms=6,mew=1.)
        ax.set(xlabel=r'Effective-energy midpoint $\tau/a+1/2$',ylabel=r'$aE_{\rm eff}$ (real-part ratio)',xlim=(5.7,22.3) if detail else (0,32))
    axes[0].legend(ncol=2)
    save(fig,'heldout_effective_energies_detail' if detail else 'fixed_source_effective_energies_71configs',
         'Fixed-source effective energies: held-out evaluation'+(' (times 6–21)' if detail else ''),
         'Training time 4. Only the two held-out source bases are shown in this detail view. Both reference diagonals and all times appear in the companion figure.\n'
         'Errors include complete delete-one source retraining. Small horizontal offsets. Black x: phase >0.05 rad or endpoint signal <3.' if detail else
         'Training time 4. Errors include complete delete-one source retraining. Small horizontal offsets. Black x: phase >0.05 rad or endpoint signal <3.\n'
         'All measured times shown. Gaps are undefined real-mean log ratios; no positive-configuration filtering. These curves do not establish ground-state energies.')


def signal_phase():
    base=read('diagonal_correlators');pr=read('projected_correlators')
    fig,axes=plt.subplots(2,4,figsize=(15,8))
    for j,r in enumerate(RADII):
        for label in ('Nv40','Nv100','pair','full'):
            if label.startswith('Nv'):
                rs=select(base,R=r,Nv=int(label[2:]));m=val(rs,'mean_re')+1j*val(rs,'mean_im');se=val(rs,'independent_se')
            else:
                rs=select(pr,R=r,basis=label,anchor=4,procedure='cross_fit');m=val(rs,'mean_re')+1j*val(rs,'mean_im');se=val(rs,'jackknife_se_re')
            x=val(rs,'tau');signal=m.real/se
            axes[0,j].plot(x,signal,ms=3.3,lw=1.,label=label,**style(label))
            axes[1,j].plot(x,np.angle(m),ms=3.3,lw=1.,**style(label))
        axes[0,j].set_title(r'$R/a={}$'.format(r));axes[0,j].set_yscale('symlog',linthresh=3)
        axes[0,j].axhline(3,ls=':',color='black',lw=.8)
        axes[1,j].axhline(.05,ls=':',color='black',lw=.8);axes[1,j].axhline(-.05,ls=':',color='black',lw=.8)
        for i in (0,1):axes[i,j].set_xlim(0,32);axes[i,j].grid(alpha=.15);axes[i,j].set_xlabel(r'$\tau/a$')
    axes[0,0].set_ylabel('Real mean / uncertainty');axes[1,0].set_ylabel('Correlator phase [rad]')
    axes[0,0].legend(ncol=2)
    save(fig,'heldout_signal_and_phase_71configs','Signal improvement and complex-phase control',
         'Training time 4. Diagonals: independent configuration SE. Held-out sources: complete-retraining jackknife SE. Dotted lines are review flags, not acceptance criteria.')


def paired_sensitivity():
    rows=read('paired_energy_contrasts');fig,axes=plt.subplots(2,4,figsize=(15,8))
    for j,r in enumerate(RADII):
        for i,(contrast,items) in enumerate([('minus_anchor4',[(6,'pair'),(8,'full')]),
                ('minus_pair_basis',[('triple','triple'),('full','full')])]):
            ax=axes[i,j]
            for item,series in items:
                rs=select(rows,R=r,basis='full',anchor=item,procedure='cross_fit',contrast=contrast) if i==0 else select(rows,R=r,basis=item,anchor=4,procedure='cross_fit',contrast=contrast)
                rs=[d for d in rs if 6<=int(d['tau'])<=21]
                label='train '+str(item)+' minus 4' if i==0 else str(item)+' minus pair'
                ax.errorbar(val(rs,'midpoint'),val(rs,'difference'),yerr=val(rs,'paired_jackknife_se'),
                            ms=3.5,lw=.95,label=label,**style(series))
            ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.15);ax.set_xlabel(r'$\tau/a+1/2$')
            if i==0:ax.set_title(r'$R/a={}$'.format(r))
    axes[0,0].set_ylabel(r'Paired $\Delta(aE_{\rm eff})$');axes[1,0].set_ylabel(r'Paired $\Delta(aE_{\rm eff})$')
    axes[0,0].legend();axes[1,0].legend()
    save(fig,'paired_source_sensitivity_71configs','Training time and basis sensitivity: paired differences',
         'Held-out procedure throughout. All 71 deletions retrain each source; covariance between candidates is retained. Pointwise error bars; no multiplicity-adjusted significance claim.')


def holdout():
    rows=read('paired_energy_contrasts');fig,axes=base_axes()
    for ax,r in zip(axes,RADII):
        for basis,shift in [('pair',-.06),('full',.06)]:
            rs=select(rows,R=r,basis=basis,anchor=4,procedure='cross_fit',contrast='minus_same_sample')
            ax.errorbar(val(rs,'midpoint')+shift,val(rs,'difference'),yerr=val(rs,'paired_jackknife_se'),
                        ms=4,lw=.95,label=basis,**style(basis))
        ax.axhline(0,color='black',lw=.8);ax.set(xlabel=r'$\tau/a+1/2$',ylabel=r'Held-out minus same-sample $aE_{\rm eff}$',xlim=(0,32))
    axes[0].legend()
    save(fig,'holdout_comparison_71configs','Does training on the same configurations change the result?',
         'Training time 4. Paired complete-retraining jackknife; the two procedures share configurations and must not be compared with independent errors.')


def window_plot():
    rows=read('constant_window_diagnostics');fig,axes=base_axes()
    starts=[6,10,14,18,22]
    for ax,r in zip(axes,RADII):
        for basis,shift in [('pair',-.12),('full',.12)]:
            rs=select(rows,R=r,basis=basis,anchor=4,procedure='cross_fit')
            x=val(rs,'start_tau');y=val(rs,'aE');se=val(rs,'jackknife_se')
            ax.errorbar(x+shift,y,yerr=se,ms=5,capsize=2,lw=1.1,label=basis+' held out',**style(basis))
            flagged=np.array([d['phase_and_signal_screen']=='False' for d in rs])
            ax.plot(x[flagged]+shift,y[flagged],'x',ms=8,color='black')
            for a,xx,yy in zip(rs,x+shift,y):
                if a['gaussian_p_approx'] and float(a['gaussian_p_approx'])<.05:ax.plot(xx,yy,'s',mfc='none',mec='black',ms=8)
        ax.set(xlabel='First effective-energy time in six-point window',ylabel='Conditional fitted constant aE')
        ax.set_xticks(starts);ax.set_xticklabels([str(s)+'–'+str(s+5) for s in starts])
    axes[0].legend()
    save(fig,'constant_window_stability_71configs','Window dependence remains a separate test from local flatness',
         'Training time 4. Black x: phase/signal screen fails. Open square: approximate local p<0.05. Estimated covariance and many tested windows limit p-value interpretation.\n'
         'At R=32 a locally flat late window still differs from the early window. These are candidate diagnostics, not a selected static potential.')


def generalized():
    rows=read('complex_pencil_diagnostics');fig,axes=base_axes()
    for ax,r in zip(axes,RADII):
        for ref,series in [(0,'pair'),(2,'triple'),(4,'full')]:
            rs=select(rows,R=r,reference_tau=ref)
            ax.scatter(val(rs,'eigenvalue_re'),val(rs,'eigenvalue_im'),s=15,marker=MARKERS[series],
                       alpha=.8,label='reference '+str(ref),edgecolors=COLORS[series],
                       facecolors=COLORS[series] if series=='full' else 'none',linewidths=.75)
        ax.axhline(0,color='black',lw=.5);ax.set(xlabel=r'Re $\lambda$',ylabel=r'Im $\lambda$')
    axes[0].legend()
    save(fig,'complex_pencil_spectra_71configs','General complex pencils: algebraic diagnostics only',
         'All four roots for each measured time beyond the reference. No Hermitian replacement, physical mode ordering or energy identification. Backward residuals <1e-10.')


if __name__=='__main__':
    metric();complexity();energy_plot(False);energy_plot(True);signal_phase();paired_sensitivity();holdout();window_plot();generalized()
    print('Rendered nine figure pairs (PNG/SVG).')
