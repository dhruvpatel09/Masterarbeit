#!/usr/bin/env python3
"""Render complete measured-time coverage and explicitly flagged fit diagnostics."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from data import HERE, LONG_RADII, CUTOFFS, FOCUS

COLORS = {40:'#0072B2',60:'#D55E00',80:'#009E73',100:'#9B51A5'}
MARKERS = {40:'o',60:'s',80:'^',100:'D'}
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,
    'legend.fontsize':9,'svg.fonttype':'none','svg.hashsalt':'qx-radius-71-2026-11-08','savefig.dpi':150,
    'axes.spines.top':False,'axes.spines.right':False})


def table(folder,name):
    with (folder/(name+'.csv')).open(newline='') as stream:
        return list(csv.DictReader(stream))


def num(row,key):
    return float(row[key]) if row[key] else np.nan


def save(fig,folder,name):
    if fig.legends:
        # Reserve physical space for title, subtitle, figure legend and axes titles.
        height=fig.get_size_inches()[1]
        fig.subplots_adjust(top=min(fig.subplotpars.top,1-1.4/height))
    for ext in ('png','svg'):
        fig.savefig(folder/(name+'.'+ext),bbox_inches='tight',metadata={'Date':None} if ext=='svg' else {})
    plt.close(fig)


def decorate(ax,xlabel=None,ylabel=None):
    ax.grid(alpha=.17)
    if xlabel:ax.set_xlabel(xlabel)
    if ylabel:ax.set_ylabel(ylabel)


def heading(fig,title,subtitle,footer):
    height=fig.get_size_inches()[1]
    fig.suptitle(title,fontsize=15,y=1-.13/height,va='top')
    fig.text(.5,1-.50/height,subtitle,ha='center',va='top',fontsize=10)
    fig.text(.5,.015,footer,ha='center',va='bottom',fontsize=9)


def cutoff_legend(fig):
    height=fig.get_size_inches()[1]
    fig.legend(handles=legend_handles(),loc='upper center',
               bbox_to_anchor=(.5,1-.82/height),ncol=4,frameon=False)


def legend_handles():
    return [Line2D([],[],color=COLORS[n],marker=MARKERS[n],lw=.8,label=r'$N_v={}$'.format(n)) for n in CUTOFFS]


def energy_points(ax,records,n,zoom=None):
    selected=sorted([r for r in records if int(r['Nv'])==n],key=lambda r:int(r['tau']))
    x=np.array([float(r['tau_midpoint']) for r in selected])
    y=np.array([num(r,'aE_eff') for r in selected]);s=np.array([num(r,'naive_se') for r in selected])
    strong=np.array([r['valid_point']=='True' and r['low_signal']=='False' for r in selected])
    defined=np.isfinite(y)
    # A NaN breaks lines. Undefined points are placed on a validity strip, not at zero energy.
    ax.plot(x,y,color=COLORS[n],lw=.7,alpha=.65)
    for mask,face in ((defined&strong,COLORS[n]),(defined&~strong,'white')):
        if mask.any():
            ax.errorbar(x[mask],y[mask],yerr=s[mask],fmt=MARKERS[n],ms=3.5,
                mfc=face,mec=COLORS[n],color=COLORS[n],elinewidth=.6,capsize=1.3)
    if (~defined).any():
        ax.plot(x[~defined],np.full(np.sum(~defined),.025+.022*CUTOFFS.index(n)),
                transform=ax.get_xaxis_transform(),marker='x',ls='none',ms=4,color=COLORS[n])
    if zoom is not None:
        lo,hi=zoom;ax.set_ylim(lo,hi)
        for mask,edge,marker in ((defined&(y<lo),lo,'v'),(defined&(y>hi),hi,'^')):
            if mask.any():ax.plot(x[mask],np.full(mask.sum(),edge),marker=marker,ls='none',color=COLORS[n],ms=4,clip_on=False)
    if np.any(y-s<0) and np.any(y+s>0):
        ax.axhline(0,color='0.65',lw=.6)


def zoom_range(records):
    vals=[]
    for row in records:
        if row['valid_point']=='True' and 6<=int(row['tau'])<=15 and row['low_signal']=='False':
            vals += [num(row,'aE_eff')-2*num(row,'naive_se'),num(row,'aE_eff')+2*num(row,'naive_se')]
    if not vals:return (-.1,.6)
    low,high=min(vals),max(vals);pad=max((high-low)*.25,.002)
    return low-pad,high+pad


def energy_atlases(folder,energies):
    # Every integer radius, using the common dense short-time scan.
    for page in range(4):
        fig,axes=plt.subplots(4,2,figsize=(12,13))
        rs=list(range(page*8+1,page*8+9))
        for ax,r in zip(axes.flat,rs):
            rows=[x for x in energies if int(x['R'])==r and int(x['Nv'])==100 and int(x['tau'])<=15]
            energy_points(ax,rows,100);ax.set_xlim(-.2,16)
            ax.set_title(r'$R/a={}$'.format(r));decorate(ax,r'$(\tau+a/2)/a$',r'$aE_{\rm eff}$')
        heading(fig,'qX: effective energies at every spatial distance',
            r'71 configurations; $p_{\rm iter}=10$; $N_v=100$; source $t_0/a=90$; full short-time range',
            'Points are log-ratios of ensemble means; bars are independent-configuration errors. No plateau selected.')
        fig.subplots_adjust(top=.92,bottom=.075,hspace=.44,wspace=.24)
        save(fig,folder,'effective_energy_dense_R{:02d}_R{:02d}'.format(rs[0],rs[-1]))
    for page in range(3):
        fig,axes=plt.subplots(4,2,figsize=(14,14))
        rs=LONG_RADII[4*page:4*page+4]
        for i,r in enumerate(rs):
            rows=[x for x in energies if int(x['R'])==r]
            limits=zoom_range(rows)
            for col in (0,1):
                ax=axes[i,col]
                for n in CUTOFFS:energy_points(ax,rows,n,limits if col else None)
                ax.set_xlim(-.5,32)
                ax.set_title(r'$R/a={}$: {}'.format(r,'full error-bar range' if col==0 else 'detail view'))
                decorate(ax,r'$(\tau+a/2)/a$',r'$aE_{\rm eff}$')
        heading(fig,'qX: complete extended-time effective energies',
            r'71 configurations; $p_{\rm iter}=10$; $\tau/a=0,\ldots,32$ for correlators; ratios end at midpoint 31.5',
            'Hollow: endpoint mean/SE < 3. Crosses on lower strip: undefined log-ratio. Detail clips bars; edge triangles mark off-scale means.')
        cutoff_legend(fig)
        fig.subplots_adjust(top=.905,bottom=.075,hspace=.44,wspace=.24)
        save(fig,folder,'effective_energy_extended_page{}'.format(page+1))


def focus_signals(folder,means):
    fig,axes=plt.subplots(5,2,figsize=(13,16))
    for i,r in enumerate(FOCUS):
        for n in CUTOFFS:
            rows=sorted([x for x in means if int(x['R'])==r and int(x['Nv'])==n],key=lambda x:int(x['tau']))
            t=[int(x['tau']) for x in rows]
            axes[i,0].errorbar(t,[num(x,'normalized_L') for x in rows],
                yerr=[num(x,'normalized_naive_se') for x in rows],color=COLORS[n],marker=MARKERS[n],ms=2.5,lw=.7,capsize=1)
            axes[i,1].plot(t,[num(x,'mean_over_se') for x in rows],color=COLORS[n],marker=MARKERS[n],ms=2.7,lw=.8)
        axes[i,0].set_yscale('symlog',linthresh=1e-6);axes[i,0].axhline(0,color='0.5',lw=.6)
        axes[i,1].set_yscale('symlog',linthresh=3)
        for v in (0,3):axes[i,1].axhline(v,color='0.5',ls='--' if v else '-',lw=.7)
        for col in (0,1):
            axes[i,col].set_title(r'$R/a={}$'.format(r));axes[i,col].set_xlim(-.3,32.3)
            decorate(axes[i,col],r'$\tau/a$',r'$\overline{C}(\tau)/\overline{C}(0)$' if col==0 else r'$\overline{C}(\tau)/{\rm SE}[\overline{C}(\tau)]$')
    heading(fig,'qX: signed signal around R/a = 24',
        'All five measured neighboring distances; normalization and errors preserve configuration pairing',
        'Symmetric log axes retain negative means. Signed mean/SE is a signal diagnostic, not an extra sample count or a plateau test.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.91,bottom=.065,hspace=.44,wspace=.25)
    save(fig,folder,'signed_signal_R20_R28')
    fig,axes=plt.subplots(3,1,figsize=(11,12))
    for n in CUTOFFS:
        rows=sorted([x for x in means if int(x['R'])==24 and int(x['Nv'])==n],key=lambda x:int(x['tau']))
        t=[int(x['tau']) for x in rows]
        axes[0].errorbar(t,[num(x,'mean_L') for x in rows],yerr=[num(x,'naive_se') for x in rows],color=COLORS[n],marker=MARKERS[n],ms=3,lw=.8,capsize=1.2)
        axes[1].errorbar(t,[num(x,'normalized_L') for x in rows],yerr=[num(x,'normalized_naive_se') for x in rows],color=COLORS[n],marker=MARKERS[n],ms=3,lw=.8,capsize=1.2)
        axes[2].plot(t,[num(x,'mean_over_se') for x in rows],color=COLORS[n],marker=MARKERS[n],ms=3,lw=.8)
    axes[0].set_yscale('symlog',linthresh=1e-14)
    axes[1].set_yscale('symlog',linthresh=1e-6)
    axes[2].set_yscale('symlog',linthresh=3);axes[2].axhline(3,color='0.5',ls='--',lw=.8)
    for ax,label in zip(axes,[r'$\overline{C}(\tau)$',r'$\overline{C}(\tau)/\overline{C}(0)$',r'Signed mean / SE']):
        ax.axhline(0,color='0.6',lw=.6);ax.set_xlim(-.3,32.3);decorate(ax,r'$\tau/a$',label)
    heading(fig,'qX: R/a = 24 signal audit',
        '71 independent configurations; all signed measurements retained; no absolute-value logarithms',
        'An apparently flat late effective energy must be checked against these endpoint signals and the correlated window scan.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.9,bottom=.08,hspace=.34)
    save(fig,folder,'R24_signal_audit')


def fit_points(ax,rows,n,xkey):
    chosen=sorted([r for r in rows if int(r['Nv'])==n],key=lambda r:int(r[xkey]))
    for row in chosen:
        if row['status']=='outside_measured_time':
            continue
        x=int(row[xkey])+.11*(CUTOFFS.index(n)-1.5)
        if row['status']=='defined':
            screen=row['diagnostic_screen']=='True'
            ax.errorbar(x,num(row,'estimate'),yerr=num(row,'naive_se'),fmt=MARKERS[n],
                color=COLORS[n],mfc=COLORS[n] if screen else 'white',ms=4,capsize=1.5,lw=.7)
        else:
            ax.plot(x,.025+.022*CUTOFFS.index(n),transform=ax.get_xaxis_transform(),marker='x',ms=4,color=COLORS[n])


def window_plots(folder,fits,contrasts):
    for page in range(2):
        fig,axes=plt.subplots(3,2,figsize=(13,12))
        rs=LONG_RADII[6*page:6*page+6]
        for ax,r in zip(axes.flat,rs):
            rows=[x for x in fits if int(x['R'])==r and int(x['tau_max'])==int(x['tau_min'])+5 and int(x['tau_min'])%2==0]
            for n in CUTOFFS:fit_points(ax,rows,n,'tau_min')
            ax.set_title(r'$R/a={}$'.format(r));decorate(ax,r'Window start $\tau_{\min}/a$ (six ratios)',r'Constant estimate $aE$');ax.set_xlim(3,27)
        heading(fig,'qX: distance-by-distance window stability',
            'Full time covariance; constant fits to six consecutive log-ratios; paired delete-one errors also exported',
            'Filled: passes exploratory signal/p/condition screen. Hollow: fails at least one screen. Cross: undefined. No window selected.')
        cutoff_legend(fig)
        fig.subplots_adjust(top=.9,bottom=.08,hspace=.39,wspace=.25)
        save(fig,folder,'constant_window_scan_page{}'.format(page+1))
    # Vary the endpoint at a fixed start, rather than only sliding fixed-width windows.
    fig,axes=plt.subplots(3,2,figsize=(13,12))
    for ax,r in zip(axes.flat,[12,20,23,24,25,28]):
        rows=[x for x in fits if int(x['R'])==r and int(x['tau_min'])==8 and int(x['tau_max']) in (13,15,19,23,27,31)]
        for n in CUTOFFS:fit_points(ax,rows,n,'tau_max')
        ax.set_title(r'$R/a={}$'.format(r));decorate(ax,r'Window end $\tau_{\max}/a$, start fixed at 8',r'Constant estimate $aE$')
    heading(fig,'qX: endpoint sensitivity of candidate plateaus',
        'Extending a fit through weak late signal need not improve its physical interpretation',
        'Filled/hollow/cross meanings as in window scans. These estimates are correlated; use paired differences for comparisons.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.9,bottom=.08,hspace=.4,wspace=.24)
    save(fig,folder,'window_endpoint_sensitivity')
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    pairs=[(6,11,14,19),(14,19,22,27),(8,19,8,31),(8,23,8,31)]
    for ax,pair in zip(axes.flat,pairs):
        for n in CUTOFFS:
            rows=[r for r in contrasts if int(r['Nv'])==n and tuple(int(r[k]) for k in ('reference_tau_min','reference_tau_max','comparison_tau_min','comparison_tau_max'))==pair]
            for row in rows:
                if row['status']!='defined':continue
                x=int(row['R'])+.13*(CUTOFFS.index(n)-1.5)
                ax.errorbar(x,num(row,'difference'),yerr=num(row,'paired_naive_se'),fmt=MARKERS[n],color=COLORS[n],mfc='white' if row['low_signal']=='True' else COLORS[n],ms=4,capsize=1,lw=.7)
        ax.axhline(0,color='0.4',lw=.8);decorate(ax,r'$R/a$',r'Later/extended minus reference $aE$')
        ax.set_title('[{},{}] minus [{},{}]'.format(pair[2],pair[3],pair[0],pair[1]))
    heading(fig,'qX: paired window changes',
        'Errors include covariance of both fit estimates on the same 71 configurations',
        'Hollow: at least one fit uses endpoint mean/SE < 3. Undefined contrasts are retained in CSV, not replaced by zero.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.88,bottom=.1,hspace=.38,wspace=.25)
    save(fig,folder,'paired_window_changes')


def spatial_plots(folder,spatial,means):
    fig,axes=plt.subplots(2,2,figsize=(13,10))
    for ax,(lo,hi) in zip(axes.flat,[(6,11),(14,19),(22,27),(8,31)]):
        selected=[r for r in spatial if int(r['tau_min'])==lo and int(r['tau_max'])==hi]
        for n in CUTOFFS:fit_points(ax,selected,n,'R')
        ax.set_title(r'Ratio window $[{},{}]$'.format(lo,hi));decorate(ax,r'$R/a$',r'Candidate constant $aE$');ax.set_xlim(.2,33)
    heading(fig,'qX: candidate energies across spatial distance',
        'Dense integer R coverage exists only in the short scan; later windows use 12 measured distances',
        'No symbol outside measured time. Cross: undefined measured window. Filled/hollow: fit screens, not a selected potential.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.87,bottom=.1,hspace=.4,wspace=.25)
    save(fig,folder,'spatial_window_comparison')
    fig,axes=plt.subplots(1,2,figsize=(13,5.5))
    for n in CUTOFFS:
        rows=sorted([r for r in means if int(r['Nv'])==n and int(r['tau'])==0],key=lambda r:int(r['R']))
        x=np.array([int(r['R']) for r in rows]);y=np.array([num(r,'mean_L') for r in rows]);s=np.array([num(r,'naive_se') for r in rows])
        # Keep the native equal-time normalization so the source-shape differences remain explicit.
        axes[0].errorbar(x,y,yerr=s,color=COLORS[n],fmt=MARKERS[n],ms=4,lw=.7,capsize=1.5)
        axes[1].errorbar(x[x>=20],y[x>=20],yerr=s[x>=20],color=COLORS[n],fmt=MARKERS[n],ms=4,lw=.7,capsize=1.5)
    for ax in axes:
        ax.set_yscale('log');decorate(ax,r'$R/a$',r'Equal-time correlator $\overline{C}(R,0)$')
    axes[0].set_title('All measured distances');axes[1].set_title('Large-distance detail')
    heading(fig,'qX: equal-time source shape near the spatial half-box',
        r'Periodic spatial extent $L/a=64$; $R/a=32$ is half-space, not half-time',
        'A source-shape minimum or symmetry-related flattening is not evidence for string breaking or a ground-state potential plateau.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.8,bottom=.2,wspace=.26)
    save(fig,folder,'half_box_source_shape')


def approach_plots(folder,energies,approach):
    fig,axes=plt.subplots(3,2,figsize=(13,12))
    for ax,r in zip(axes.flat,[4,12,20,23,24,25]):
        records=[x for x in energies if int(x['R'])==r]
        limits=zoom_range(records)
        for n in CUTOFFS:
            energy_points(ax,records,n,limits)
            row=next(x for x in approach if int(x['R'])==r and int(x['Nv'])==n and int(x['tau_min'])==8)
            if row['status']=='defined':
                t=np.linspace(8,31,160)
                value=num(row,'A')+num(row,'B')*np.exp(-num(row,'c')*(t-8))
                ax.plot(t+.5,value,color=COLORS[n],lw=1.4,ls='-' if row['diagnostic_screen']=='True' else '--')
        ax.set_title(r'$R/a={}$'.format(r));ax.set_xlim(3.5,32);decorate(ax,r'$(\tau+a/2)/a$',r'$aE_{\rm eff}$')
    heading(fig,'qX: leading exponential approach checks',
        r'$A+B\,e^{-c(\tau-\tau_{\min})/a}$ fitted to ratios 8,...,31 with their full covariance',
        'Dashed models fail at least one diagnostic. A is not automatically E0; errors are conditional, and the model is not an exact two-state formula.')
    cutoff_legend(fig)
    fig.subplots_adjust(top=.9,bottom=.08,hspace=.4,wspace=.25)
    save(fig,folder,'approach_model_checks')
    fig,axes=plt.subplots(4,2,figsize=(14,14))
    for i,n in enumerate(CUTOFFS):
        rows=[r for r in approach if int(r['Nv'])==n and int(r['tau_max'])==31]
        for col in (0,1):
            ax=axes[i,col]
            for start,offset,color in [(4,-.18,'#0072B2'),(8,0,'#D55E00'),(12,.18,'#009E73')]:
                for row in rows:
                    if int(row['tau_min'])!=start or row['status']!='defined':continue
                    a,se=num(row,'A'),num(row,'A_conditional_se')
                    if np.isfinite(se):
                        ax.errorbar(int(row['R'])+offset,a,yerr=se,fmt='o',mfc=color if row['diagnostic_screen']=='True' else 'white',color=color,ms=4,capsize=1,lw=.7)
                        if col and (a<-.05 or a>.65):
                            ax.plot(int(row['R'])+offset,-.05 if a<-.05 else .65,marker='v' if a<-.05 else '^',color=color,ms=4,clip_on=False)
            if col:ax.set_ylim(-.05,.65)
            ax.set_title(r'$N_v={}$: {}'.format(n,'full conditional-error range' if not col else 'detail view'))
            decorate(ax,r'$R/a$',r'Approach asymptote $A$ (conditional error)')
    heading(fig,'qX: approach-asymptote start-time sensitivity',
        'Windows end at ratio 31; start = 4 (blue), 8 (orange), 12 (green)',
        'Hollow: unresolved/weak/model flag. Detail clips bars; edge triangles mark off-scale means. Missing curvature/undefined cases remain in CSV.')
    fig.subplots_adjust(top=.9,bottom=.08,hspace=.42,wspace=.27)
    save(fig,folder,'approach_asymptote_sensitivity')


def axis_plot(folder,axisrows):
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    labels={(0,1):'axis 0 minus 1',(0,2):'axis 0 minus 2',(1,2):'axis 1 minus 2'}
    palette=['#0072B2','#D55E00','#009E73']
    for ax,n in zip(axes.flat,CUTOFFS):
        for color,pair in zip(palette,labels):
            rows=sorted([r for r in axisrows if int(r['R'])==24 and int(r['Nv'])==n and (int(r['axis']),int(r['reference_axis']))==pair],key=lambda r:int(r['tau']))
            x=[int(r['tau'])+.5 for r in rows];y=[num(r,'difference') for r in rows];e=[num(r,'paired_naive_se') for r in rows]
            ax.errorbar(x,y,yerr=e,color=color,fmt='o',ms=3,lw=.6,capsize=1,label=labels[pair])
        ax.axhline(0,color='0.5',lw=.7);ax.set_title(r'$R/a=24$, $N_v={}$'.format(n));decorate(ax,r'$(\tau+a/2)/a$',r'Axis energy difference')
    heading(fig,'qX: paired spatial-axis differences at R/a = 24',
        'Axes are measurements on the same configurations; their covariance is retained',
        'Full error-bar range shown; undefined axis ratios are gaps. Pointwise differences are exploratory, not a global rotational-symmetry test.')
    axes[0,0].legend(frameon=False,fontsize=8)
    fig.subplots_adjust(top=.88,bottom=.1,hspace=.4,wspace=.25)
    save(fig,folder,'R24_paired_axis_checks')


def coverage_plot(folder,coverage):
    fig,ax=plt.subplots(figsize=(12,4.7))
    for i,n in enumerate(CUTOFFS):
        rows=[r for r in coverage if int(r['Nv'])==n]
        for row in rows:
            ax.scatter(int(row['R']),i,s=120,marker='s' if int(row['max_correlator_tau'])==32 else 'o',
                color=COLORS[n],edgecolor='white')
    ax.set_yticks(range(4));ax.set_yticklabels([str(n) for n in CUTOFFS]);ax.set_xticks(range(1,33))
    decorate(ax,r'$R/a$',r'Mode cutoff $N_v$');ax.set_xlim(.4,32.6);ax.set_ylim(-.5,3.5)
    heading(fig,'qX: available measurement coverage',
        r'Every displayed point contains the full set of 71 configurations at $p_{\rm iter}=10$',
        r'Squares: correlators through $\tau/a=32$. Circles: through 16. No symbol means unmeasured, not zero.')
    fig.subplots_adjust(top=.8,bottom=.22)
    save(fig,folder,'measurement_coverage')


def render(folder):
    folder=Path(folder)
    means=table(folder,'correlator_means');energies=table(folder,'effective_energies')
    energy_atlases(folder,energies)
    focus_signals(folder,means)
    window_plots(folder,table(folder,'constant_fit_catalogue'),table(folder,'paired_window_contrasts'))
    spatial_plots(folder,table(folder,'spatial_window_estimates'),means)
    approach_plots(folder,energies,table(folder,'approach_fits'))
    axis_plot(folder,table(folder,'axes_near_R24'))
    coverage_plot(folder,table(folder,'data_coverage'))


if __name__=='__main__':render(HERE/'results')
