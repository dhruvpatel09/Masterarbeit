#!/usr/bin/env python3
"""Scientific figures for joint source-operator energy diagnostics."""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter
from analysis import RESULTS,CUTOFFS
from reproduce import load_rows

STYLES={'pair_dense':('#17639d','o','Nv 40/100, every 2 times'),
        'pair_matched':('#7d8992','s','Nv 40/100, every 4 times'),
        'all_four':('#bd6925','D','All four cutoffs, every 4 times')}
RADII=[4,8,12,16,24,32]


def style():
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.titleweight':'semibold','axes.spines.top':False,
        'axes.spines.right':False,'axes.edgecolor':'#b0bcc4','axes.labelcolor':'#243746','text.color':'#243746',
        'xtick.color':'#243746','ytick.color':'#243746','grid.color':'#dce2e6','grid.alpha':.7,
        'legend.frameon':False,'svg.fonttype':'none'})


def save(fig,name):
    for extension in ('png','svg'):
        path=RESULTS/(name+'.'+extension);tmp=type(path)(str(path)+'.pending')
        fig.savefig(str(tmp),format=extension,dpi=190,facecolor='white',metadata={'Creator':'qX joint source analysis'})
        with tmp.open('rb') as f:os.fsync(f.fileno())
        tmp.replace(path)
    plt.close(fig)


def window_scan(rows):
    fig,axes=plt.subplots(2,3,figsize=(13.4,8.4));fig.subplots_adjust(top=.81,bottom=.16,wspace=.30,hspace=.34)
    for ax,R in zip(axes.flat,RADII):
        for i,(design,(color,marker,label)) in enumerate(STYLES.items()):
            chosen=sorted((r for r in rows if r['R']==R and r['design']==design and r['model']=='shared_one'),key=lambda r:r['tau_min'])
            ax.errorbar(np.array([r['tau_min'] for r in chosen])+(i-1)*.16,[r['E0'] for r in chosen],
                yerr=[r['E0_naive_se'] for r in chosen],fmt=marker+'-',color=color,ms=4,capsize=2,lw=1,label=label)
        ax.set(title=r'$R/a=%d$'%R,xticks=[8,12,16],xlabel=r'Fit start $\tau_{\min}/a$',ylabel=r'Common fitted $aE$')
        ax.grid(axis='y');ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.90),ncol=3,fontsize=10)
    fig.suptitle('Shared one-state energy: source set and time-grid controls',fontsize=17,y=.99)
    fig.text(.5,.94,r'71 matched configurations  |  independent source amplitudes  |  full source/time covariance  |  $\tau_{\max}/a=32$',ha='center',fontsize=10.5)
    fig.text(.07,.045,'Bars: conditional independent-configuration errors. Every fit retains all 71 configurations.\nA precise constrained estimate does not establish source compatibility or a ground-state plateau; use the separate fit and constraint tests.',fontsize=10)
    save(fig,'shared_energy_window_scan_71configs')


def compatibility():
    with (RESULTS/'bootstrap_compatibility_71configs.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    fig,axes=plt.subplots(2,1,figsize=(12.6,8.3),sharex=True);fig.subplots_adjust(top=.79,bottom=.20,hspace=.26,left=.10,right=.97)
    for ax,prefix,title in zip(axes,('delta','shared'),('Source compatibility: common versus separate one-state energies',
                                                     'Absolute goodness of the common one-state model')):
        for i,(design,(color,marker,label)) in enumerate(STYLES.items()):
            chosen=sorted((r for r in rows if r['design']==design),key=lambda r:int(r['R']))
            y=np.array([float(r[prefix+'_p_add_one']) for r in chosen])
            lower=np.array([float(r[prefix+'_conditional_tail_low95']) for r in chosen]);upper=np.array([float(r[prefix+'_conditional_tail_high95']) for r in chosen])
            ax.errorbar(np.arange(6)+(i-1)*.16,y,yerr=np.vstack([y-lower,upper-y]),fmt=marker,color=color,ms=5,capsize=2,lw=1,label=label)
        ax.set_yscale('log');ax.set_ylim(.001,1.1);ax.axhline(.05,color='#9aa5ad',ls=':',lw=1)
        ax.text(5.4,.055,'0.05 reference',ha='right',fontsize=9);ax.set_title(title,fontsize=12)
        ax.set_ylabel('Bootstrap upper-tail\nestimate');ax.grid(axis='y',which='major')
    axes[-1].set_xticks(np.arange(6));axes[-1].set_xticklabels(RADII);axes[-1].set_xlabel(r'$R/a$')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.89),ncol=3,fontsize=10)
    fig.suptitle('Does one energy describe the different source cutoffs?',fontsize=17,y=.99)
    fig.text(.5,.94,r'$\tau_{\min}/a=12$, $\tau_{\max}/a=32$  |  399 recentered configuration resamples per case',ha='center')
    fig.text(.10,.045,'Bars: exact 95% binomial intervals for finite-draw tail probability, conditional on the empirical ensemble.\nCovariance is re-estimated in every draw. Finite-sample calibration bias and model uncertainty are not included. No automatic acceptance rule is applied.',fontsize=10)
    save(fig,'source_compatibility_bootstrap_71configs')


def model_order(rows):
    fig,axes=plt.subplots(2,3,figsize=(13.4,8.7));fig.subplots_adjust(top=.80,bottom=.16,wspace=.30,hspace=.36)
    for ax,R in zip(axes.flat,RADII):
        for i,design in enumerate(('pair_dense','all_four')):
            color=STYLES[design][0];shift=(-.24,.24)[i]
            for kind,marker in (('shared_one','s'),('shared_two','o')):
                chosen=sorted((r for r in rows if r['design']==design and r['R']==R and r['model']==kind),key=lambda r:r['tau_min'])
                shown=[r for r in chosen if kind=='shared_one' or r['local_checks_pass']]
                if shown:ax.errorbar([r['tau_min']+shift for r in shown],[r['E0'] for r in shown],
                    yerr=[r['E0_naive_se'] for r in shown],fmt=marker,color=color,ms=4.5,capsize=2,lw=1,
                    mfc=color if kind=='shared_one' else 'white')
                for r in chosen:
                    if kind=='shared_two' and not r['local_checks_pass']:
                        ax.plot(r['tau_min']+shift,.035+.04*i,marker='x',color=color,ms=6,transform=ax.get_xaxis_transform())
        ax.set(title=r'$R/a=%d$'%R,xticks=[8,12,16],xlabel=r'$\tau_{\min}/a$',ylabel=r'Common fitted $aE_0$');ax.grid(axis='y')
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    handles=[Line2D([],[],color=STYLES[d][0],marker='s',ls='none',label=STYLES[d][2]) for d in ('pair_dense','all_four')]
    handles += [Line2D([],[],color='#536b7a',marker='o',mfc='white',ls='none',label='Two states: local checks pass'),
                Line2D([],[],color='#536b7a',marker='x',ls='none',label='Two states unresolved (panel foot)')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,fontsize=10)
    fig.suptitle('Adding a shared excited state usually leaves parameters unresolved',fontsize=16.5,y=.99)
    fig.text(.5,.94,'71 configurations  |  common E0 and gap, source-specific amplitudes  |  no extra-state likelihood-ratio p-value',ha='center',fontsize=10)
    fig.text(.07,.045,'Filled squares: one state. Open circles: two states with all local parameter checks passing. Crosses mark unresolved models, not zero energy.\nOnly 5 of all 54 two-state fits pass the listed local checks. None of the selected start-12 two-state audits is treated as resolved.',fontsize=9.5)
    save(fig,'shared_one_vs_two_state_71configs')


def cutoff_energies(rows):
    fig,axes=plt.subplots(2,3,figsize=(13.4,8.2));fig.subplots_adjust(top=.83,bottom=.16,wspace=.30,hspace=.35)
    for ax,R in zip(axes.flat,RADII):
        one=next(r for r in rows if r['design']=='all_four' and r['R']==R and r['tau_min']==12 and r['model']=='shared_one')
        ind=next(r for r in rows if r['design']=='all_four' and r['R']==R and r['tau_min']==12 and r['model']=='independent_one')
        ax.axhspan(one['E0']-one['E0_naive_se'],one['E0']+one['E0_naive_se'],color='#d6e4ed',alpha=.9)
        ax.axhline(one['E0'],color='#17639d',lw=1)
        ax.errorbar(CUTOFFS,[ind['energy_Nv%d'%n] for n in CUTOFFS],yerr=[ind['energy_se_Nv%d'%n] for n in CUTOFFS],
            fmt='o',color='#bd6925',ms=4.5,capsize=2,lw=1)
        ax.set(title=r'$R/a=%d$'%R,xticks=CUTOFFS,xlabel=r'Source cutoff $N_v$',ylabel=r'Fitted $aE$');ax.grid(axis='y')
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    fig.suptitle('Allowing source-dependent energies reveals residual source effects',fontsize=16.5,y=.99)
    fig.text(.5,.94,r'All four cutoffs fitted together  |  $\tau/a=12,16,20,24,28,32$  |  same full covariance in both models',ha='center')
    fig.legend(handles=[Line2D([],[],color='#bd6925',marker='o',ls='none',label='Independent energy for each source'),
        Patch(facecolor='#d6e4ed',edgecolor='#17639d',label='Constrained common energy ± naive SE')],loc='upper center',bbox_to_anchor=(.5,.90),ncol=2,fontsize=10)
    fig.text(.07,.045,'All markers are correlated. Overlapping individual bars do not test equality; paired energy contrasts retain the source covariance.\nThese are finite-window model estimates, not four different physical static potentials.',fontsize=10)
    save(fig,'source_dependent_energies_71configs')


def source_shapes():
    with (RESULTS/'source_shape_diagnostics_71configs.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    fig,axes=plt.subplots(1,2,figsize=(13.4,6.4));fig.subplots_adjust(top=.76,bottom=.23,wspace=.25,left=.075,right=.97)
    colors=('#17639d','#16857c','#d57b16','#7854a0')
    for n,color in zip(CUTOFFS,colors):
        for ax,tau in zip(axes,(0,24)):
            chosen=sorted((r for r in rows if int(r['Nv'])==n and int(r['numerator_tau'])==tau),key=lambda r:int(r['R']))
            ax.errorbar([int(r['R']) for r in chosen],[float(r['ratio_of_means']) for r in chosen],
                yerr=[float(r['paired_naive_se']) for r in chosen],fmt='o-',color=color,ms=3.8,capsize=1.8,lw=1.1,label=r'$N_v=%d$'%n)
    axes[0].set_yscale('log');axes[0].set(title='Equal-time correlator shape',ylabel=r'$\langle C_N(R,0)\rangle/\langle C_N(1,0)\rangle$')
    axes[1].set_yscale('symlog',linthresh=1e-5,linscale=.7)
    axes[1].set(title=r'Signed correlator at $\tau/a=24$',ylabel=r'$\langle C_N(R,24)\rangle/\langle C_N(R,0)\rangle$')
    axes[1].axhline(0,color='#9aa5ad',lw=.8)
    for ax in axes:ax.set_xlabel(r'$R/a$');ax.set_xticks([1,4,8,12,16,20,24,28,32]);ax.grid(axis='y')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.88),ncol=4)
    fig.suptitle('Source cutoff changes spatial shape and long-time signal',fontsize=17,y=.99)
    fig.text(.5,.93,'All 71 configurations  |  twelve measured distances  |  ratios of ensemble means with paired naive errors',ha='center',fontsize=10.5)
    fig.text(.075,.055,'Equal-time shape is an operator diagnostic, not a direct measurement of ground-state overlap. Lines only connect measured distances.\nThe right panel retains signs and uncertainty: a source-shape minimum is not evidence for a feature in the physical potential.',fontsize=10)
    save(fig,'source_shape_and_long_time_signal_71configs')


def main():
    style();rows=load_rows();window_scan(rows);compatibility();model_order(rows);cutoff_energies(rows);source_shapes()
    print('PASS five shared-energy diagnostic figures in PNG and SVG')


if __name__=='__main__':main()
