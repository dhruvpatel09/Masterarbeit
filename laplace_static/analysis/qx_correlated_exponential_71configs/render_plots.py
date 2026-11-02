#!/usr/bin/env python3
"""Render fit diagnostics without selecting a potential or a plateau."""
import json
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import fit_models as fm
from run_analysis import RESULTS,CUTOFFS,load_data,parameters
from derived_results import load_rows,window_for

COLORS={40:'#17639d',60:'#16857c',80:'#d57b16',100:'#7854a0'}
INK='#243746'


def save(fig,name):
    for extension in ('png','svg'):
        path=RESULTS/(name+'.'+extension);tmp=path.with_suffix(path.suffix+'.pending')
        fig.savefig(str(tmp),format=extension,dpi=190,facecolor='white',metadata={'Creator':'qX correlator analysis'})
        with tmp.open('rb') as f:os.fsync(f.fileno())
        tmp.replace(path)
    plt.close(fig)


def style():
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.labelsize':11,'axes.edgecolor':'#b4bec5',
        'axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,'grid.color':'#dce2e6',
        'grid.alpha':.65,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none',
        'legend.frameon':False,'axes.titleweight':'semibold'})


def windows(rows):
    fig,axes=plt.subplots(2,3,figsize=(13.4,8.2));fig.subplots_adjust(top=.83,bottom=.15,hspace=.36,wspace=.30)
    for ax,R in zip(axes.flat,(4,8,12,16,24,32)):
        for i,n in enumerate(CUTOFFS):
            chosen=[r for r in rows if r['states']==1 and r['Nv']==n and r['R']==R and r['tau_max']==32]
            chosen.sort(key=lambda r:r['tau_min'])
            x=np.array([r['tau_min'] for r in chosen])+(i-1.5)*.18
            ax.errorbar(x,[r['E0'] for r in chosen],yerr=[r['E0_naive_se'] for r in chosen],
                color=COLORS[n],marker=('o','s','^','D')[i],ms=4,lw=1.1,capsize=2,label=r'$N_v=%d$'%n)
        ax.set(title=r'$R/a=%d$'%R,xticks=[4,8,12,16],xlabel=r'Fit start $\tau_{\min}/a$',ylabel=r'Fitted $aE$')
        ax.grid(axis='y');ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False));ax.ticklabel_format(axis='y',style='plain')
    handles,labels=axes.flat[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.93),ncol=4)
    fig.suptitle('One exponential: sensitivity to time window and mode cutoff',fontsize=17,x=.5,y=.99)
    fig.text(.5,.945,r'71 matched configurations  |  $C(\tau)=A\,e^{-E(\tau-\tau_{\min})}$  |  $\tau_{\max}/a=32$',ha='center',fontsize=11)
    fig.text(.07,.045,'Bars: conditional independent-configuration errors with the full time covariance.\nCurves connect fit choices; they do not select a plateau. The large-distance source dependence remains visible.',fontsize=10)
    save(fig,'single_exponential_window_scan_71configs')


def states(rows):
    fig,axes=plt.subplots(2,4,figsize=(16.2,8.2));fig.subplots_adjust(top=.82,bottom=.16,hspace=.34,wspace=.31)
    for axes_row,n in zip(axes,(40,100)):
        for ax,R in zip(axes_row,(4,12,24,32)):
            for s,color,marker,shift in ((1,'#17639d','s',-.14),(2,'#b36127','o',.14)):
                chosen=sorted((r for r in rows if r['states']==s and r['Nv']==n and r['R']==R and r['tau_max']==32),key=lambda r:r['tau_min'])
                plotted=[r for r in chosen if s==1 or r['identified']]
                if plotted:
                    ax.errorbar([r['tau_min']+shift for r in plotted],[r['E0'] for r in plotted],
                        yerr=[r['E0_naive_se'] for r in plotted],fmt=marker,color=color,ms=5,capsize=2,lw=1)
                for r in chosen:
                    if s==2 and not r['identified']:
                        ax.plot(r['tau_min']+shift,.04,marker='x',color='#8b9399',ms=7,mew=1.8,transform=ax.get_xaxis_transform())
            ax.set(title=r'$N_v=%d,\ R/a=%d$'%(n,R),xticks=[4,8,12,16],xlabel=r'$\tau_{\min}/a$',ylabel=r'Fitted $aE_0$')
            ax.set_xlim(3,17);ax.grid(axis='y');ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    handles=[Line2D([],[],color='#17639d',marker='s',ls='none',label='One exponential'),
        Line2D([],[],color='#b36127',marker='o',ls='none',label='Two exponentials: local checks pass'),
        Line2D([],[],color='#8b9399',marker='x',ls='none',label='Extra state unresolved (marks at panel foot)')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.93),ncol=3)
    fig.suptitle('Does adding a second exponential constrain the lower energy?',fontsize=17,y=.99)
    fig.text(.5,.947,r'71 configurations  |  $\tau_{\max}/a=32$  |  all fitted parameters remain in the numerical exports',ha='center')
    fig.text(.055,.045,'A local parameter check is separate from goodness of fit and physical ground-state isolation.\nAt Nv=100, R/a=24, start 8, the two-state full-sample check passes but 5/71 jackknife fits fail local checks.',fontsize=10)
    save(fig,'one_vs_two_state_71configs')


def signed_checks(rows,a):
    fig=plt.figure(figsize=(13.5,10.1));outer=fig.add_gridspec(2,2,left=.075,right=.97,top=.87,bottom=.12,wspace=.24,hspace=.32)
    lookup={(r['Nv'],r['R'],r['states']):r for r in rows if r['tau_min']==8 and r['tau_max']==32}
    for cell,(n,R) in zip(outer,((40,12),(100,12),(40,24),(100,24))):
        inner=cell.subgridspec(2,1,height_ratios=[2.1,1],hspace=.05)
        ax=fig.add_subplot(inner[0]);ar=fig.add_subplot(inner[1],sharex=ax)
        row=lookup[n,R,1];w=window_for(row,a)
        ax.errorbar(w.times,w.mean,yerr=w.se,fmt='o',color=INK,ms=3,capsize=1.4,lw=.7,label='Signed mean ± naive SE')
        fine=np.linspace(0,24,250)
        for s,color,ls in ((1,'#17639d','-'),(2,'#b36127','--')):
            r=lookup[n,R,s];p=parameters(r)
            label='%d exponential%s'%(s,'' if s==1 else 's')
            if s==2 and not r['identified']:label+=' (parameters unresolved)'
            ax.plot(fine+8,fm.model(p,fine,s),color=color,ls=ls,lw=1.6,label=label)
            ar.plot(w.times,w.residual(p,s),color=color,ls=ls,marker='.',ms=4,lw=.9)
        ax.set_yscale('symlog',linthresh=1e-3,linscale=.7);ax.axhline(0,color='#8b9399',lw=.7)
        ax.set_title(r'$N_v=%d,\ R/a=%d$'%(n,R));ax.set_ylabel(r'$\langle C(\tau)\rangle/|\langle C(8)\rangle|$')
        ax.tick_params(labelbottom=False);ax.grid(axis='y');ax.legend(fontsize=8,loc='upper right')
        ar.axhline(0,color='#8b9399',lw=.8);ar.set(xlabel=r'$\tau/a$',ylabel='Whitened\nresidual',xticks=[8,12,16,20,24,28,32]);ar.grid(axis='y')
    fig.suptitle('Direct correlator fits retain signed measurements',fontsize=17,y=.99)
    fig.text(.5,.943,r'One- and two-exponential fits over $\tau/a=8\ldots32$  |  all 71 configurations  |  symmetric-log ordinate',ha='center')
    fig.text(.075,.035,'Residuals are Cholesky-whitened model minus mean, using the full time covariance.\nEach residual combines time slices; these are not independent point-by-point error bars. A fitted curve can be well determined while its two energies are not.',fontsize=10)
    save(fig,'signed_correlator_model_checks_71configs')


def calibration():
    import csv
    with (RESULTS/'bootstrap_goodness_71configs.csv').open(newline='') as f:rows=list(csv.DictReader(f))
    fig,ax=plt.subplots(figsize=(12.5,6.6));fig.subplots_adjust(top=.79,bottom=.25,left=.09,right=.97)
    x=np.arange(len(rows))
    for shift,field,color,marker,label in ((-.19,'p_gaussian_reference','#8b9399','s','Known-covariance Gaussian reference'),
            (0,'p_hotelling_linearized_reference','#17639d','D','Linearized Hotelling reference'),
            (.19,'p_recentered_bootstrap','#b36127','o','Recentered configuration bootstrap')):
        y=[float(r[field]) for r in rows]
        ax.errorbar(x+shift,y,yerr=[float(r['monte_carlo_se']) for r in rows] if 'bootstrap' in field else None,
            fmt=marker,color=color,ms=6,capsize=2,label=label)
    ax.set_yscale('log');ax.set_ylim(3e-8,1.5);ax.set_xticks(x)
    ax.set_xticklabels(['Nv %s\nR/a %s'%(r['Nv'],r['R']) for r in rows])
    ax.axhline(.05,color='#9aa4ac',ls=':',lw=1);ax.text(7.55,.055,'0.05',ha='right',fontsize=9)
    ax.set_ylabel('Upper-tail probability');ax.grid(axis='y',which='major')
    fig.suptitle('Goodness of fit depends on covariance calibration',fontsize=17,y=.99)
    fig.text(.5,.93,r'Same one-exponential fits  |  $\tau/a=8\ldots32$  |  $N=71$, 25 fitted time slices',ha='center')
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.89),ncol=3,fontsize=9)
    fig.text(.09,.06,'399 resamples per fit; the full covariance is re-estimated in every resample. Bars show finite-draw uncertainty only.\nFinite-sample bootstrap bias is not quantified. These are exploratory checks, not automatic model-acceptance thresholds.',fontsize=10)
    save(fig,'goodness_covariance_calibration_71configs')


def resolution(rows):
    fig,axes=plt.subplots(2,2,figsize=(12.3,7.7));fig.subplots_adjust(top=.82,bottom=.18,wspace=.23,hspace=.42)
    windows=[(lo,hi) for lo in (4,8,12,16) for hi in (24,32)];radii=[4,12,16,24,32]
    cmap=ListedColormap(['#32877d','#e2b468','#b7c0c7']);norm=BoundaryNorm([-.5,.5,1.5,2.5],3)
    for ax,n in zip(axes.flat,CUTOFFS):
        chosen={(r['R'],r['tau_min'],r['tau_max']):r for r in rows if r['Nv']==n and r['states']==2}
        grid=np.zeros((5,8),int)
        for i,R in enumerate(radii):
            for j,(lo,hi) in enumerate(windows):
                r=chosen[R,lo,hi]
                grid[i,j]=0 if r['identified'] else 2 if 'parameter_boundary' in r['flags'] or 'unidentified_curvature' in r['flags'] else 1
        ax.imshow(grid,cmap=cmap,norm=norm,aspect='auto')
        ax.set(title=r'$N_v=%d$'%n,xticks=np.arange(8),yticks=np.arange(5),ylabel=r'$R/a$')
        ax.set_xticklabels(['%d–%d'%v for v in windows],fontsize=9);ax.set_yticklabels(radii)
        for boundary in (1.5,3.5,5.5):ax.axvline(boundary,color='white',lw=1)
        ax.set_xlabel(r'Fitted time interval $[\tau_{\min},\tau_{\max}]/a$')
    fig.suptitle('Two-exponential parameter resolution across the scan',fontsize=17,y=.99)
    fig.text(.5,.94,'71 matched configurations  |  160 two-state fits  |  50 pass the listed local parameter checks',ha='center')
    handles=[Patch(color=cmap(i),label=s) for i,s in enumerate(('Local checks pass','Weak amplitude or gap','Boundary or singular curvature'))]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.90),ncol=3,fontsize=10)
    fig.text(.07,.045,'Convergence alone does not establish two resolved energies. The numerical tables retain every fitted result and flag.\nThe colours describe parameter diagnostics; they do not test goodness of fit or certify a ground-state plateau.',fontsize=10)
    save(fig,'two_state_parameter_resolution_71configs')


def main():
    style();rows=load_rows();a=load_data()
    windows(rows);states(rows);signed_checks(rows,a);calibration();resolution(rows)
    print('PASS rendered five exponential-fit diagnostic figures in PNG and SVG')


if __name__=='__main__':main()
