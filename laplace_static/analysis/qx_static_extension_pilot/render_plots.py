#!/usr/bin/env python3
"""Render signed pilot signals and explicitly flagged effective energies."""
import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
RESULTS=HERE/'results'
CUTOFFS=[40,60,80,100]
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]
COLORS=['#0072B2','#E69F00','#009E73','#CC79A7']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.alpha':.18,'svg.hashsalt':'qx-static-extension-pilot',
    'savefig.facecolor':'white'})


def load(name):
    with (RESULTS/name).open(newline='') as stream:return list(csv.DictReader(stream))


def save(fig,name):
    for extension in ('png','svg'):
        destination=RESULTS/(name+'.'+extension)
        temporary=RESULTS/(name+'.'+extension+'.pending')
        options=dict(dpi=180) if extension=='png' else dict(metadata={'Date':None})
        fig.savefig(str(temporary),format=extension,**options)
        with temporary.open('rb') as stream:os.fsync(stream.fileno())
        temporary.replace(destination)
    plt.close(fig)


def selected(table,n,r):
    return sorted([row for row in table if int(row['Nv'])==n and int(row['R'])==r],key=lambda x:int(x['tau']))


def signals():
    table=load('correlator_diagnostics_pilot.csv')
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.3));axes=axes.ravel()
    for ax,r in zip(axes,[4,12,24,32]):
        for n,color in zip(CUTOFFS,COLORS):
            rows=selected(table,n,r);t=np.array([int(x['tau']) for x in rows])
            y=np.array([float(x['normalized_mean']) for x in rows]);err=np.array([float(x['normalized_se']) for x in rows])
            ax.plot(t,y,color=color,lw=1.5,label=r'$N_v=%d$'%n)
            ax.errorbar(t[::4],y[::4],yerr=err[::4],fmt='o',ms=3,color=color,capsize=2,lw=.8)
            negative=y<=0
            if negative.any():ax.scatter(t[negative],y[negative],s=25,marker='x',color=color,zorder=5)
        ax.axvline(16,color='.5',ls='--',lw=1)
        if r in (4,12):
            ax.set_yscale('log');ax.set_ylim(.01 if r==4 else .0008,1.3)
        else:
            ax.axhline(0,color='.3',lw=.8);ax.set_yscale('symlog',linthresh=1e-5)
        ax.set_xlim(0,32)
        ax.set_title(r'$R/a=%d$'%r);ax.set_xlabel(r'$τ/a$');ax.set_ylabel(r'$⟨L(τ)⟩/⟨L(0)⟩$')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.932),ncol=4,frameon=False)
    fig.suptitle('qX temporal extension: normalized correlators at four mode cutoffs',fontsize=15,y=.98)
    fig.text(.5,.934,'n170, n640, n650 · piter=10 · source time 90 · all spatial sources and axes',ha='center',fontsize=10)
    fig.text(.5,.018,'Technical pilot, N=3: illustrative paired naive errors every four times. Signed log at R=24,32; crosses mark nonpositive means.\n'
        'Dashed line: previous time limit. Every configuration is retained; no absolute-value transformation.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.09,right=.98,bottom=.12,top=.85,hspace=.32,wspace=.3)
    save(fig,'normalized_correlators_mode_cutoffs_pilot')


def energies():
    table=load('effective_energies_pilot.csv')
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.3));axes=axes.ravel()
    for ax,r in zip(axes,[4,12,16,24]):
        for n,color in zip(CUTOFFS,COLORS):
            rows=selected(table,n,r);t=np.array([int(x['tau'])+.5 for x in rows])
            y=np.array([float(x['effective_energy']) if x['effective_energy'] else np.nan for x in rows])
            err=np.array([float(x['naive_se']) if x['naive_se'] else np.nan for x in rows])
            strong=np.array([x['diagnostic']=='positive_endpoint_means' for x in rows])
            ax.plot(t,y,color=color,lw=1.1,label=r'$N_v=%d$'%n)
            show=np.isfinite(y)&(np.arange(len(y))%2==0)
            ax.errorbar(t[show],y[show],yerr=err[show],fmt='none',color=color,lw=.7,capsize=2)
            ax.scatter(t[strong],y[strong],s=14,color=color,zorder=4)
            weak=np.isfinite(y)&~strong
            ax.scatter(t[weak],y[weak],s=20,facecolors='white',edgecolors=color,zorder=4)
        ax.axvline(16,color='.5',ls='--',lw=1);ax.axhline(0,color='.5',lw=.7)
        ax.set_xlim(0,32);ax.set_title(r'$R/a=%d$'%r)
        ax.set_xlabel(r'$(τ+1/2)/a$');ax.set_ylabel(r'$aE_{eff}=log[⟨L(τ)⟩/⟨L(τ+1)⟩]$')
        if r==4:ax.set_ylim(.09,.15)
        elif r==12:ax.set_ylim(-.1,.65)
        elif r==16:ax.set_ylim(-.5,1.1)
        else:ax.set_ylim(-2.5,3)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.932),ncol=4,frameon=False)
    fig.suptitle('qX temporal extension: effective-energy diagnostics',fontsize=15,y=.98)
    fig.text(.5,.934,'n170, n640, n650 · piter=10 · ensemble means formed before logarithms',ha='center',fontsize=10)
    fig.text(.5,.018,'Technical pilot, N=3: paired naive errors every two times. Gaps: a nonpositive endpoint mean; hollow markers: endpoint mean/SE < 3.\n'
        'Display ranges emphasize the central trends; all numerical values and flags are exported. No plateau fits or potential claims.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.09,right=.98,bottom=.12,top=.85,hspace=.32,wspace=.3)
    save(fig,'effective_energies_mode_cutoffs_pilot')


def support():
    table=load('correlator_diagnostics_pilot.csv')
    fig,axes=plt.subplots(2,2,figsize=(11.8,8.3));axes=axes.ravel()
    for ax,n in zip(axes,CUTOFFS):
        matrix=np.array([[float(row['mean_over_se']) for row in selected(table,n,r)] for r in RADII])
        ax.grid(False)
        im=ax.imshow(matrix,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-5,vmax=5,
            extent=(-.5,32.5,-.5,11.5),interpolation='nearest')
        bad=np.where(matrix<=0)
        ax.scatter(bad[1],bad[0],s=17,marker='x',color='black',lw=.8)
        ax.axvline(16.5,color='black',ls='--',lw=1)
        ax.set_xticks(range(0,33,4));ax.set_yticks(range(12));ax.set_yticklabels(RADII)
        ax.set_xlabel(r'$τ/a$');ax.set_ylabel(r'selected $R/a$');ax.set_title(r'$N_v=%d$'%n)
    fig.suptitle('qX temporal extension: signed mean / standard error',fontsize=15,y=.98)
    fig.text(.5,.94,'n170, n640, n650 · piter=10 · one statistical sample per configuration',ha='center',fontsize=10)
    fig.subplots_adjust(left=.075,right=.865,bottom=.12,top=.88,hspace=.35,wspace=.25)
    colorbar=fig.colorbar(im,cax=fig.add_axes([.895,.22,.022,.56]),extend='both')
    colorbar.set_label('Mean / standard error (clipped at +/-5)')
    fig.text(.5,.025,'Technical pilot, N=3: ratios can be unstable and are not significance or acceptance tests.\n'
        'Crosses mark nonpositive ensemble means; negative individual measurements remain included.',ha='center',fontsize=9)
    save(fig,'signed_signal_support_pilot')


if __name__=='__main__':signals();energies();support()
