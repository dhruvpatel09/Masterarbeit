#!/usr/bin/env python3
"""Render the uncertainty and stability diagnostics from archived CSV exports."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import HERE,RADII,csv_rows,read

OUT=HERE/'results'
COLORS={'pair':'#0072B2','triple':'#AA4499','full':'#D55E00'}
MARKERS={'pair':'o','triple':'v','full':'s'}
LINES={'pair':'-','triple':':','full':'--'}
LABELS={'pair':'pair (40,100)','triple':'triple (40,60,100)','full':'full (40,60,80,100)'}
WINDOWS=[6,10,14,18,22]
WLABELS=['6–12','10–16','14–20','18–24','22–28']
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,
    'legend.fontsize':9,'savefig.dpi':180,'svg.fonttype':'none','axes.spines.top':False,
    'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'axes.axisbelow':True})


def load(name):return csv_rows(OUT/(name+'.csv'))
def select(rows,**kw):return [r for r in rows if all(str(r[k])==str(v) for k,v in kw.items())]
def values(rows,name):return np.array([float(r[name]) if r[name]!='' else np.nan for r in rows])
def ordered(rows,**kw):return sorted(select(rows,**kw),key=lambda d:int(d['start_tau']))
def style(basis):
    return dict(color=COLORS[basis],marker=MARKERS[basis],linestyle=LINES[basis],
        markerfacecolor='white' if basis!='full' else COLORS[basis],markersize=4,
        linewidth=1.25,capsize=2,label=LABELS[basis])


def save(fig,name,title,caption,handles=None):
    fig.suptitle(title,fontsize=15,y=.98)
    if handles:fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.94),ncol=len(handles),frameon=False)
    fig.text(.5,.022,caption,ha='center',va='bottom',fontsize=9,linespacing=1.35)
    fig.tight_layout(rect=(.015,.10,.99,.88 if handles else .93),h_pad=1.8,w_pad=2)
    for ext in ('png','svg'):fig.savefig(str(OUT/(name+'.'+ext)),facecolor='white')
    plt.close(fig)


def basis_handles():return [Line2D([0],[0],**{k:v for k,v in style(b).items() if k!='capsize'}) for b in ('pair','full')]


def complex_family():
    rows=load('complex_family_summary');summary=read(OUT/'SUMMARY.json')['complex_component_tests']
    fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True,sharey=True)
    for ax,r in zip(axes.flat,RADII):
        a=sorted(select(rows,R=r),key=lambda d:int(d['tau']))
        ax.plot(values(a,'tau'),values(a,'maximum_absolute_pointwise_t'),'o-',color='#0072B2',label='largest component |t| at this R,time')
        ax.axhline(summary['bootstrap_95percent_max_t'],color='#D55E00',ls='--',label='95% global max-|t| threshold')
        ax.set(title=r'$R/a={}$'.format(r),xlabel=r'correlator separation $\tau/a$',ylabel='maximum absolute component t',ylim=(0,4.5))
        ax.set_xticks([4,8,12,16,24,32])
    handles,labels=axes[0,0].get_legend_handles_labels()
    save(fig,'antihermitian_family_71configs','Anti-Hermitian components: a joint uncertainty check',
        '71 configurations; 324 prescribed components; 1,999 recentered, studentized configuration draws.\n'
        'Global bootstrap tail = {:.4f}; no family-adjusted rejection at 5%. Non-rejection does not prove Hermiticity.'.format(summary['joint_null_tail']),handles)


def phases():
    rows=load('projected_correlators');fig,axes=plt.subplots(2,4,figsize=(14,7),sharex=True)
    for col,r in enumerate(RADII):
        for b in ('pair','full'):
            a=sorted(select(rows,R=r,basis=b,anchor=4,variant='alternating'),key=lambda d:int(d['tau']))
            t=values(a,'tau')
            axes[0,col].errorbar(t,values(a,'phase'),yerr=values(a,'phase_jackknife_se'),**style(b))
            if r!=32:axes[1,col].plot(t,values(a,'imaginary_pointwise_t'),**{k:v for k,v in style(b).items() if k!='capsize'})
        for sign in (-1,1):axes[0,col].axhline(sign*.05,color='#777777',lw=.7,ls=':')
        axes[0,col].set(title=r'$R/a={}$'.format(r),ylabel='phase (radians)')
        axes[1,col].axhline(0,color='black',lw=.7)
        axes[1,col].set(xlabel=r'$\tau/a$',ylabel=r'Im $Q$ / paired jackknife SE',xticks=[0,8,16,24,32])
        if r==32:axes[1,col].text(.5,.5,'Imaginary part is\nroundoff only at half box;\nno significance test',ha='center',va='center',transform=axes[1,col].transAxes)
    save(fig,'projected_phase_uncertainty_71configs','Projected correlator phases and imaginary uncertainties',
        'Original alternating held-out split; training time 4a; all 71 configurations with complete source retraining in every deletion.\n'
        'Dotted ±0.05 rad lines are visual guides. Pointwise imaginary t values have no family-adjusted interpretation.',basis_handles())


def split_shifts():
    rows=load('paired_fit_contrasts');fig,axes=plt.subplots(2,2,figsize=(11,8))
    colors=['#0072B2','#D55E00','#AA4499','#444444'];markers=['o','s','^','D'];handles=[]
    for j,seed in enumerate((1701,1702,1703,1704)):
        handles.append(Line2D([0],[0],color=colors[j],marker=markers[j],label='split seed '+str(seed),ls=['-','--',':','-.'][j]))
        for ax,r in zip(axes.flat,RADII):
            a=ordered(rows,R=r,basis='full',anchor=4,variant='seed'+str(seed),contrast='minus_alternating')
            ax.errorbar(np.arange(5)+(j-1.5)*.05,values(a,'difference'),yerr=values(a,'paired_jackknife_se'),
                color=colors[j],marker=markers[j],ls=['-','--',':','-.'][j],ms=4,capsize=2,lw=1)
    for ax,r in zip(axes.flat,RADII):
        ax.axhline(0,color='black',lw=.8);ax.set(title=r'$R/a={}$'.format(r),xticks=range(5),xticklabels=WLABELS,
            xlabel='seven-point correlator fit window',ylabel=r'$aE_{\rm split}-aE_{\rm alternating}$')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    save(fig,'heldout_split_energy_shifts_71configs','Changing the held-out split: paired energy shifts',
        'Full source basis; training time 4a. Error bars pair the same deleted configuration in both fits.\n'
        'All five splits reuse the same 71 configurations; their spread is not an independent-sample uncertainty.',handles)


def window_scan():
    rows=load('correlator_fits');fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,r in zip(axes.flat,RADII):
        for j,b in enumerate(('pair','full')):
            a=ordered(rows,R=r,basis=b,anchor=4,variant='alternating');x=np.arange(5)+(j-.5)*.09;y=values(a,'aE')
            ax.errorbar(x,y,yerr=values(a,'energy_jackknife_se'),**style(b))
            questionable=(values(a,'minimum_endpoint_signal')<3)|(values(a,'maximum_endpoint_phase')>.05)
            poor=values(a,'gaussian_p_approx')<.05
            ax.plot(x[questionable],y[questionable],'x',color='black',ms=9,mew=1.2,ls='none')
            ax.plot(x[poor],y[poor],'s',mfc='none',mec='black',ms=10,ls='none')
        ax.set(title=r'$R/a={}$'.format(r),xticks=range(5),xticklabels=WLABELS,xlabel='correlator fit window',ylabel=r'fitted $aE$')
    handles=basis_handles()+[Line2D([0],[0],marker='x',color='black',ls='none',label='phase >0.05 rad or real signal <3'),
        Line2D([0],[0],marker='s',mfc='none',mec='black',ls='none',label='approx. Gaussian tail <0.05')]
    save(fig,'correlator_fit_window_scan_71configs','Direct correlated one-exponential fits across time windows',
        'Original alternating held-out split; training time 4a; 71-configuration complete-retraining jackknife errors.\n'
        'Fit tails are uncalibrated diagnostics. A satisfactory local tail does not establish a time-independent energy.',handles)


def estimator_difference():
    rows=load('correlator_fits');fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,r in zip(axes.flat,RADII):
        for j,b in enumerate(('pair','full')):
            a=ordered(rows,R=r,basis=b,anchor=4,variant='alternating')
            ax.errorbar(np.arange(5)+(j-.5)*.08,values(a,'direct_minus_constant'),yerr=values(a,'paired_direct_minus_constant_se'),**style(b))
        ax.axhline(0,color='black',lw=.8);ax.set(title=r'$R/a={}$'.format(r),xticks=range(5),xticklabels=WLABELS,
            xlabel='correlator fit window',ylabel=r'$aE_{\rm correlator}-aE_{\rm effective\ constant}$')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    save(fig,'correlator_vs_effective_energy_71configs','Direct fits versus effective-energy constants',
        'Seven correlators at t…t+6 match six log ratios at t…t+5. Both estimators use the same measurements.\n'
        'Differences and errors are paired through all 71 deletions; a tiny paired error is not a second energy determination.',basis_handles())


def basis_training():
    rows=load('paired_fit_contrasts');fig,axes=plt.subplots(2,4,figsize=(14,7.5),sharex=True)
    for col,r in enumerate(RADII):
        for j,b in enumerate(('triple','full')):
            a=ordered(rows,R=r,basis=b,anchor=4,variant='alternating',contrast='minus_pair_basis')
            axes[0,col].errorbar(np.arange(5)+(j-.5)*.06,values(a,'difference'),yerr=values(a,'paired_jackknife_se'),**style(b))
        for j,anchor in enumerate((6,8)):
            a=ordered(rows,R=r,basis='full',anchor=anchor,variant='alternating',contrast='minus_anchor4')
            axes[1,col].errorbar(np.arange(5)+(j-.5)*.06,values(a,'difference'),yerr=values(a,'paired_jackknife_se'),
                color=['#0072B2','#D55E00'][j],marker=['o','s'][j],ls=['-','--'][j],ms=4,capsize=2,label='training '+str(anchor)+'a minus 4a')
        axes[0,col].set(title=r'$R/a={}$'.format(r),ylabel=r'$aE_{\rm basis}-aE_{\rm pair}$')
        axes[1,col].set(ylabel=r'$aE_{\rm training}-aE_{4a}$',xlabel='correlator fit window')
        for ax in axes[:,col]:
            ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(5));ax.set_xticklabels(WLABELS,rotation=30)
            ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    handles=[Line2D([0],[0],**{k:v for k,v in style(b).items() if k!='capsize'}) for b in ('triple','full')]
    handles.extend(axes[1,0].get_legend_handles_labels()[0])
    save(fig,'paired_basis_and_training_checks_71configs','Source basis and training-time sensitivity',
        'Original alternating held-out split. Upper row: basis change at training 4a; lower row: training change in the full basis.\n'
        'Every contrast is paired over identical configuration deletions. These are source controls, not independent ensembles.',handles)


def signed_checks():
    corrs=load('projected_correlators');fits=load('correlator_fits');fig,axes=plt.subplots(2,2,figsize=(11,8))
    colors=['#0072B2','#D55E00','#AA4499'];lines=['-','--',':'];handles=[Line2D([0],[0],color='black',marker='o',ls='none',label='Re Q, jackknife error')]
    handles.extend(Line2D([0],[0],color=colors[i],ls=lines[i],label='fit '+str(t)+'–'+str(t+6)) for i,t in enumerate((6,14,22)))
    for ax,r in zip(axes.flat,RADII):
        a=sorted(select(corrs,R=r,basis='full',anchor=4,variant='alternating'),key=lambda d:int(d['tau']))
        reference=float(a[6]['mean_re']);a=a[4:]
        ax.errorbar(values(a,'tau'),values(a,'mean_re')/reference,yerr=values(a,'jackknife_se_re')/abs(reference),
            color='black',fmt='o',ms=3,capsize=2,lw=.8)
        for i,t in enumerate((6,14,22)):
            fit=select(fits,R=r,basis='full',anchor=4,variant='alternating',start_tau=t)[0]
            times=np.linspace(t,t+6,80);model=float(fit['amplitude'])*np.exp(-float(fit['aE'])*(times-t))
            ax.plot(times,model/reference,color=colors[i],ls=lines[i],lw=1.8)
        ax.set_yscale('symlog',linthresh=.001);ax.axhline(0,color='#777777',lw=.7)
        ax.set(title=r'$R/a={}$'.format(r),xlabel=r'$\tau/a$',ylabel=r'Re $Q$ / fixed full-sample Re $Q(6)$',xlim=(3,33))
    save(fig,'signed_correlator_fit_checks_71configs','Signed correlators and fitted one-state shapes',
        'Full basis; training 4a; original alternating held-out split. Signed symlog axis retains zero and negative values.\n'
        'A fixed full-sample display normalization scales both points and errors; model lines appear only within their fit windows.',handles)


def time_changes():
    rows=load('paired_fit_contrasts');fig,axes=plt.subplots(1,2,figsize=(11,5.5))
    for ax,start in zip(axes,(14,22)):
        for j,b in enumerate(('pair','full')):
            a=select(rows,basis=b,anchor=4,variant='alternating',start_tau=start,contrast='minus_window6_12')
            a=sorted(a,key=lambda d:int(d['R']))
            ax.errorbar(np.arange(4)+(j-.5)*.08,values(a,'difference'),yerr=values(a,'paired_jackknife_se'),**style(b))
        ax.axhline(0,color='black',lw=.8);ax.set(title='{}–{} minus 6–12'.format(start,start+6),
            xticks=range(4),xticklabels=RADII,xlabel=r'$R/a$',ylabel=r'paired time-window change in $aE$')
    save(fig,'paired_time_window_changes_71configs','Time-window dependence survives held-out source training',
        'Training time 4a; original alternating held-out split; errors pair the same 71 configuration deletions.\n'
        'At R/a=32 the late-minus-early full-basis shift is −0.1359 ± 0.0396; this is a diagnostic, not a potential curve.',basis_handles())


if __name__=='__main__':
    for function in (complex_family,phases,split_shifts,window_scan,estimator_difference,basis_training,signed_checks,time_changes):function()
    print('Rendered eight uncertainty/stability figures in PNG and SVG.')
