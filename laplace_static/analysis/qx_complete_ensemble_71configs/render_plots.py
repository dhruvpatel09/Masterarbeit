#!/usr/bin/env python3
"""Render the complete qX cohort with explicit rounding scales and bin support."""
from pathlib import Path

PITERS = (4, 10, 20, 40)
PC = {4:'#cd5519', 10:'#9c5f9b', 20:'#278260', 40:'#176b8e'}
PM = {4:'x', 10:'s', 20:'^', 40:'o'}
TC = {2:'#3269a8', 4:'#db7c21', 6:'#168774'}
FOOT = '71 matched configurations, n100–800 including n640; combined spatial axes'
ERRORS = 'Preliminary paired configuration jackknife; no autocorrelation correction'


def cut(records, tau, r=4, longitudinal=2):
    return sorted([v for v in records if tuple(int(v[k]) for k in
        ('tau','axis','R','long','tr2')) == (tau,-1,r,longitudinal,0)], key=lambda v:int(v['tr1']))


def render(raw, trace, windows, temporal, covariance, controls, binned, plans, influence, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    output = Path(output)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
        'axes.spines.top':False, 'axes.spines.right':False, 'axes.titlesize':12,
        'svg.hashsalt':'qx-complete-71-configurations'})

    def values(rows, field, scale=1):
        return [float(r[field])*scale for r in rows]

    def save(fig, name):
        fig.savefig(output/(name+'.png'), dpi=180, facecolor='white')
        fig.savefig(output/(name+'.svg'), metadata={'Date':None}, facecolor='white')
        plt.close(fig)

    def finish(fig, title, subtitle, footer, bottom=.12, top=.83, wspace=.16, hspace=.25):
        fig.suptitle(title, fontsize=15, y=.98)
        fig.text(.5,.89 if fig.get_figheight()<5.6 else .925,subtitle,ha='center',fontsize=10)
        fig.text(.5,.02,footer,ha='center',fontsize=9,color='#404040')
        fig.subplots_adjust(left=.10,right=.98,top=top,bottom=bottom,wspace=wspace,hspace=hspace)

    def style(ax, zero=False):
        ax.grid(alpha=.18)
        ax.set_xticks(range(-2,3))
        if zero:
            ax.axhline(0,color='#666666',lw=.8)

    # Separate the p4 and rounding-scale p10/p20 differences instead of hiding them.
    fig, axes = plt.subplots(3,3,figsize=(12.8,10),sharex=True,sharey='row')
    for col,tau in enumerate((2,4,6)):
        rows=cut(raw,tau); x=values(rows,'tr1')
        for p in reversed(PITERS):
            axes[0,col].plot(x,values(rows,'C_p%d'%p,1e3),color=PC[p],marker=PM[p],
                markerfacecolor='none',ms=7 if p==40 else 5,lw=1.3,ls='-' if p==40 else '--',label='piter=%d'%p)
        axes[1,col].plot(x,values(rows,'delta_C_p4_minus_p40',1e9),'x-',color=PC[4],label='piter=4')
        for p in (10,20):
            axes[2,col].plot(x,values(rows,'delta_C_p%d_minus_p40'%p,1e16),marker=PM[p],
                            color=PC[p],ls='--',label='piter=%d'%p)
        axes[0,col].set_title(r'$\tau/a=%d$'%tau)
        for k,ax in enumerate(axes[:,col]):
            style(ax,k>0); ax.legend(fontsize=8)
        axes[2,col].set_xlabel(r'$r_\perp/a$')
    axes[0,0].set_ylabel(r'$10^3 C_p$')
    axes[1,0].set_ylabel(r'$10^9(C_4-C_{40})$')
    axes[2,0].set_ylabel(r'$10^{16}(C_p-C_{40})$')
    finish(fig,'Raw profiles across four projection prescriptions',r'$R/a=4$, longitudinal midpoint; $p$ denotes piter',
        FOOT+'\nTop curves overlap; p40 defines zero difference. Bottom row resolves rounding-scale differences; no error bars.',top=.86,bottom=.09)
    save(fig,'profiles_all_piters_71configs')

    fig, axes = plt.subplots(2,3,figsize=(12.8,7.6),sharex=True,sharey='row')
    for col,tau in enumerate((2,4,6)):
        rows=cut(raw,tau); x=values(rows,'tr1'); y=values(rows,'C_p40',1e3); se=values(rows,'C_p40_jackknife_se',1e3)
        axes[0,col].fill_between(x,np.array(y)-se,np.array(y)+se,color=PC[40],alpha=.12)
        axes[0,col].plot(x,y,'o-',color=PC[40],label='piter=40, ±1 SE shading')
        axes[0,col].plot(x,values(rows,'C_p4',1e3),'x--',color=PC[4],label='piter=4')
        axes[1,col].errorbar(x,values(rows,'delta_C_p4_minus_p40',1e9),
            yerr=values(rows,'delta_C_paired_jackknife_se',1e9),fmt='o-',color=PC[4],capsize=3)
        axes[0,col].set_title(r'$\tau/a=%d$'%tau); axes[0,col].legend(fontsize=8)
        axes[1,col].set_xlabel(r'$r_\perp/a$')
        for k,ax in enumerate(axes[:,col]):style(ax,k>0)
    axes[0,0].set_ylabel(r'$10^3 C_p$');axes[1,0].set_ylabel(r'$10^9(C_4-C_{40})$')
    finish(fig,'Raw profiles and paired projection differences',r'$R/a=4$, longitudinal midpoint; $p$ denotes piter',
           FOOT+'\n'+ERRORS)
    save(fig,'profiles_preliminary_errors_71configs')

    configs=list(range(100,801,10))
    window_lookup={(int(r['configuration']),int(r['piter'])):[] for r in windows}
    for r in windows:window_lookup[int(r['configuration']),int(r['piter'])].append(r)
    fig,axes=plt.subplots(2,2,figsize=(12.8,7.7),sharex=True)
    for col,(field,title) in enumerate((('unit_max',r'$\|U^\dagger U-I\|_F$'),('det_max',r'$|\det U-1|$'))):
        for p in PITERS:
            y=[max(float(r[field]) for r in window_lookup[c,p]) for c in configs]
            axes[0,col].semilogy(configs,y,marker=PM[p],color=PC[p],ms=4,lw=1,label='piter=%d'%p)
            if p!=4:
                axes[1,col].plot(configs,np.array(y)*1e16,marker=PM[p],color=PC[p],
                                 ms=4,lw=.8,ls='--',label='piter=%d'%p)
        axes[0,col].set_title(title);axes[0,col].legend(fontsize=8)
        axes[1,col].legend(fontsize=8,ncol=3);axes[1,col].set_xlabel('Configuration n')
        for ax in axes[:,col]:
            ax.grid(alpha=.18);ax.axvline(640,color='#777777',ls=':',lw=.8)
            ax.set_xticks((100,200,300,400,500,600,640,700,800));ax.tick_params(axis='x',rotation=45)
    axes[0,0].set_ylabel('Largest central-window maximum')
    axes[1,0].set_ylabel(r'$10^{16}$ × defect (rounding-scale detail)')
    finish(fig,'Central-window SU(3) accuracy across the complete cohort',
        r'Largest maximum over $\tau/a=2,4,6$ per configuration; midpoint $t/a=98$',
        FOOT+'\nAll four curves are present above; lower panels separate p10/p20/p40 near roundoff. Dotted line: n640.',bottom=.14)
    save(fig,'su3_central_windows_all_piters_71configs')

    fig,axes=plt.subplots(1,2,figsize=(12.7,5.1))
    for ax,kind in zip(axes,('unit','det')):
        for p,suffix,label,color,ls,marker in ((4,'max','p4 maximum',PC[4],'-','o'),
            (4,'rms','p4 RMS','#278260','--','s'),(40,'max','p40 maximum',PC[40],':','^')):
            y=[max(float(r[kind+'_'+suffix]) for r in window_lookup[c,p]) for c in configs]
            ax.semilogy(configs,y,color=color,ls=ls,marker=marker,ms=3,lw=1,label=label)
        ax.axhline(.01,color='#888888',ls='--',lw=.7);ax.legend(fontsize=8)
        ax.set_title('Unitarity defect' if kind=='unit' else 'Determinant defect')
        ax.grid(alpha=.18);ax.set_xlabel('Configuration n');ax.set_ylim(1e-16,.05)
    axes[0].set_ylabel('Largest value over the three measured windows')
    finish(fig,'Central-link maxima and RMS',r'Measured temporal transporter windows only',
        FOOT+'\n0.01 is a diagnostic counter threshold, not an agreed physics acceptance limit.',bottom=.16,top=.79)
    save(fig,'su3_central_windows_71configs')

    fig,axes=plt.subplots(2,3,figsize=(12.8,7.5),sharex=True,sharey='row')
    for col,tau in enumerate((2,4,6)):
        rows=cut(trace,tau);x=values(rows,'tr1');color=TC[tau]
        for definition,marker,ls,label in (('RAW','s','--','Raw'),('TL','o','-','Traceless')):
            field='C_%s_p40'%definition
            axes[0,col].errorbar(x,values(rows,field,1e3),yerr=values(rows,field+'_jackknife_se',1e3),
                color=color,marker=marker,ls=ls,mfc='white' if definition=='RAW' else color,capsize=3,label=label)
        axes[1,col].errorbar(x,values(rows,'C_TRACE_p40',1e6),yerr=values(rows,'C_TRACE_p40_jackknife_se',1e6),
            fmt='o-',color=color,capsize=3)
        axes[0,col].set_title(r'$\tau/a=%d$'%tau);axes[0,col].legend(fontsize=9)
        axes[1,col].set_xlabel(r'$r_\perp/a$')
        for ax in axes[:,col]:style(ax)
    axes[0,0].set_ylabel(r'$10^3 C$');axes[1,0].set_ylabel(r'$10^6(C_{\rm RAW}-C_{\rm TL})$')
    finish(fig,'Raw and colour-traceless clover',r'$piter=40$, $R/a=4$, longitudinal midpoint; original-field clover',FOOT+'\n'+ERRORS)
    save(fig,'clover_raw_vs_traceless_71configs')

    fig,axes=plt.subplots(1,3,figsize=(12.7,4.8),sharey=True)
    for ax,(r,l) in zip(axes,((1,0),(2,1),(4,2))):
        for tau in (2,4,6):
            rows=cut(trace,tau,r,l)
            ax.errorbar(values(rows,'tr1'),values(rows,'trace_fraction_percent_p40'),
                yerr=values(rows,'trace_fraction_percent_p40_jackknife_se'),fmt='o-',color=TC[tau],capsize=3,label=r'$\tau/a=%d$'%tau)
        ax.set_title(r'$R/a=%d$, $l/a=%d$'%(r,l));ax.set_xlabel(r'$r_\perp/a$');style(ax)
    axes[0].set_ylabel(r'$100(C_{\rm RAW}-C_{\rm TL})/C_{\rm RAW}$ [%]');axes[2].legend(fontsize=9)
    finish(fig,'Relative change from colour-trace subtraction',r'$piter=40$; paired errors retain covariance between definitions',FOOT+'\n'+ERRORS,bottom=.19,top=.76)
    save(fig,'clover_trace_fraction_71configs')

    def midcut(row):return tuple(int(row[k]) for k in ('axis','R','long','tr2'))==(-1,4,2,0)

    fig,axes=plt.subplots(2,2,figsize=(11.7,8),sharex=True,sharey='row')
    for col,definition in enumerate(('RAW','TL')):
        rows=sorted([r for r in temporal if midcut(r) and r['definition']==definition and int(r['piter'])==40],key=lambda r:int(r['tr1']))
        x=values(rows,'tr1')
        for t in (2,4,6):
            field='C_tau%d'%t
            axes[0,col].errorbar(x,values(rows,field,1e3),yerr=values(rows,field+'_jackknife_se',1e3),
                fmt='o-',color=TC[t],capsize=3,label=r'$\tau/a=%d$'%t)
        for high,low,marker in ((4,2,'s'),(6,4,'o')):
            field='delta_tau%d_minus_tau%d'%(high,low)
            axes[1,col].errorbar(x,values(rows,field,1e5),yerr=values(rows,field+'_jackknife_se',1e5),
                color=TC[high],marker=marker,capsize=3,label=r'$C_{\tau=%d}-C_{\tau=%d}$'%(high,low))
        axes[0,col].set_title('Raw clover' if definition=='RAW' else 'Traceless clover')
        for k,ax in enumerate(axes[:,col]):style(ax,k>0);ax.legend(fontsize=9)
        axes[1,col].set_xlabel(r'$r_\perp/a$')
    axes[0,0].set_ylabel(r'$10^3 C$');axes[1,0].set_ylabel(r'$10^5(C_{\tau_{\rm high}}-C_{\tau_{\rm low}})$')
    finish(fig,'Paired temporal-separation comparison',r'$piter=40$, $R/a=4$, longitudinal midpoint; subscripts here denote $\tau/a$',FOOT+'\n'+ERRORS)
    save(fig,'temporal_separation_profiles_71configs')

    fig,axes=plt.subplots(1,2,figsize=(10.6,5.2))
    for ax,definition in zip(axes,('RAW','TL')):
        rows=[r for r in covariance if midcut(r) and int(r['tr1'])==0 and r['definition']==definition and int(r['piter'])==40]
        lookup={(int(r['tau_i']),int(r['tau_j'])):float(r['correlation']) for r in rows}
        matrix=np.array([[lookup[t,u] for u in (2,4,6)] for t in (2,4,6)])
        image=ax.imshow(matrix,cmap='RdBu_r',vmin=-1,vmax=1)
        for i in range(3):
            for j in range(3):ax.text(j,i,'%.3f'%matrix[i,j],ha='center',va='center',fontsize=13,color='white' if abs(matrix[i,j])>.65 else '#202020')
        ax.set_xticks(range(3));ax.set_yticks(range(3));ax.set_xticklabels((2,4,6));ax.set_yticklabels((2,4,6))
        ax.set_xlabel(r'$\tau/a$');ax.set_ylabel(r'$\tau/a$');ax.set_title('Raw clover' if definition=='RAW' else 'Traceless clover')
    fig.suptitle('Correlation between temporal profile estimates',y=.98,fontsize=15)
    fig.text(.5,.90,r'71 configurations; $piter=40$, $R/a=4$, central probe; combined axes',ha='center',fontsize=10)
    fig.text(.5,.025,'Jackknife correlations of ratio estimates; these are not Monte Carlo autocorrelations.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.80,bottom=.16,left=.07,right=.84,wspace=.3)
    cax=fig.add_axes((.89,.20,.018,.54));fig.colorbar(image,cax=cax,label='Correlation')
    save(fig,'temporal_separation_correlations_71configs')

    fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),sharey=True)
    for ax,definition in zip(axes,('RAW','TL')):
        for high,low,marker in ((4,2,'s'),(6,4,'o')):
            rows=sorted([r for r in controls if midcut(r) and r['definition']==definition
                and int(r['tau_high'])==high and int(r['tau_low'])==low],key=lambda r:int(r['tr1']))
            ax.errorbar(values(rows,'tr1'),values(rows,'contrast_p4_minus_p40',1e9),
                yerr=values(rows,'contrast_p4_minus_p40_jackknife_se',1e9),color=TC[high],marker=marker,capsize=3,
                label=r'$\Delta=C_{\tau=%d}-C_{\tau=%d}$'%(high,low))
        ax.set_title('Raw clover' if definition=='RAW' else 'Traceless clover');style(ax,True)
        ax.set_xlabel(r'$r_\perp/a$');ax.legend(fontsize=9)
    axes[0].set_ylabel(r'$10^9[\Delta(piter=4)-\Delta(piter=40)]$')
    finish(fig,'Projection control for temporal differences',r'$R/a=4$, longitudinal midpoint; paired across both $\tau$ and piter',FOOT+'\n'+ERRORS,bottom=.18,top=.77)
    save(fig,'temporal_separation_projection_control_71configs')

    fields=[('C_tau%d'%t,r'$C_{\tau=%d}$'%t,TC[t]) for t in (2,4,6)]
    differences=[('delta_tau4_minus_tau2',r'$C_{\tau=4}-C_{\tau=2}$',TC[4]),('delta_tau6_minus_tau4',r'$C_{\tau=6}-C_{\tau=4}$',TC[6])]
    fig,axes=plt.subplots(2,2,figsize=(11.7,7.8),sharex=True,sharey='row')
    for col,definition in enumerate(('RAW','TL')):
        rows=sorted([r for r in binned if midcut(r) and int(r['tr1'])==0 and r['definition']==definition],key=lambda r:int(r['bin_width']))
        for k,channels in enumerate((fields,differences)):
            ax=axes[k,col]
            for field,label,color in channels:
                ax.plot(values(rows,'bin_width'),values(rows,field+'_se_ratio'),'o-',color=color,label=label)
            ax.axhline(1,color='#666666',ls='--',lw=.8);ax.grid(alpha=.18);ax.legend(fontsize=8)
            ax.set_xscale('log',base=2);ax.set_xticks((1,2,4,8));ax.set_xticklabels(['%d\nN=%d'%(p['bin_width'],p['configurations']) for p in plans])
        axes[0,col].set_title('Raw clover' if definition=='RAW' else 'Traceless clover')
        axes[1,col].set_xlabel('Configurations per bin; retained configuration count N')
    axes[0,0].set_ylabel('Profile SE ratio\nblocked / unblocked, same cohort')
    axes[1,0].set_ylabel('Temporal-difference SE ratio\nblocked / unblocked, same cohort')
    finish(fig,'Conditional binning sensitivity at the central probe',r'$piter=40$, $R/a=4$; increasing-label ordering remains an assumption',
        'Full71 main estimates are unchanged. Complete equal-size bins retain N=71/70/68/64; bin counts 71/35/17/8.\nSame-cohort denominators; no MDU conversion, independence claim or final autocorrelation correction.',bottom=.17,top=.84,hspace=.18)
    save(fig,'binning_error_ratios_midpoint_71configs')

    support=np.array([[1 if c in p['selected_configurations'] else 2 for c in configs] for p in plans])
    fig,ax=plt.subplots(figsize=(11.8,3.6));palette=['#b7bdc5','#168774','#e59a35']
    ax.imshow(support,aspect='auto',interpolation='nearest',cmap=ListedColormap(palette),vmin=0,vmax=2,
              extent=(95,805,len(plans)-.5,-.5))
    ax.set_yticks(range(4));ax.set_yticklabels(['b=%d: N=%d, B=%d'%(p['bin_width'],p['configurations'],p['bins']) for p in plans])
    ax.set_xticks((100,200,300,400,500,600,640,700,800));ax.set_xlabel('Configuration n; n640 is included')
    ax.set_title('Configuration support for each bin width',pad=13)
    ax.legend(handles=[Patch(facecolor=palette[1],label='Included'),Patch(facecolor=palette[2],label='Incomplete end bin omitted')],
              loc='upper center',bbox_to_anchor=(.5,-.27),ncol=2,frameon=False)
    fig.text(.5,.02,'All 71 samples enter the main estimates. Larger-bin subsets retain only complete equal-size bins; every omitted label is exported.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.19,right=.98,top=.80,bottom=.37)
    save(fig,'binning_configuration_support_71configs')

    fig,axes=plt.subplots(2,3,figsize=(12.8,7.5),sharex=True,sharey='row')
    for col,tau in enumerate((2,4,6)):
        rows=cut(influence,tau);x=values(rows,'tr1')
        axes[0,col].plot(x,values(rows,'C_RAW_p40_70configs',1e3),'s--',mfc='white',color='#888888',label='70: excluding n640')
        axes[0,col].errorbar(x,values(rows,'C_RAW_p40_71configs',1e3),
            yerr=values(rows,'C_RAW_p40_71configs_jackknife_se',1e3),fmt='o-',color=PC[40],capsize=3,label='71: including n640, ±1 SE')
        axes[1,col].plot(x,values(rows,'C_RAW_p40_change_over_full71_se'),'o-',color=PC[40])
        axes[0,col].set_title(r'$\tau/a=%d$'%tau);axes[0,col].legend(fontsize=8)
        axes[1,col].set_xlabel(r'$r_\perp/a$')
        for k,ax in enumerate(axes[:,col]):style(ax,k>0)
    axes[0,0].set_ylabel(r'$10^3 C_{\rm RAW}$')
    axes[1,0].set_ylabel(r'$(C_{71}-C_{70})/\mathrm{SE}_{71}$')
    finish(fig,'Effect of adding the recovered n640 configuration',r'$piter=40$, $R/a=4$, longitudinal midpoint',
        FOOT+'\nCorrelated cohorts: the normalized change is descriptive, not an independent-sample significance test. '+ERRORS,bottom=.13)
    save(fig,'n640_inclusion_71configs')
