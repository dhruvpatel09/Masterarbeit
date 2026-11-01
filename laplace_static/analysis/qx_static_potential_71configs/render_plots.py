#!/usr/bin/env python3
"""Scientific figures for the complete fixed-source qX static study."""


def render(means,energies,fits,stability,cornell,spatial,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np
    from pathlib import Path
    output=Path(output)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.hashsalt':'qx-static-full71'})
    colors=('#3269a8','#df8422','#168774','#945b9e','#c94f59','#67634b')
    def yes(v):return v is True or str(v)=='True'
    def save(fig,name):
        fig.savefig(output/(name+'.png'),dpi=200,facecolor='white')
        fig.savefig(output/(name+'.svg'),metadata={'Date':None},facecolor='white');plt.close(fig)
    def selected(rows,R=None):return sorted([r for r in rows if int(r['axis'])==-1 and (R is None or int(r['R'])==R)],key=lambda r:int(r.get('tau',r.get('tau_min',0))))
    fig,axes=plt.subplots(1,2,figsize=(12,5.3),gridspec_kw={'width_ratios':[1.05,1]})
    for R,c in zip((1,4,8,16,24,32),colors):
        rows=selected(means,R);axes[0].errorbar([int(r['tau']) for r in rows],[float(r['L_RE']) for r in rows],yerr=[float(r['L_RE_naive_se']) for r in rows],color=c,marker='o',ms=3,capsize=2,lw=1.3,label=r'$R/a=%d$'%R)
    axes[0].set_yscale('log');axes[0].set_xlabel(r'$\tau/a$');axes[0].set_ylabel(r'$\mathrm{Re}\,\overline{L}(R,\tau)$')
    axes[0].grid(alpha=.2,which='both');axes[0].legend(ncol=2);axes[0].set_title('Correlators with independent-configuration errors')
    grid=np.full((17,32),np.nan)
    for r in selected(means):grid[int(r['tau']),int(r['R'])-1]=float(r['mean_over_se'])
    image=axes[1].imshow(grid,origin='lower',aspect='auto',extent=(.5,32.5,-.5,16.5),norm=LogNorm(vmin=1,vmax=10000),cmap='viridis')
    axes[1].contour(np.arange(1,33),np.arange(17),grid,levels=[5],colors='white',linewidths=1)
    axes[1].set_xlabel(r'$R/a$');axes[1].set_ylabel(r'$\tau/a$');axes[1].set_xticks((1,8,16,24,32));axes[1].set_title('Mean / standard error across 71 configurations')
    fig.colorbar(image,ax=axes[1],fraction=.046,pad=.025,extend='max',label=r'$\overline{L}/\mathrm{SE}_{71}$')
    fig.suptitle('qX static correlators: complete 71-configuration ensemble',fontsize=15,y=.98)
    fig.text(.5,.895,r'$piter=10$, $N_v=100$, $t_0/a=90$; combined spatial axes',ha='center')
    fig.text(.5,.025,'White contour: mean / SE = 5 (diagnostic signal threshold). Source sites and axes are averaged within configurations.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.78,bottom=.15,left=.08,right=.97,wspace=.28);save(fig,'correlator_signals_71configs')

    fig,axes=plt.subplots(3,2,figsize=(11.5,10.5),sharex=True)
    for ax,Rs in zip(axes.flat,((1,2,4),(6,8,10),(12,14,16),(18,20,22),(23,24),(26,28,32))):
        ax.axvspan(6,11,color='#dce5ef',alpha=.55,zorder=0)
        for R,c in zip(Rs,colors):
            rows=selected(energies,R);valid=[r for r in rows if yes(r['valid_point'])]
            ax.errorbar([int(r['tau']) for r in valid],[float(r['aE_eff']) for r in valid],yerr=[float(r['naive_se']) for r in valid],color=c,marker='o',ms=3,capsize=2,lw=1.1,label=r'$R/a=%d$'%R)
        ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=len(Rs));ax.set_xticks((0,4,8,12,15))
    for ax in axes[-1]:ax.set_xlabel(r'Lower temporal separation $\tau/a$')
    for ax in axes[:,0]:ax.set_ylabel(r'$aE_{\rm eff}=\log(\overline{L}_{\tau}/\overline{L}_{\tau+1})$')
    fig.suptitle('Effective energies from consecutive times',y=.98,fontsize=15)
    fig.text(.5,.935,'71 configurations; full paired covariance retained; shaded region is the reference fit window 6–11',ha='center')
    fig.text(.5,.025,'All spatial separations are exported. Selected curves illustrate plateau approach, later-time drift and weak long-distance signals.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.89,bottom=.09,left=.10,right=.98,hspace=.18,wspace=.18);save(fig,'effective_energies_71configs')

    fig,axes=plt.subplots(1,2,figsize=(11.8,5.2))
    for window,c,marker in zip(((4,9),(6,11),(8,13),(10,15)),colors,('s','o','^','D')):
        rows=sorted([r for r in fits if int(r['axis'])==-1 and (int(r['tau_min']),int(r['tau_max']))==window and r['status']=='fit'],key=lambda r:int(r['R']))
        axes[0].errorbar([int(r['R']) for r in rows],[float(r['aV']) for r in rows],yerr=[float(r['naive_se']) for r in rows],color=c,marker=marker,ms=3,capsize=2,lw=1.1,label=r'$\tau/a=%d\ldots%d$'%window)
    axes[0].set_xlabel(r'$R/a$');axes[0].set_ylabel('Correlated constant-fit energy');axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[0].set_title('Reference-window estimates across the measured range')
    for window,c,marker in zip(((4,9),(8,13)),(colors[0],colors[2]),('s','^')):
        rows=sorted([r for r in stability if int(r['axis'])==-1 and (int(r['comparison_min']),int(r['comparison_max']))==window and r['z']!=''],key=lambda r:int(r['R']))
        axes[1].plot([int(r['R']) for r in rows],[float(r['z']) for r in rows],color=c,marker=marker,ms=3,label='6–11 minus %d–%d'%window)
    for z in (-2,0,2):axes[1].axhline(z,color='.5',ls='--',lw=.8)
    axes[1].set_xlabel(r'$R/a$');axes[1].set_ylabel('Window shift / paired standard error');axes[1].legend(fontsize=8);axes[1].grid(alpha=.2);axes[1].set_title('Correlations between fit windows included')
    fig.suptitle('Temporal-window stability of the candidate static energies',fontsize=14,y=.98)
    fig.text(.5,.025,'A compatible constant fit in one window does not establish absence of later-time bias. Horizontal ±2 lines are descriptive.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.83,bottom=.16,left=.08,right=.98,wspace=.23);save(fig,'temporal_window_stability_71configs')

    fig,axes=plt.subplots(1,2,figsize=(11.8,5.2))
    rows=sorted([r for r in fits if int(r['axis'])==-1 and (int(r['tau_min']),int(r['tau_max']))==(6,11) and r['status']=='fit'],key=lambda r:int(r['R']))
    axes[0].errorbar([int(r['R']) for r in rows],[float(r['aV']) for r in rows],yerr=[float(r['naive_se']) for r in rows],color=colors[0],marker='o',ms=4,capsize=2,lw=1,label='Reference window 6–11')
    late=sorted([r for r in fits if int(r['axis'])==-1 and (int(r['tau_min']),int(r['tau_max']))==(10,15) and r['status']=='fit'],key=lambda r:int(r['R']))
    axes[0].errorbar([int(r['R']) for r in late],[float(r['aV']) for r in late],yerr=[float(r['naive_se']) for r in late],color=colors[2],marker='^',ms=3,capsize=2,lw=1,label='Later window 10–15')
    axes[0].set_xlabel(r'$R/a$');axes[0].set_ylabel('Candidate energy in lattice units');axes[0].legend(fontsize=9);axes[0].grid(alpha=.2)
    for interval,c in zip(((4,12),(5,12),(6,16)),colors):
        rows=sorted([r for r in cornell if (int(r['R_min']),int(r['R_max']))==interval and r['status']=='fit'],key=lambda r:int(r['tau_min']))
        x=[int(r['tau_min']) for r in rows];y=[float(r['sigma_a2']) for r in rows];err=[float(r['sigma_a2_se']) for r in rows]
        axes[1].errorbar(x,y,yerr=err,color=c,lw=1.2,capsize=3,label=r'$R/a=%d\ldots%d$'%interval)
        for xx,yy,r in zip(x,y,rows):axes[1].plot(xx,yy,'o',color=c,mfc=c if float(r['p_gaussian_approx'])>=.05 else 'white',ms=6)
    axes[1].set_xticks((4,6,8,10));axes[1].set_xticklabels(('4–9','6–11','8–13','10–15'));axes[1].set_xlabel(r'Temporal fit window $\tau/a$')
    axes[1].set_ylabel(r'Cornell candidate $\sigma a^2$');axes[1].legend(fontsize=9);axes[1].grid(alpha=.2)
    fig.suptitle('Potential and Cornell fit-range diagnostics',fontsize=15,y=.98)
    fig.text(.5,.025,'Cornell: V0 + sigma R − alpha/R. Filled symbols have approximate p ≥ 0.05; temporal and range systematics remain.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.83,bottom=.16,left=.08,right=.98,wspace=.24);save(fig,'potential_and_cornell_diagnostics_71configs')

    windows=sorted({(int(r['tau_min']),int(r['tau_max'])) for r in fits});ranges=sorted({(int(r['R_min']),int(r['R_max'])) for r in cornell});cw=sorted({(int(r['tau_min']),int(r['tau_max'])) for r in cornell})
    grid=np.full((len(windows),32),np.nan);cg=np.full((len(ranges),len(cw)),np.nan)
    for r in fits:
        if int(r['axis'])==-1 and r['status']=='fit':grid[windows.index((int(r['tau_min']),int(r['tau_max']))),int(r['R'])-1]=max(1e-8,float(r['p_gaussian_approx']))
    for r in cornell:
        if r['status']=='fit':cg[ranges.index((int(r['R_min']),int(r['R_max']))),cw.index((int(r['tau_min']),int(r['tau_max'])))]=max(1e-8,float(r['p_gaussian_approx']))
    fig,axes=plt.subplots(1,2,figsize=(12,7.6),gridspec_kw={'width_ratios':[1.7,1]})
    for ax,data in zip(axes,(grid,cg)):image=ax.imshow(data,aspect='auto',origin='upper',norm=LogNorm(vmin=1e-8,vmax=1),cmap='viridis')
    axes[0].set_yticks(range(len(windows)));axes[0].set_yticklabels(['%d–%d'%w for w in windows]);axes[0].set_ylabel(r'Temporal window $\tau/a$')
    axes[0].set_xticks((0,7,15,23,31));axes[0].set_xticklabels(('1','8','16','24','32'));axes[0].set_xlabel(r'$R/a$');axes[0].set_title('Correlated constant-energy fits')
    axes[1].set_yticks(range(len(ranges)));axes[1].set_yticklabels(['%d–%d'%w for w in ranges]);axes[1].set_ylabel(r'Cornell spatial range $R/a$')
    axes[1].set_xticks(range(len(cw)));axes[1].set_xticklabels(['%d–%d'%w for w in cw],rotation=45);axes[1].set_xlabel(r'Temporal window $\tau/a$');axes[1].set_title('Correlated Cornell fits')
    cax=fig.add_axes((.91,.24,.014,.56));fig.colorbar(image,cax=cax,label='Approximate Gaussian p')
    fig.suptitle('Fit-quality scans on the complete qX sample',fontsize=15,y=.98)
    fig.text(.5,.025,'Color floor: 1e−8. Estimated covariance, 71 samples. Scans are exploratory; p values are not multiple-comparison corrected.',ha='center',fontsize=8.5)
    fig.subplots_adjust(top=.88,bottom=.20,left=.08,right=.86,wspace=.38);save(fig,'fit_quality_scans_71configs')

    fig,ax=plt.subplots(figsize=(9,5.4))
    for tau,c in zip((0,4,8,16),colors):
        rows=sorted([r for r in spatial if int(r['axis'])==-1 and int(r['tau'])==tau],key=lambda r:int(r['R']))
        ax.plot([int(r['R']) for r in rows],[float(r['correct_over_ignore_regional_covariance']) for r in rows],marker='o',ms=5,color=c,label=r'$\tau/a=%d$'%tau)
    ax.axhline(1,color='.4',ls='--',lw=1);ax.grid(alpha=.2);ax.legend(ncol=2);ax.set_xlabel(r'$R/a$')
    ax.set_ylabel('Configuration SE / SE ignoring regional covariance')
    fig.suptitle('Spatial averaging and correlations between regions',y=.98,fontsize=14)
    fig.text(.5,.90,'71 configuration samples; 64 spatial cubes × 3 axes averaged within each sample',ha='center',fontsize=10)
    fig.text(.5,.025,'Ratios can lie above or below one. Regional channels are correlated diagnostics, not extra independent configurations.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.80,bottom=.15,left=.12,right=.98);save(fig,'regional_covariance_diagnostics_71configs')


if __name__=='__main__':
    import reproduce as a
    render(*[a.read_csv(a.RESULTS/name) for name in (a.EXACT[0],a.EXACT[1],a.FLOAT[0],a.FLOAT[1],a.FLOAT[3],a.FLOAT[4])],output=a.RESULTS)
