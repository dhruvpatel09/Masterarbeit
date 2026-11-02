#!/usr/bin/env python3
"""Static scientific figures for the matched full71 temporal extension."""
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
SHOW_R=[4,8,12,16,24,32]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.alpha':.18,'svg.hashsalt':'qx-static-extension-full71',
    'savefig.facecolor':'white'})


def read(name):
    with (RESULTS/name).open(newline='') as stream:return list(csv.DictReader(stream))


def number(row,key):return float(row[key]) if row[key] else np.nan


def save(fig,name):
    for extension in ('png','svg'):
        target=RESULTS/(name+'.'+extension);pending=RESULTS/(name+'.'+extension+'.pending')
        options=dict(dpi=180) if extension=='png' else dict(metadata={'Date':None})
        fig.savefig(str(pending),format=extension,**options)
        with pending.open('rb') as stream:os.fsync(stream.fileno())
        pending.replace(target)
    plt.close(fig)


def panel_figure(title,subtitle,footnote):
    fig,axes=plt.subplots(2,3,figsize=(14,9))
    fig.suptitle(title,fontsize=16,y=.98)
    fig.text(.5,.938,subtitle,ha='center',fontsize=10)
    fig.text(.5,.02,footnote,ha='center',fontsize=9)
    fig.subplots_adjust(left=.07,right=.98,bottom=.115,top=.855,hspace=.34,wspace=.3)
    return fig,axes.ravel()


def legend(fig,axes):
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',bbox_to_anchor=(.5,.93),ncol=4,frameon=False)


def capped_points(ax,x,y,error,color,strong=None,every=1,offset=0.,label=None):
    x=np.asarray(x,float)+offset;y=np.asarray(y,float);error=np.asarray(error,float)
    finite=np.isfinite(y)
    if strong is None:strong=np.ones(len(y),bool)
    strong=np.asarray(strong,bool)&finite
    low,high=ax.get_ylim();inside=finite&(y>=low)&(y<=high)
    ax.plot(x,y,color=color,lw=1.2,label=label)
    chosen=inside&(np.arange(len(y))%every==0)
    ax.errorbar(x[chosen],y[chosen],yerr=error[chosen],fmt='none',lw=.8,capsize=2,color=color)
    ax.scatter(x[inside&strong],y[inside&strong],s=15,color=color,zorder=4)
    weak=inside&~strong
    ax.scatter(x[weak],y[weak],s=22,facecolors='white',edgecolors=color,zorder=4)
    ax.scatter(x[finite&(y<low)],np.full(np.sum(finite&(y<low)),low+.018*(high-low)),
        marker='v',s=22,color=color,zorder=5)
    ax.scatter(x[finite&(y>high)],np.full(np.sum(finite&(y>high)),high-.018*(high-low)),
        marker='^',s=22,color=color,zorder=5)


def normalized_signals():
    data=read('correlator_means_71configs.csv')
    fig,axes=panel_figure('qX: normalized correlators across four mode cutoffs',
        '71 configurations · piter=10 · source time 90 · all spatial sources and axes',
        'Paired naive errors every four times. Signed logarithmic scale at R=24,32; crosses: nonpositive means.\n'
        'Dashed line: previous time limit. The different cutoffs use the same 71 configuration samples.')
    for ax,R in zip(axes,SHOW_R):
        for n,color in zip(CUTOFFS,COLORS):
            rows=sorted([r for r in data if int(r['Nv'])==n and int(r['R'])==R],key=lambda r:int(r['tau']))
            x=np.array([int(r['tau']) for r in rows]);y=np.array([number(r,'normalized_L') for r in rows])
            se=np.array([number(r,'normalized_naive_se') for r in rows])
            ax.plot(x,y,color=color,lw=1.5,label=r'$N_v=%d$'%n)
            ax.errorbar(x[::4],y[::4],yerr=se[::4],fmt='o',ms=3,color=color,capsize=2,lw=.8)
            bad=y<=0
            if bad.any():ax.scatter(x[bad],y[bad],marker='x',s=25,color=color,zorder=5)
        if R in (24,32):ax.set_yscale('symlog',linthresh=1e-5);ax.axhline(0,color='.3',lw=.7)
        else:ax.set_yscale('log')
        ax.axvline(16,color='.5',ls='--',lw=1);ax.set_xlim(0,32)
        ax.set_title(r'$R/a=%d$'%R);ax.set_xlabel(r'$τ/a$');ax.set_ylabel(r'$⟨L(τ)⟩/⟨L(0)⟩$')
    legend(fig,axes);save(fig,'normalized_correlators_71configs')


def effective_energies():
    data=read('effective_energies_71configs.csv')
    limits={4:(.10,.17),8:(.135,.22),12:(.14,.28),16:(.14,.36),24:(-.5,1.),32:(-1.,1.5)}
    fig,axes=panel_figure('qX: effective energies across four mode cutoffs',
        '71 configurations · piter=10 · logarithms formed after configuration averaging',
        'Paired naive errors every two times. Hollow: endpoint mean/SE <3; triangles: values beyond the display range.\n'
        'Gaps: a nonpositive endpoint mean. Every numerical value and flag is retained in the CSV; these are not plateau fits.')
    for ax,R in zip(axes,SHOW_R):
        ax.set_ylim(*limits[R])
        for n,color in zip(CUTOFFS,COLORS):
            rows=sorted([r for r in data if int(r['Nv'])==n and int(r['R'])==R],key=lambda r:int(r['tau']))
            x=[int(r['tau'])+.5 for r in rows]
            capped_points(ax,x,[number(r,'aE_eff') for r in rows],[number(r,'naive_se') for r in rows],
                color,strong=[r['diagnostic']=='positive_endpoint_means' for r in rows],every=2,label=r'$N_v=%d$'%n)
        ax.axvline(16,color='.5',ls='--',lw=1);ax.set_xlim(0,32);ax.set_title(r'$R/a=%d$'%R)
        ax.set_xlabel(r'$(τ+1/2)/a$');ax.set_ylabel(r'$aE_{eff}=log[⟨L(τ)⟩/⟨L(τ+1)⟩]$')
    legend(fig,axes);save(fig,'effective_energies_71configs')


def signed_support():
    data=read('correlator_means_71configs.csv')
    fig,axes=plt.subplots(2,2,figsize=(12,8.5));axes=axes.ravel()
    for ax,n in zip(axes,CUTOFFS):
        selected={(int(r['R']),int(r['tau'])):r for r in data if int(r['Nv'])==n}
        matrix=np.array([[number(selected[R,t],'mean_over_se') for t in range(33)] for R in RADII])
        ax.grid(False);im=ax.imshow(matrix,origin='lower',aspect='auto',cmap='RdBu_r',vmin=-10,vmax=10,
            extent=(-.5,32.5,-.5,11.5),interpolation='nearest')
        bad=np.where(matrix<=0);ax.scatter(bad[1],bad[0],marker='x',s=18,color='black',lw=.8)
        ax.axvline(16.5,color='black',ls='--',lw=1);ax.set_xticks(range(0,33,4))
        ax.set_yticks(range(12));ax.set_yticklabels(RADII)
        value=matrix[RADII.index(24),24]
        ax.text(24,RADII.index(24),'%.1f'%value,ha='center',va='center',fontsize=8,
            color='white' if value>7 else 'black',fontweight='bold')
        ax.set_title(r'$N_v=%d$'%n);ax.set_xlabel(r'$τ/a$');ax.set_ylabel(r'selected $R/a$')
    fig.suptitle('qX: signed mean / standard error on all 71 configurations',fontsize=15,y=.98)
    fig.text(.5,.938,'piter=10 · labels highlight R/a=24, τ/a=24 · one statistical sample per configuration',ha='center',fontsize=10)
    fig.subplots_adjust(left=.075,right=.86,bottom=.12,top=.87,hspace=.34,wspace=.25)
    bar=fig.colorbar(im,cax=fig.add_axes([.895,.22,.022,.55]),extend='both')
    bar.set_label('Mean / standard error (clipped at +/-10)')
    fig.text(.5,.025,'Crosses: nonpositive ensemble means. All negative individual measurements remain included.\n'
        'Signal ratios diagnose logarithm reliability; they do not establish a ground-state plateau.',ha='center',fontsize=9)
    save(fig,'signed_signal_support_71configs')


def fit_windows():
    data=read('constant_fits_71configs.csv')
    limits={4:(.115,.1185),8:(.156,.174),12:(.18,.215),16:(.175,.24),24:(.05,.62),32:(0.,.53)}
    fig,axes=panel_figure('qX: stability of correlated effective-energy fits',
        '71 matched configurations · six adjacent effective energies per fit · full time covariance',
        'Naive covariance errors; paired jackknife errors are exported. Hollow: low endpoint signal; crosses: approximate p <0.05.\n'
        'Triangles: estimates outside the display range; gaps: undefined input means. Window scans are exploratory, not automatic plateau selection.')
    for ax,R in zip(axes,SHOW_R):
        ax.set_ylim(*limits[R])
        for i,(n,color) in enumerate(zip(CUTOFFS,COLORS)):
            rows=sorted([r for r in data if int(r['Nv'])==n and int(r['R'])==R],key=lambda r:int(r['tau_min']))
            x=np.array([int(r['tau_min']) for r in rows]);y=np.array([number(r,'estimate') for r in rows])
            capped_points(ax,x,y,[number(r,'naive_se') for r in rows],color,
                strong=[r['low_signal_flag']=='False' for r in rows],offset=(i-1.5)*.13,label=r'$N_v=%d$'%n)
            bad=np.array([bool(r['p_gaussian_approx']) and number(r,'p_gaussian_approx')<.05 for r in rows])
            ax.scatter(x[bad]+(i-1.5)*.13,y[bad],marker='x',s=25,color=color,zorder=5)
        ax.set_xlim(3,27);ax.set_xticks([4,8,12,16,20,24]);ax.set_title(r'$R/a=%d$'%R)
        ax.set_xlabel('First τ/a in six-point fit window');ax.set_ylabel('Fitted aE')
    legend(fig,axes);save(fig,'fit_window_stability_71configs')


def fitted_distances():
    data=read('constant_fits_71configs.csv')
    fig,axes=plt.subplots(1,3,figsize=(14,6.5));axes=axes.ravel()
    for ax,start in zip(axes,[6,14,22]):
        ax.set_ylim(0,.65)
        for i,(n,color) in enumerate(zip(CUTOFFS,COLORS)):
            rows=sorted([r for r in data if int(r['Nv'])==n and int(r['tau_min'])==start],key=lambda r:int(r['R']))
            capped_points(ax,[int(r['R']) for r in rows],[number(r,'estimate') for r in rows],
                [number(r,'naive_se') for r in rows],color,strong=[r['low_signal_flag']=='False' for r in rows],
                offset=(i-1.5)*.10,label=r'$N_v=%d$'%n)
        ax.set_xlim(0,33);ax.set_xticks([1,8,16,24,32]);ax.set_xlabel(r'$R/a$');ax.set_ylabel('Fitted aE')
        ax.set_title('τ/a=%d,...,%d'%(start,start+5))
    fig.suptitle('qX: fitted energies at selected spatial distances',fontsize=16,y=.985)
    fig.text(.5,.912,'71 configurations · piter=10 · correlated six-point fits to effective energies',ha='center',fontsize=10)
    h,l=axes[0].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',bbox_to_anchor=(.5,.89),ncol=4,frameon=False)
    fig.text(.5,.035,'Hollow: low endpoint signal. Gaps: undefined input means; triangles: values outside the display range.\n'
        'The time-window and cutoff dependence remain part of the analysis; these curves are not a finalized ground-state potential.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.065,right=.985,bottom=.18,top=.77,wspace=.27)
    save(fig,'fitted_energies_vs_distance_71configs')


def paired_cutoffs():
    data=read('cutoff_fit_contrasts_71configs.csv')
    fig,axes=plt.subplots(2,3,figsize=(14,8.5))
    for column,start in enumerate([6,14,22]):
        for row,(lo,hi) in enumerate([(1,20),(23,32)]):
            ax=axes[row,column]
            if row==0:
                selected=[r for r in data if int(r['tau_min'])==start and lo<=int(r['R'])<=hi and r['difference']]
                lower=min([0.]+[number(r,'difference')-1.3*number(r,'paired_naive_se') for r in selected])
                upper=max([0.]+[number(r,'difference')+1.3*number(r,'paired_naive_se') for r in selected])
                margin=.08*(upper-lower);ax.set_ylim(lower-margin,upper+margin)
            else:ax.set_ylim(-.4,.4)
            for i,(n,color) in enumerate(zip(CUTOFFS[:-1],COLORS[:-1])):
                rows=sorted([r for r in data if int(r['Nv'])==n and int(r['tau_min'])==start and lo<=int(r['R'])<=hi],key=lambda r:int(r['R']))
                capped_points(ax,[int(r['R']) for r in rows],[number(r,'difference') for r in rows],
                    [number(r,'paired_naive_se') for r in rows],color,strong=[r['low_signal_flag']=='False' for r in rows],
                    offset=(i-1)*.10,label=r'$N_v=%d$ vs 100'%n)
            ax.axhline(0,color='.35',lw=1);ax.set_xlim(lo-1,hi+1)
            ax.set_xlabel(r'$R/a$');ax.set_ylabel('Fitted aE(Nv) - aE(100)')
            ax.set_title('τ/a=%d,...,%d; R/a=%d...%d'%(start,start+5,lo,hi))
    fig.suptitle('qX: paired cutoff differences resolve overlapping energy curves',fontsize=15,y=.98)
    fig.text(.5,.936,'71 configurations · covariance between the cutoffs retained in every difference',ha='center',fontsize=10)
    legend(fig,axes.ravel())
    fig.text(.5,.02,'Paired naive errors. Hollow: low endpoint signal; triangles: values outside display limits; gaps: undefined input fits.\n'
        'Small differences can be resolved because the measurements share configurations. No independent-error addition is used.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.075,right=.98,bottom=.12,top=.85,hspace=.35,wspace=.3)
    save(fig,'paired_cutoff_differences_71configs')


def cornell_windows():
    data=read('cornell_diagnostics_71configs.csv')
    fig,axes=plt.subplots(2,2,figsize=(12,8.5))
    for column,tag in enumerate(['4,8,12,16,20','8,12,16,20,24,28']):
        ax=axes[0,column];ax.set_ylim(.0025,.0085)
        paxis=axes[1,column];paxis.set_yscale('log');paxis.set_ylim(1e-12,1.4)
        for i,(n,color) in enumerate(zip(CUTOFFS,COLORS)):
            rows=sorted([r for r in data if int(r['Nv'])==n and r['R_values']==tag],key=lambda r:int(r['tau_min']))
            x=np.array([int(r['tau_min']) for r in rows])
            capped_points(ax,x,[number(r,'sigma_a2') for r in rows],[number(r,'sigma_naive_se') for r in rows],
                color,strong=[r['low_signal_flag']=='False' for r in rows],offset=(i-1.5)*.13,label=r'$N_v=%d$'%n)
            p=np.array([number(r,'p_gaussian_approx') for r in rows]);paxis.plot(x,np.maximum(p,1e-12),'o-',ms=3,color=color,lw=1)
            below=np.isfinite(p)&(p<1e-12);paxis.scatter(x[below],np.full(np.sum(below),1.2e-12),marker='v',s=25,color=color)
        ax.set_title('Selected R/a = '+tag);ax.set_ylabel(r'Cornell slope $σ a^2$')
        paxis.axhline(.05,color='.4',ls='--',lw=1);paxis.set_ylabel('Approximate Gaussian p value')
        for target in (ax,paxis):target.set_xlim(3,27);target.set_xticks([4,8,12,16,20,24]);target.set_xlabel('First τ/a in six-point fit window')
    fig.suptitle('qX: Cornell parameters remain window and range diagnostics',fontsize=15,y=.98)
    fig.text(.5,.936,'aE(R) = V0 - α/R + (σa²)R · full covariance across distances · 71 configurations',ha='center',fontsize=10)
    legend(fig,axes.ravel())
    fig.text(.5,.025,'Naive covariance errors. Hollow slopes: low endpoint signal; triangles: slope outside the display range or p below 10⁻¹².\n'
        'P values are approximate and unadjusted across exploratory scans. Sparse ranges and time dependence prevent a final string-tension determination.',ha='center',fontsize=9)
    fig.subplots_adjust(left=.085,right=.98,bottom=.12,top=.85,hspace=.34,wspace=.3)
    save(fig,'cornell_window_diagnostics_71configs')


if __name__=='__main__':
    normalized_signals();effective_energies();signed_support();fit_windows()
    fitted_distances();paired_cutoffs();cornell_windows()
