#!/usr/bin/env python3
"""Plot the independent-configuration baseline and paired temporal contrasts."""


def render(profiles,temporal,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    output=Path(output)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.hashsalt':'qx-independent-full71'})
    colors=('#3269a8','#df8422','#158471')
    def selected(r):return tuple(int(r[k]) for k in ('axis','R','long','tr2'))==(-1,4,2,0)
    def save(fig,name):
        fig.savefig(output/(name+'.png'),dpi=200,facecolor='white')
        fig.savefig(output/(name+'.svg'),metadata={'Date':None},facecolor='white');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7.6),sharex=True)
    for col,kind in enumerate(('RAW','TL')):
        for tau,color in zip((2,4,6),colors):
            rows=sorted([r for r in profiles if selected(r) and int(r['tau'])==tau],key=lambda r:int(r['tr1']))
            x=[int(r['tr1']) for r in rows];field='C_'+kind+'_p40'
            axes[0,col].errorbar(x,[1e3*float(r[field]) for r in rows],
                yerr=[1e3*float(r[field+'_naive_se']) for r in rows],marker='o',ms=4,capsize=3,
                color=color,label=r'$\tau/a=%d$'%tau)
            axes[1,col].plot(x,[1e4*(float(r[field+'_naive_over_jackknife'])-1) for r in rows],
                marker='o',ms=4,color=color,label=r'$\tau/a=%d$'%tau)
        axes[0,col].set_title('Raw clover' if kind=='RAW' else 'Traceless clover')
        axes[0,col].legend();axes[1,col].axhline(0,color='.4',ls='--',lw=.8)
        axes[1,col].set_xlabel(r'Transverse coordinate $r_\perp/a$')
        for ax in axes[:,col]:ax.grid(alpha=.2);ax.set_xticks((-2,-1,0,1,2))
    axes[0,0].set_ylabel(r'$10^3 C$ with independent errors')
    axes[1,0].set_ylabel(r'$10^4\,(\sigma_{\rm naive}/\sigma_{\rm JK}-1)$')
    fig.suptitle('Independent-configuration errors on the complete qX sample',y=.98,fontsize=14)
    fig.text(.5,.928,r'71 configurations, $piter=40$, $R/a=4$, longitudinal midpoint; combined spatial axes',ha='center')
    fig.text(.5,.025,'Full covariance of the ratio of ensemble means; source sites and axes are averaged within each configuration.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.87,bottom=.11,hspace=.17,wspace=.17,left=.10,right=.98)
    save(fig,'independent_errors_profiles_71configs')
    fig,axes=plt.subplots(1,2,figsize=(11,4.9),sharey=True)
    for ax,kind in zip(axes,('RAW','TL')):
        rows=sorted([r for r in temporal if selected(r) and r['definition']==kind and int(r['piter'])==40],key=lambda r:int(r['tr1']))
        x=[int(r['tr1']) for r in rows]
        for (high,low),color,marker in zip(((4,2),(6,4)),colors[1:],('s','o')):
            field='delta_tau%d_minus_tau%d'%(high,low)
            ax.errorbar(x,[1e5*float(r[field]) for r in rows],yerr=[1e5*float(r[field+'_naive_se']) for r in rows],
                color=color,marker=marker,ms=4,capsize=3,label=r'$C_%d-C_%d$'%(high,low))
        ax.axhline(0,color='.4',ls='--',lw=.8);ax.grid(alpha=.2);ax.legend()
        ax.set_title('Raw clover' if kind=='RAW' else 'Traceless clover')
        ax.set_xlabel(r'Transverse coordinate $r_\perp/a$');ax.set_xticks((-2,-1,0,1,2))
    axes[0].set_ylabel(r'$10^5\,(C_{\tau_{\rm high}}-C_{\tau_{\rm low}})$')
    fig.suptitle('Paired temporal differences with independent-configuration errors',y=.98,fontsize=14)
    fig.text(.5,.88,r'71 configurations, $piter=40$, $R/a=4$; covariance between temporal separations retained',ha='center')
    fig.text(.5,.025,'Three temporal separations provide a consistency check; a ground-state plateau has not been established.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.77,bottom=.17,left=.09,right=.98,wspace=.13)
    save(fig,'temporal_differences_independent_errors_71configs')


if __name__=='__main__':
    import reproduce as a
    render(a.read_csv(a.RESULTS/a.EXPORTS[0]),a.read_csv(a.RESULTS/a.EXPORTS[1]),a.RESULTS)
