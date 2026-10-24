#!/usr/bin/env python3
"""Plot central-window accuracy and matched observable convergence."""
import csv
import json
from pathlib import Path

PLOTS=['projection_accuracy_and_stability_4configs','profiles_raw_piter_spectrum_4configs',
       'profiles_traceless_piter_spectrum_4configs','per_configuration_projection_stability']
PITERS=[4,6,8,10,20,40,80]
TL=[4,6,8,40,80]
COLORS={4:'#cf5a18',6:'#c93442',8:'#24845e',10:'#0e8d9c',20:'#4768aa',40:'#222222',80:'#9760a0'}
MARKERS={4:'x',6:'v',8:'s',10:'^',20:'d',40:'o',80:'+'}
COHORT='Matched diagnostic sample: n100, n170, n670, n740'


def render(records,summary,individual,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'qx-clover-spectrum-four-configs'})
    def save(fig,name):
        for ext in ('png','svg'):
            fig.savefig(str(output/(name+'.'+ext)),dpi=190,
                        metadata={'Date':None} if ext=='svg' else {'Software':'Matplotlib'})
        plt.close(fig)
    def val(rows,key,scale=1):return [scale*float(r[key]) for r in rows]
    locations={p:i for i,p in enumerate(PITERS)}

    fig,axes=plt.subplots(2,1,figsize=(9.8,8.1),sharex=True)
    defects=summary['central_window_defects']
    x=[locations[int(r['piter'])] for r in defects]
    axes[0].semilogy(x,val(defects,'unit_max'),'o-',color='#176b8e',label=r'Maximum $\|U^\dagger U-I\|_F$')
    axes[0].semilogy(x,val(defects,'det_max'),'s--',color='#d06727',label=r'Maximum $|\det U-1|$')
    axes[0].set_ylabel('Central-link defect')
    axes[0].legend(fontsize=10)
    for definition,marker,color in [('RAW','o','#176b8e'),('TL','s','#9760a0')]:
        selected=[r for r in summary['comparison'] if r['definition']==definition and int(r['piter'])!=40]
        measured={int(r['piter']):float(r['max_abs_delta_C_vs_p40']) for r in selected}
        axes[1].semilogy(range(len(PITERS)),[measured.get(p,float('nan')) for p in PITERS],
                         marker+'--',color=color,label='Raw clover' if definition=='RAW' else 'Traceless clover')
    axes[1].set_ylabel(r'Maximum $|C_p-C_{40}|$')
    axes[1].set_xlabel('Projection iterations piter')
    axes[1].legend(loc='lower left',fontsize=10)
    axes[1].text(.98,.96,'p40 is the zero-difference reference',transform=axes[1].transAxes,
                 ha='right',va='top',fontsize=10,color='#444444')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.axvline(locations[40],color='#999999',linestyle=':',linewidth=1)
        ax.set_xticks(range(len(PITERS)))
        ax.set_xticklabels(PITERS)
    fig.suptitle('SU(3) accuracy and observable iteration stability\nCentral windows at '+r'$\tau/a=2,4,6$',fontsize=14,y=.97)
    fig.text(.5,.017,COHORT+'\nLower panel: maxima over 249 combined-axis ensemble positions; p40 zero omitted',ha='center',fontsize=10)
    fig.subplots_adjust(left=.12,right=.97,top=.85,bottom=.12,hspace=.22)
    save(fig,PLOTS[0])

    for definition,piters,name in [('RAW',PITERS,PLOTS[1]),('TL',TL,PLOTS[2])]:
        fig,axes=plt.subplots(2,3,figsize=(13.1,7.9),sharex=True,sharey='row')
        for col,tau in enumerate((2,4,6)):
            selected=sorted([r for r in records if tuple(int(r[k]) for k in ('tau','axis','R','long','tr2'))==(tau,-1,4,2,0)],
                            key=lambda r:int(r['tr1']))
            x=val(selected,'tr1')
            center=val(selected,'C_%s_p40'%definition,1e3)
            err=val(selected,'C_%s_p40_jackknife_se'%definition,1e3)
            axes[0,col].fill_between(x,[a-b for a,b in zip(center,err)],[a+b for a,b in zip(center,err)],color='#555555',alpha=.10)
            for p in piters:
                axes[0,col].plot(x,val(selected,'C_%s_p%d'%(definition,p),1e3),
                    marker=MARKERS[p],linestyle='-' if p==40 else '--',color=COLORS[p],
                    markerfacecolor='white' if p not in (4,40,80) else COLORS[p],
                    linewidth=1.8 if p==40 else 1.1,alpha=.9)
                if p!=40:
                    field='delta_C_%s_p%d_minus_p40'%(definition,p)
                    axes[1,col].errorbar(x,val(selected,field),yerr=val(selected,field+'_jackknife_se'),
                        marker=MARKERS[p],linestyle='--',color=COLORS[p],capsize=2,linewidth=1)
            axes[0,col].set_title(r'$\tau/a=%d$'%tau)
            axes[1,col].set_yscale('symlog',linthresh=1e-14)
            axes[1,col].axhline(0,color='#666666',linewidth=.7)
            axes[1,col].set_xlabel(r'$r_\perp/a$')
            for ax in axes[:,col]:
                ax.grid(alpha=.18)
                ax.set_xticks(range(-2,3))
        axes[0,0].set_ylabel(r'$C_{\rm %s}\;(\times10^{-3})$'%definition)
        axes[1,0].set_ylabel(r'$C_p-C_{40}$ (symmetric log scale)')
        handles=[Line2D([0],[0],color=COLORS[p],marker=MARKERS[p],linestyle='-' if p==40 else '--',label='piter=%d'%p) for p in piters]
        fig.legend(handles=handles,ncol=len(piters),loc='upper center',bbox_to_anchor=(.52,.89),fontsize=10,frameon=False)
        fig.suptitle(('Raw' if definition=='RAW' else 'Colour-traceless')+' clover: matched projection spectrum\n'+
                     r'$R/a=4$, longitudinal midpoint; four-configuration diagnostic',fontsize=14,y=.98)
        fig.text(.5,.016,COHORT+'\nPreliminary four-configuration jackknife: p40 shading; paired difference bars\n'+
                 r'No autocorrelation correction; difference scale is linear for $|\Delta C|\leq10^{-14}$',
                 ha='center',fontsize=9.6)
        fig.subplots_adjust(left=.095,right=.98,top=.75,bottom=.15,wspace=.12,hspace=.24)
        save(fig,name)

    fig,axes=plt.subplots(2,1,figsize=(10.1,7.8),sharex=True)
    configs=[100,170,670,740]
    for ax,definition,piters in zip(axes,('RAW','TL'),(PITERS,TL)):
        for p in piters:
            if p==40:continue
            y=[max(float(r['max_absolute_delta_C']) for r in individual if
                   int(r['configuration'])==c and r['definition']==definition and int(r['piter'])==p) for c in configs]
            ax.semilogy(range(4),y,marker=MARKERS[p],linestyle='--',color=COLORS[p],label='piter=%d'%p)
        ax.set_ylabel(('Raw' if definition=='RAW' else 'Traceless')+r': max $|C_p-C_{40}|$')
        ax.grid(alpha=.18)
        ax.legend(ncol=3,fontsize=9,loc='center right',framealpha=.95)
    axes[1].set_xticks(range(4));axes[1].set_xticklabels(['n%d'%c for c in configs])
    axes[1].set_xlabel('Configuration')
    fig.suptitle('Projection stability on each configuration\nMaxima over the three temporal separations and measured combined-axis positions',fontsize=13,y=.97)
    fig.text(.5,.02,COHORT+'\nSingle-configuration diagnostics; these maxima are descriptive, without ensemble cancellation',ha='center',fontsize=10)
    fig.subplots_adjust(left=.13,right=.97,top=.84,bottom=.13,hspace=.22)
    save(fig,PLOTS[3])


if __name__=='__main__':
    here=Path(__file__).resolve().parent
    def read(name):
        with (here/'results'/name).open(newline='') as f:return list(csv.DictReader(f))
    render(read('profiles_ensemble_4configs.csv'),json.loads((here/'results/SUMMARY.json').read_text()),
           read('per_configuration_comparison.csv'),here/'results')
