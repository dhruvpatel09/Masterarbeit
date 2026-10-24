#!/usr/bin/env python3
"""Render the matched clover comparison from preserved numerical exports."""
import csv
from pathlib import Path

PLOTS = ['clover_raw_vs_traceless_70configs', 'clover_trace_fraction_70configs',
         'profiles_p4_vs_p40_traceless_70configs']
COLORS = {2:'#176b8e',4:'#278260',6:'#a05a9c'}
FOOTER = '70 matched configurations; preliminary paired jackknife errors; no autocorrelation correction'


def cut(rows, tau, r=4, longitudinal=2):
    return sorted([row for row in rows if tuple(int(row[k]) for k in
                  ('tau','axis','R','long','tr2')) == (tau,-1,r,longitudinal,0)],
                  key=lambda row:int(row['tr1']))


def render(rows, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.hashsalt':'qx-clover-trace-full70','axes.titlepad':12})

    def values(selected, name, factor=1):
        return [factor*float(row[name]) for row in selected]

    def style(ax):
        ax.grid(alpha=.18)
        ax.set_xticks(range(-2,3))

    def save(fig, name):
        for ext in ('png','svg'):
            metadata = {'Date':None} if ext == 'svg' else {'Software':'Matplotlib'}
            fig.savefig(str(output/(name+'.'+ext)),dpi=190,metadata=metadata)
        plt.close(fig)

    fig, axes = plt.subplots(2,3,figsize=(13.0,7.5),sharex=True,sharey='row')
    for col,tau in enumerate((2,4,6)):
        selected = cut(rows,tau)
        x = values(selected,'tr1')
        color = COLORS[tau]
        for field,marker,linestyle,label in [('C_RAW_p40','s','--','Raw clover'),
                                             ('C_TL_p40','o','-','Traceless clover')]:
            axes[0,col].errorbar(x,values(selected,field,1e3),
                  yerr=values(selected,field+'_jackknife_se',1e3),
                  marker=marker,linestyle=linestyle,color=color,alpha=.85 if marker=='s' else 1,
                  markerfacecolor='white' if marker=='s' else color,capsize=3,label=label)
        axes[0,col].set_title(r'$\tau/a=%d$' % tau)
        axes[1,col].errorbar(x,values(selected,'C_TRACE_p40',1e6),
                  yerr=values(selected,'C_TRACE_p40_jackknife_se',1e6),
                  fmt='o-',color=color,capsize=3)
        axes[1,col].set_xlabel(r'$r_\perp/a$')
        for ax in axes[:,col]:
            style(ax)
    axes[0,0].legend(loc='upper left',fontsize=10)
    axes[0,0].set_ylabel(r'Connected profile $C\;(\times10^{-3})$')
    axes[1,0].set_ylabel(r'$C_{\rm RAW}-C_{\rm TL}\;(\times10^{-6})$')
    fig.suptitle('Raw and colour-traceless clover: matched ensemble comparison\n'+
                 r'$piter=40$, $R/a=4$, longitudinal midpoint, original gauge-field probe',fontsize=14,y=.975)
    fig.text(.5,.018,FOOTER,ha='center',fontsize=10,color='#404040')
    fig.subplots_adjust(left=.09,right=.98,top=.82,bottom=.115,wspace=.12,hspace=.24)
    save(fig,PLOTS[0])

    fig,axes = plt.subplots(1,3,figsize=(12.5,4.9),sharey=True)
    for col,(r,longitudinal) in enumerate(((1,0),(2,1),(4,2))):
        for tau in (2,4,6):
            selected = cut(rows,tau,r,longitudinal)
            axes[col].errorbar(values(selected,'tr1'),values(selected,'trace_fraction_percent_p40'),
                    yerr=values(selected,'trace_fraction_percent_p40_jackknife_se'),
                    fmt='o-',color=COLORS[tau],capsize=3,label=r'$\tau/a=%d$' % tau)
        axes[col].set_title(r'$R/a=%d$, longitudinal offset $l/a=%d$' % (r,longitudinal),fontsize=11)
        axes[col].set_xlabel(r'$r_\perp/a$')
        style(axes[col])
    axes[0].set_ylabel(r'$100\,(C_{\rm RAW}-C_{\rm TL})/C_{\rm RAW}$ [%]')
    axes[2].legend(fontsize=10)
    fig.suptitle('Relative change from colour-trace subtraction\n'+
                 r'$piter=40$; ratios and errors retain the covariance of both definitions',fontsize=14,y=.975)
    fig.text(.5,.018,FOOTER,ha='center',fontsize=10,color='#404040')
    fig.subplots_adjust(left=.08,right=.98,top=.73,bottom=.17,wspace=.12)
    save(fig,PLOTS[1])

    fig,axes = plt.subplots(2,1,figsize=(8.4,7.8),sharex=True)
    for tau in (2,4,6):
        selected = cut(rows,tau)
        x = values(selected,'tr1')
        y = values(selected,'C_TL_p40',1e3)
        errors = values(selected,'C_TL_p40_jackknife_se',1e3)
        axes[0].plot(x,y,'o-',color=COLORS[tau],label=r'$\tau/a=%d$' % tau)
        axes[0].fill_between(x,[a-b for a,b in zip(y,errors)],[a+b for a,b in zip(y,errors)],
                             color=COLORS[tau],alpha=.12)
        axes[0].plot(x,values(selected,'C_TL_p4',1e3),'x--',color=COLORS[tau])
        axes[1].errorbar(x,values(selected,'delta_C_TL_p4_minus_p40',1e9),
                 yerr=values(selected,'delta_C_TL_p4_minus_p40_jackknife_se',1e9),
                 fmt='o-',color=COLORS[tau],capsize=3)
    axes[0].legend(loc='upper left',fontsize=10)
    handles=[Line2D([0],[0],color='#303030',marker='o',label='piter=40'),
             Line2D([0],[0],color='#303030',marker='x',linestyle='--',label='piter=4')]
    axes[0].add_artist(axes[0].get_legend())
    axes[0].legend(handles=handles,loc='upper right',fontsize=10)
    axes[0].set_ylabel(r'$C_{\rm TL}\;(\times10^{-3})$')
    axes[1].set_ylabel(r'$C_{{\rm TL},4}-C_{{\rm TL},40}\;(\times10^{-9})$')
    axes[1].set_xlabel(r'$r_\perp/a$')
    axes[1].axhline(0,color='#555555',linewidth=.8)
    for ax in axes:
        style(ax)
    fig.suptitle('Projection sensitivity with the traceless clover\n'+
                 r'$R/a=4$, longitudinal midpoint; paired differences below',fontsize=14,y=.975)
    fig.text(.5,.02,'70 matched configurations; preliminary jackknife errors\n'+
             'Shading: p40 profile error; bars: paired p4-p40 difference error',
             ha='center',fontsize=10,color='#404040')
    fig.subplots_adjust(left=.16,right=.97,top=.85,bottom=.13,hspace=.23)
    save(fig,PLOTS[2])


if __name__ == '__main__':
    here = Path(__file__).resolve().parent
    with (here/'results/profiles_ensemble_70configs.csv').open(newline='') as stream:
        render(list(csv.DictReader(stream)),here/'results')
