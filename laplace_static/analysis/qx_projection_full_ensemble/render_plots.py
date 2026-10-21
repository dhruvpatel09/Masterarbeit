#!/usr/bin/env python3
"""Render full70 ensemble figures without changing the preserved estimator."""

def plot_results(records, windows, output, analysis):
    CONFIGS, PITERS, PLOTS = analysis.CONFIGS, analysis.PITERS, analysis.PLOTS
    COLORS, COHORT = analysis.COLORS, analysis.COHORT
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
        'svg.hashsalt': 'qx-projection-full70'})

    def save(fig, name):
        fig.savefig(str(output/(name+'.png')), dpi=180)
        fig.savefig(str(output/(name+'.svg')), metadata={'Date': None})
        path = output/(name+'.svg')
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        plt.close(fig)

    def cut(tau):
        return sorted([r for r in records if (
            r['tau'], r['axis'], r['R'], r['long'], r['tr2']) == (tau, -1, 4, 2, 0)],
            key=lambda r: r['tr1'])

    def decorate(ax):
        ax.grid(alpha=.18)
        ax.set_xticks([-2, -1, 0, 1, 2])

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, layout='constrained')
    for column, tau in enumerate((2, 4, 6)):
        rows = cut(tau)
        x = [r['tr1'] for r in rows]
        for p in PITERS:
            axes[0, column].plot(x, [float(r['C_p'+str(p)]) for r in rows],
                                 'o-', ms=4, color=COLORS[p], label='piter='+str(p))
            if p != 40:
                axes[1, column].plot(x, [float(r['delta_C_p'+str(p)+'_minus_p40']) for r in rows],
                                     'o-', ms=4, color=COLORS[p], label='piter='+str(p))
        axes[0, column].set_title(r'$\tau/a=%d$' % tau)
        axes[0, column].legend(fontsize=9)
        axes[1, column].legend(fontsize=9)
        axes[1, column].set_xlabel(r'$r_\perp/a$')
        axes[1, column].axhline(0, color='#87919b', lw=.7)
        for ax in axes[:, column]:
            decorate(ax)
    axes[0, 0].set_ylabel(r'$C_p$')
    axes[1, 0].set_ylabel(r'$C_p-C_{40}$')
    fig.suptitle(COHORT+'\n'+r'Four measured prescriptions; $R/a=4$, longitudinal midpoint', fontsize=14)
    fig.supxlabel('Raw clover convention; no final statistical error bars.', fontsize=10)
    save(fig, PLOTS[0])

    tau_colors = {2: '#176b8e', 4: '#cd5519', 6: '#437c45'}
    for errors, name in ((False, PLOTS[2]), (True, PLOTS[3])):
        fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.3), sharex=True, layout='constrained')
        for tau in (2, 4, 6):
            rows = cut(tau)
            x = [r['tr1'] for r in rows]
            color = tau_colors[tau]
            if errors:
                axes[0].errorbar(x, [float(r['C_p40']) for r in rows],
                                 yerr=[float(r['C_p40_jackknife_se']) for r in rows],
                                 fmt='o-', color=color, capsize=4, label=r'$\tau/a=%d$' % tau)
                axes[1].errorbar(x, [float(r['delta_C_p4_minus_p40']) for r in rows],
                                 yerr=[float(r['delta_C_paired_jackknife_se']) for r in rows],
                                 fmt='o-', color=color, capsize=4, label=r'$\tau/a=%d$' % tau)
            else:
                for p, style in ((4, '-'), (40, '--')):
                    axes[0].plot(x, [float(r['C_p'+str(p)]) for r in rows],
                                 style, marker='o', ms=4, color=color,
                                 label=r'$\tau/a=%d$; piter=%d' % (tau, p))
                axes[1].plot(x, [float(r['delta_C_p4_minus_p40']) for r in rows],
                             'o-', color=color, label=r'$\tau/a=%d$' % tau)
        for ax in axes:
            decorate(ax)
            ax.axhline(0, color='#87919b', lw=.7)
            ax.legend(ncol=3, fontsize=9)
        axes[0].set_ylabel(r'$C_{40}$ with preliminary SE' if errors else r'$C_p$')
        axes[1].set_ylabel(r'$C_4-C_{40}$'+(' with paired SE' if errors else ''))
        axes[1].set_xlabel(r'$r_\perp/a$')
        fig.suptitle(COHORT+'\n'+r'Raw profiles, $R/a=4$, longitudinal midpoint', fontsize=14)
        fig.supxlabel('Matched delete-one-configuration jackknife; no autocorrelation correction; preliminary.'
                      if errors else 'Paired piter=4/40 comparison; raw clover convention.', fontsize=10)
        save(fig, name)

    def window_axes():
        fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), layout='constrained')
        for ax in axes:
            ax.set_yscale('log')
            ax.grid(axis='y', alpha=.2)
            ax.set_xticks(CONFIGS[::4]+[CONFIGS[-1]])
            ax.tick_params(axis='x', labelrotation=45, labelsize=8)
            ax.axvspan(635, 645, color='#87919b', alpha=.12)
            ax.set_xlabel('Configuration n (n640 unavailable)')
        return fig, axes

    def window_curve(ax, values, label, color, style='-', marker='o'):
        for start, stop in ((0, 54), (54, len(CONFIGS))):
            ax.plot(CONFIGS[start:stop], values[start:stop], style, marker=marker,
                    ms=4, color=color, label=label if start == 0 else '_nolegend_')

    fig, axes = window_axes()
    for ax, field, title in zip(axes, ('unit_max', 'det_max'),
                               ('Frobenius unitarity defect', 'Determinant defect')):
        for p in PITERS:
            values = [max(float(r[field]) for r in windows
                          if int(r['configuration']) == cfg and int(r['piter']) == p)
                      for cfg in CONFIGS]
            window_curve(ax, values, 'piter='+str(p), COLORS[p],
                         style='--' if p in (10, 20) else '-')
        ax.set_title(title)
        ax.set_ylabel('Largest maximum over the three central windows')
        ax.legend(fontsize=9)
    fig.suptitle(COHORT+'\n'+r'One temporal HYP step; midpoint $t/a=98$; $\tau/a=2,4,6$', fontsize=14)
    fig.supxlabel('Measured windows only; Frobenius norm differs from the n100 entry-norm diagnostic.', fontsize=10)
    save(fig, PLOTS[1])

    fig, axes = window_axes()
    for ax, maximum, rms, title in zip(axes, ('unit_max', 'det_max'), ('unit_rms', 'det_rms'),
                                      ('Frobenius unitarity defect', 'Determinant defect')):
        for p, field, style, color, marker in (
                (4, maximum, '-', COLORS[4], 'o'),
                (4, rms, '--', '#437c45', 's'),
                (40, maximum, ':', COLORS[40], '^')):
            values = [max(float(r[field]) for r in windows
                          if int(r['configuration']) == cfg and int(r['piter']) == p)
                      for cfg in CONFIGS]
            label = 'piter=%d: largest window %s' % (p, 'RMS' if field == rms else 'maximum')
            window_curve(ax, values, label, color, style, marker)
        ax.axhline(.01, color='#87919b', ls='--', lw=.8)
        ax.set_ylim(1e-16, .05)
        ax.set_title(title)
        ax.legend(fontsize=9, loc='center right', bbox_to_anchor=(1, .3))
    fig.suptitle(COHORT+'\nCentral transporter windows: maxima and RMS', fontsize=14)
    fig.supxlabel('0.01 is a diagnostic counter threshold, not an agreed physics acceptance limit.', fontsize=10)
    save(fig, PLOTS[4])
