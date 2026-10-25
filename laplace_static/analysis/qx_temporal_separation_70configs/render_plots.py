#!/usr/bin/env python3
"""Render the matched qX temporal comparisons from reproducible numerical exports."""


def render(records, covariance_rows, controls, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from pathlib import Path

    output = Path(output)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 12, 'axes.labelsize': 11,
                         'legend.fontsize': 9, 'svg.hashsalt': 'qx-temporal-70configs',
                         'axes.spines.top': False, 'axes.spines.right': False})
    colors = ('#3269a8', '#db7c21', '#168774')
    titles = {'RAW': 'Raw clover', 'TL': 'Traceless clover'}

    def select(row):
        return (int(row['axis']), int(row['R']), int(row['long']), int(row['tr2'])) == (-1, 4, 2, 0)

    def finish(fig, name):
        fig.savefig(output / (name + '.png'), dpi=210, facecolor='white')
        fig.savefig(output / (name + '.svg'), metadata={'Date': None}, facecolor='white')
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True, sharey='row')
    for column, definition in enumerate(('RAW', 'TL')):
        rows = sorted([r for r in records if select(r) and r['definition'] == definition
                       and int(r['piter']) == 40], key=lambda r: int(r['tr1']))
        x = [int(r['tr1']) for r in rows]
        for t, color in zip((2, 4, 6), colors):
            axes[0, column].errorbar(x, [1e3 * float(r['C_tau%d' % t]) for r in rows],
                yerr=[1e3 * float(r['C_tau%d_jackknife_se' % t]) for r in rows],
                color=color, marker='o', markersize=4, capsize=3, lw=1.4,
                label=r'$\tau/a=%d$' % t)
        for (high, low), color, marker in zip(((4, 2), (6, 4)), colors[1:], ('s', 'o')):
            field = 'delta_tau%d_minus_tau%d' % (high, low)
            axes[1, column].errorbar(x, [1e5 * float(r[field]) for r in rows],
                yerr=[1e5 * float(r[field + '_jackknife_se']) for r in rows],
                color=color, marker=marker, markersize=4, capsize=3, lw=1.4,
                label=r'$C_{%d}-C_{%d}$' % (high, low))
        axes[0, column].set_title(titles[definition])
        axes[0, column].legend(loc='lower right')
        axes[1, column].legend(loc='upper right')
        axes[1, column].axhline(0, color='#555555', lw=.8, ls='--')
        for row in (0, 1):
            axes[row, column].grid(alpha=.18)
            axes[row, column].set_xticks((-2, -1, 0, 1, 2))
        axes[1, column].set_xlabel(r'Transverse coordinate $r_\perp/a$')
    axes[0, 0].set_ylabel(r'$10^3 C$')
    axes[1, 0].set_ylabel(r'$10^5\,(C_{\tau_\mathrm{high}}-C_{\tau_\mathrm{low}})$')
    fig.suptitle('Temporal-separation comparison on 70 configurations', y=.975, fontsize=15)
    fig.text(.5, .928, r'$piter=40$, $R/a=4$, longitudinal midpoint; combined spatial axes',
             ha='center', fontsize=10)
    fig.text(.5, .025, 'Paired configuration jackknife errors; preliminary, without autocorrelation correction.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(top=.875, bottom=.10, hspace=.15, wspace=.12, left=.09, right=.98)
    finish(fig, 'temporal_separation_profiles_70configs')

    fig, axes = plt.subplots(1, 2, figsize=(10.3, 5.1))
    for ax, definition in zip(axes, ('RAW', 'TL')):
        rows = [r for r in covariance_rows if select(r) and int(r['tr1']) == 0
                and r['definition'] == definition and int(r['piter']) == 40]
        lookup = {(int(r['tau_i']), int(r['tau_j'])): float(r['correlation']) for r in rows}
        matrix = np.array([[lookup[t, u] for u in (2, 4, 6)] for t in (2, 4, 6)])
        image = ax.imshow(matrix, vmin=-1, vmax=1, cmap='RdBu_r')
        for i in range(3):
            for j in range(3):
                ax.text(j, i, '%.3f' % matrix[i, j], ha='center', va='center', fontsize=13,
                        color='white' if abs(matrix[i, j]) > .65 else '#202020')
        ax.set_xticks((0, 1, 2))
        ax.set_yticks((0, 1, 2))
        ax.set_xticklabels(('2', '4', '6'))
        ax.set_yticklabels(('2', '4', '6'))
        ax.set_xlabel(r'$\tau/a$')
        ax.set_ylabel(r'$\tau/a$')
        ax.set_title(titles[definition])
    fig.suptitle('Correlation between temporal separations', y=.97, fontsize=15)
    fig.text(.5, .89, r'70 configurations, $piter=40$, $R/a=4$, central probe; combined spatial axes',
             ha='center', fontsize=10)
    fig.text(.5, .03, 'Jackknife correlations of the ratio estimates; no autocorrelation correction.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(top=.81, bottom=.16, left=.07, right=.84, wspace=.32)
    cax = fig.add_axes((.89, .20, .018, .56))
    fig.colorbar(image, cax=cax, label='Correlation')
    finish(fig, 'temporal_separation_correlations_70configs')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
    for ax, definition in zip(axes, ('RAW', 'TL')):
        for (high, low), color, marker in zip(((4, 2), (6, 4)), colors[1:], ('s', 'o')):
            rows = sorted([r for r in controls if select(r) and r['definition'] == definition
                           and int(r['tau_high']) == high and int(r['tau_low']) == low],
                          key=lambda r: int(r['tr1']))
            ax.errorbar([int(r['tr1']) for r in rows],
                [1e9 * float(r['contrast_p4_minus_p40']) for r in rows],
                yerr=[1e9 * float(r['contrast_p4_minus_p40_jackknife_se']) for r in rows],
                marker=marker, color=color, markersize=4, capsize=3, lw=1.4,
                label=r'$\Delta=C_{%d}-C_{%d}$' % (high, low))
        ax.set_title(titles[definition])
        ax.axhline(0, color='#555555', lw=.8, ls='--')
        ax.set_xlabel(r'Transverse coordinate $r_\perp/a$')
        ax.set_xticks((-2, -1, 0, 1, 2))
        ax.grid(alpha=.18)
        ax.legend(loc='upper right')
    axes[0].set_ylabel(r'$10^9\,[\Delta(piter=4)-\Delta(piter=40)]$')
    fig.suptitle('Projection control for the temporal differences', y=.98, fontsize=15)
    fig.text(.5, .89, r'70 configurations, $R/a=4$, longitudinal midpoint; combined spatial axes',
             ha='center', fontsize=10)
    fig.text(.5, .025, 'Paired across both temporal separation and projection prescription; preliminary jackknife errors.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(top=.79, bottom=.16, left=.09, right=.98, wspace=.12)
    finish(fig, 'temporal_separation_projection_control_70configs')
