#!/usr/bin/env python3
"""Render conditional binning diagnostics with explicit cohort support."""


def render(records, plans, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    from pathlib import Path

    output = Path(output)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 12, 'axes.labelsize': 10, 'legend.fontsize': 9,
        'svg.hashsalt': 'qx-error-binning-70configs',
        'axes.spines.top': False, 'axes.spines.right': False})
    colors = ('#3269a8', '#db7c21', '#168774')
    titles = {'RAW': 'Raw clover', 'TL': 'Traceless clover'}

    def save(fig, name):
        fig.savefig(output / (name + '.png'), dpi=210, facecolor='white')
        fig.savefig(output / (name + '.svg'), metadata={'Date': None}, facecolor='white')
        plt.close(fig)

    def midpoint(r):
        return tuple(int(r[k]) for k in ('axis', 'R', 'long', 'tr1', 'tr2')) == (-1, 4, 2, 0, 0)

    panels = [('C_tau%d' % t, r'$C_%d$' % t, color) for t, color in zip((2, 4, 6), colors)]
    differences = [('delta_tau4_minus_tau2', r'$C_4-C_2$', colors[1]),
                   ('delta_tau6_minus_tau4', r'$C_6-C_4$', colors[2])]
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.5), sharex=True, sharey='row')
    for column, definition in enumerate(('RAW', 'TL')):
        selected = sorted([r for r in records if midpoint(r) and r['definition'] == definition],
                          key=lambda r: int(r['bin_width']))
        widths = [int(r['bin_width']) for r in selected]
        for row, fields in enumerate((panels, differences)):
            ax = axes[row, column]
            for field, label, color in fields:
                ax.plot(widths, [float(r[field + '_se_ratio']) for r in selected],
                        marker='o', color=color, lw=1.5, markersize=4, label=label)
            ax.axhline(1, color='#555555', lw=.9, ls='--')
            ax.grid(alpha=.18)
            ax.legend(loc='best')
            ax.set_xscale('log', base=2)
            ax.set_xticks((1, 2, 4, 8))
            ax.set_xticklabels(('1\nN=70', '2\nN=70', '4\nN=68', '8\nN=64'))
        axes[0, column].set_title(titles[definition])
        axes[1, column].set_xlabel('Configurations per bin; selected configuration count N')
    axes[0, 0].set_ylabel('Profile error ratio\nblocked / unblocked, same cohort')
    axes[1, 0].set_ylabel('Temporal-difference error ratio\nblocked / unblocked, same cohort')
    fig.suptitle('Binning sensitivity at the central probe', y=.98, fontsize=15)
    fig.text(.5, .928, r'Conditional configuration-label ordering; $piter=40$, $R/a=4$, combined axes',
             ha='center', fontsize=10)
    fig.text(.5, .025, 'The n640 gap is preserved. Bin counts: 70, 35, 17, 8. These diagnostics do not establish independence.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.10, right=.98, top=.86, bottom=.16, hspace=.14, wspace=.14)
    save(fig, 'binning_error_ratios_midpoint')

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.5), sharex=True, sharey='row')
    for column, definition in enumerate(('RAW', 'TL')):
        selected = sorted([r for r in records if r['definition'] == definition and int(r['bin_width']) == 2
                           and tuple(int(r[k]) for k in ('axis', 'R', 'long', 'tr2')) == (-1, 4, 2, 0)],
                          key=lambda r: int(r['tr1']))
        x = [int(r['tr1']) for r in selected]
        for row, fields in enumerate((panels, differences)):
            ax = axes[row, column]
            for field, label, color in fields:
                ax.plot(x, [float(r[field + '_se_ratio']) for r in selected],
                        color=color, marker='o', lw=1.5, markersize=4, label=label)
            ax.axhline(1, color='#555555', lw=.9, ls='--')
            ax.grid(alpha=.18)
            ax.legend(loc='best')
            ax.set_xticks((-2, -1, 0, 1, 2))
        axes[0, column].set_title(titles[definition])
        axes[1, column].set_xlabel(r'Transverse coordinate $r_\perp/a$')
    axes[0, 0].set_ylabel('Profile error ratio\nblocked / unblocked')
    axes[1, 0].set_ylabel('Temporal-difference error ratio\nblocked / unblocked')
    fig.suptitle('Two-configuration bins: all 70 samples retained', y=.98, fontsize=15)
    fig.text(.5, .928, r'Conditional configuration-label ordering; $piter=40$, $R/a=4$, longitudinal midpoint',
             ha='center', fontsize=10)
    fig.text(.5, .025, '35 bins; identical full70 point estimates. This is an error-sensitivity check, without an MDU assignment.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.10, right=.98, top=.86, bottom=.10, hspace=.14, wspace=.14)
    save(fig, 'binning_two_configs_full70_profiles')

    nominal = list(range(100, 801, 10))
    support = np.zeros((len(plans), len(nominal)), dtype=int)
    for i, plan in enumerate(plans):
        for j, c in enumerate(nominal):
            if c in plan['selected_configurations']:
                support[i, j] = 1
            elif c in plan['dropped_incomplete_bin_configurations']:
                support[i, j] = 2
    fig, ax = plt.subplots(figsize=(11.3, 3.7))
    palette = ['#b7bdc5', '#168774', '#e59a35']
    ax.imshow(support, aspect='auto', interpolation='nearest', cmap=ListedColormap(palette),
              vmin=0, vmax=2, extent=(95, 805, len(plans) - .5, -.5))
    ax.set_yticks(range(len(plans)))
    ax.set_yticklabels(['b=%d: N=%d, B=%d' % (p['bin_width'], p['configurations'], p['bins']) for p in plans])
    ax.set_xticks((100, 200, 300, 400, 500, 600, 640, 700, 800))
    ax.set_xlabel('Configuration label n; no MDU conversion assigned')
    ax.set_title('Configuration support for each bin width', pad=15)
    ax.legend(handles=[Patch(facecolor=palette[1], label='Included'),
                       Patch(facecolor=palette[2], label='Incomplete end bin omitted'),
                       Patch(facecolor=palette[0], label='n640 absent')],
              loc='upper center', bbox_to_anchor=(.5, -.28), ncol=3, frameon=False)
    fig.text(.5, .02, 'No bin crosses the missing n640 label. Larger-bin subsets are diagnostics; the full70 results remain unchanged.',
             ha='center', fontsize=9, color='#444444')
    fig.subplots_adjust(left=.19, right=.985, top=.79, bottom=.38)
    save(fig, 'binning_configuration_support')
