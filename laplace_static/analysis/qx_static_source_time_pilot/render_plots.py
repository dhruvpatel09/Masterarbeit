"""Descriptive pilot figures. Each configuration remains visible; no error bars."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import HERE, CONFIGS, RADII, POSITIONS, csv_read

COLORS = {58: '#0072B2', 90: '#D55E00', 122: '#7B3294'}
MARKERS = {58: 'o', 90: 's', 122: '^'}
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'svg.fonttype': 'none', 'svg.hashsalt': 'qx-source-time-pilot-6036903'})


def number(row, key):
    return float(row[key]) if row[key] != '' else np.nan


def matching(rows, **conditions):
    return sorted([r for r in rows if all(r[k] == str(v) for k, v in conditions.items())],
                  key=lambda r: int(r.get('tau', 0)))


def positions():
    return [Line2D([], [], color=COLORS[t], marker=MARKERS[t], lw=1.1,
                   label=r'$t_0/a={}$'.format(t)) for t in POSITIONS]


def panel_grid():
    fig, axes = plt.subplots(3, 4, figsize=(15.6, 10), sharex=True)
    for i, cfg in enumerate(CONFIGS):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            ax.set_title(r'n{}, $R/a={}$'.format(cfg, radius))
            ax.grid(alpha=.18)
            ax.set_xlim(-.4, 32.4)
            ax.set_xticks([0, 8, 16, 24, 32])
            if i == 2:
                ax.set_xlabel(r'Correlator separation $\tau/a$')
    return fig, axes


def save(fig, output, name, title, subtitle, footer, handles):
    height = fig.get_size_inches()[1]
    fig.suptitle(title, fontsize=15, y=1 - .12 / height, va='top')
    fig.text(.5, 1 - .48 / height, subtitle, ha='center', va='top', fontsize=10)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, 1 - .77 / height),
               ncol=min(len(handles), 5), frameon=False, fontsize=9)
    fig.text(.5, .10 / height, footer, ha='center', va='bottom', fontsize=9, linespacing=1.3)
    bottom = 1.20 if name == 'frozen_basis_training_controls' else .90
    fig.subplots_adjust(top=1 - 1.30 / height, bottom=bottom / height,
                        left=.075, right=.985, hspace=.42, wspace=.32)
    for extension in ('png', 'svg'):
        fig.savefig(output / (name + '.' + extension), dpi=170, facecolor='white')
    plt.close(fig)


def projected_curves(rows, output, own=False):
    fig, axes = panel_grid()
    for i, cfg in enumerate(CONFIGS):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            for t0 in POSITIONS:
                for basis, line, alpha in [('pair', '--', .65), ('full', '-', 1.)]:
                    a = matching(rows, configuration=cfg, R=radius, basis=basis, anchor=4, t0=t0)
                    field = 'own_equal_time_normalized_re' if own else 'common_t90_normalized_re'
                    ax.plot([int(r['tau']) for r in a], [number(r, field) for r in a],
                            color=COLORS[t0], marker=MARKERS[t0], ls=line, lw=1.05,
                            ms=2.8, markevery=2, alpha=alpha)
            ax.set_yscale('symlog', linthresh=1e-4)
            ax.axhline(0, color='.5', lw=.65)
            if j == 0:
                ax.set_ylabel('Signed normalized Re Q')
    handles = positions() + [Line2D([], [], color='.2', ls='--', label='pair basis'),
                             Line2D([], [], color='.2', ls='-', label='full basis')]
    name = 'own_time_normalized_correlators' if own else 'source_position_signed_correlators'
    normalization = 'Each curve divided by its own-position Re Q(0)' if own else 'Common denominator: baseline Re Q(t0=90, tau=0) for each configuration and basis'
    save(fig, output, name, 'Frozen-source correlators at three temporal positions', normalization,
         'Pair (40,100) and full (40,60,80,100); frozen training at separation 4a on the other 68 configurations.\n'
         'Signed symlog axes retain negative measurements. Three technical cases; no ensemble error or potential fit.', handles)


def phases(rows, output):
    fig, axes = panel_grid()
    for i, cfg in enumerate(CONFIGS):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            for t0 in POSITIONS:
                a = matching(rows, configuration=cfg, R=radius, basis='full', anchor=4, t0=t0)
                x = np.array([int(r['tau']) for r in a])
                y = np.array([number(r, 'phase') for r in a])
                line = y.copy()
                for k in range(1, len(line)):
                    if abs(y[k] - y[k - 1]) > np.pi:
                        line[k] = np.nan
                ax.plot(x, line, color=COLORS[t0], lw=.9)
                ax.plot(x, y, color=COLORS[t0], marker=MARKERS[t0], ms=2.8, ls='none')
            ax.axhline(0, color='.5', lw=.6)
            for sign in (-1, 1):
                ax.axhline(sign * .05, color='.6', ls=':', lw=.7)
            ax.set_ylim(-3.4, 3.4)
            if j == 0:
                ax.set_ylabel('Phase of Q [radians]')
    save(fig, output, 'source_position_phases', 'Complex phase and real-sign changes across positions',
         'Full basis; frozen training at separation 4a; fixed coefficients at every position',
         'At R/a=32 the imaginary part is roundoff: phase 0 means positive real Q; +/-pi are the same negative-real direction.\n'
         '+/-0.05 rad guides are descriptive. No phase standard errors or significance tests from three cases.', positions())


def differences(rows, output):
    fig, axes = panel_grid()
    for i, cfg in enumerate(CONFIGS):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            for t0 in (58, 122):
                a = matching(rows, configuration=cfg, R=radius, basis='full', anchor=4, t0=t0)
                ax.plot([int(r['tau']) for r in a], [number(r, 'common_t90_normalized_delta_re') for r in a],
                        color=COLORS[t0], marker=MARKERS[t0], ms=2.7, lw=.95)
            ax.set_yscale('symlog', linthresh=1e-4)
            ax.axhline(0, color='.45', lw=.75)
            if j == 0:
                ax.set_ylabel('Normalized Re Q shift')
    save(fig, output, 'source_position_paired_changes', 'Within-configuration source-position changes',
         'Full basis; training separation 4a; same frozen source at both positions in every difference',
         'Y = [Re Q(position,tau) - Re Q(90,tau)] / Re Q(90,0). Signed symlog axes; panel scales may differ.\n'
         'Points are paired measurements, not independent time samples; no statistical significance assigned.',
         [Line2D([], [], color=COLORS[t], marker=MARKERS[t], label='%d minus 90' % t) for t in (58, 122)])


def matrices(rows, output):
    fig, axes = panel_grid()
    for i, cfg in enumerate(CONFIGS):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            for t0 in (58, 122):
                a = matching(rows, configuration=cfg, R=radius, t0=t0)
                ax.plot([int(r['tau']) for r in a], [number(r, 'relative_delta_to_baseline') for r in a],
                        color=COLORS[t0], marker=MARKERS[t0], ms=2.7, lw=.95)
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 1e2)
            ax.axhline(1, color='.5', ls=':', lw=.7)
            if j == 0:
                ax.set_ylabel(r'$\|\widetilde C_{t_0}-\widetilde C_{90}\|_F/\|\widetilde C_{90}\|_F$')
    save(fig, output, 'matrix_position_changes', 'Changes of the full complex source/sink matrix',
         'Fixed diagonal scaling from the same configuration at t0=90, tau=0; all 16 ordered entries retained',
         'Identical vertical ranges. A ratio near one is a matrix change comparable to the baseline matrix norm.\n'
         'Weak late-time baselines can amplify relative changes. This is not a boundary-effect or potential significance test.',
         [Line2D([], [], color=COLORS[t], marker=MARKERS[t], label='%d versus 90' % t) for t in (58, 122)])


def source_sensitivity(rows, output):
    fig, axes = plt.subplots(2, 4, figsize=(15.6, 8.5))
    candidates = [(b, a) for b in ('pair', 'triple', 'full') for a in (4, 6, 8)]
    cfgstyles = {170: '-', 640: '--', 650: ':'}
    for i, tau in enumerate((16, 24)):
        for j, radius in enumerate(RADII):
            ax = axes[i, j]
            for cfg in CONFIGS:
                for t0 in (58, 122):
                    y = [number(matching(rows, configuration=cfg, R=radius, basis=b, anchor=a, t0=t0)[tau],
                                'common_t90_normalized_delta_re') for b, a in candidates]
                    ax.plot(np.arange(9) + (-.06 if t0 == 58 else .06), y,
                            color=COLORS[t0], ls=cfgstyles[cfg], marker=MARKERS[t0], ms=3, lw=.9, alpha=.8)
            ax.axhline(0, color='.5', lw=.65)
            ax.grid(alpha=.18)
            ax.set_title(r'$R/a={},\ \tau/a={}$'.format(radius, tau))
            ax.set_xticks(range(9))
            ax.set_xticklabels(['%s %d' % (b, a) for b, a in candidates], rotation=65, fontsize=8)
            if j == 0:
                ax.set_ylabel('Shift minus baseline / baseline Re Q(0)')
    handles = [Line2D([], [], color=COLORS[t], marker=MARKERS[t], label='%d minus 90' % t) for t in (58, 122)]
    handles += [Line2D([], [], color='.3', ls=cfgstyles[c], label='n%d' % c) for c in CONFIGS]
    save(fig, output, 'frozen_basis_training_controls', 'Position changes across all frozen source candidates',
         'Candidate labels give basis and training separation; no source is retuned or selected by shifted data',
         'Two display checkpoints; the complete 36-candidate, 33-time inventory is in position_projections.csv.\n'
         'Each line is one configuration and shift. Visual source stability is separate from ensemble significance.', handles)


def controls(rows, output):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    palette = ['#0072B2', '#D55E00', '#7B3294', '#333333']
    handles = []
    for radius, color, marker in zip(RADII, palette, ('o', 's', '^', 'D')):
        a = matching(rows, R=radius)
        x = [int(r['tau']) for r in a]
        axes[0].plot(x, [number(r, 'maximum_absolute_component_difference') / 1e-24 for r in a], color=color, marker=marker, ms=3, lw=.9)
        axes[1].plot(x, [number(r, 'relative_frobenius_difference') / 1e-15 for r in a], color=color, marker=marker, ms=3, lw=.9)
        handles.append(Line2D([], [], color=color, marker=marker, label=r'$R/a={}$'.format(radius)))
    axes[0].set_ylabel(r'Maximum component difference [$10^{-24}$]')
    axes[1].set_ylabel(r'Relative matrix-norm difference [$10^{-15}$]')
    for ax in axes:
        ax.set_xlabel(r'Correlator separation $\tau/a$')
        ax.grid(alpha=.18)
    save(fig, output, 'shifted_projection_control', 'Shifted piter=10 versus piter=40 control',
         'n650 at t0/a=122; axis-averaged full complex matrices at identical physical endpoints',
         'All 8,448 native ordered records pass the original component check; these panels show the axis-averaged subset.\n'
         'This checks projection iteration accuracy, not source-position independence.', handles)


def timings(rows, output):
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    kinds = [('audit', '#444444'), ('reference90_p10', '#999999'),
             ('shift58_p10', COLORS[58]), ('shift122_p10', COLORS[122]), ('shift122_p40', '#D55E00')]
    bottom = np.zeros(3)
    handles = []
    for stem, color in kinds:
        values = [sum(number(r, 'seconds') / 60 for r in rows if r['configuration'] == str(cfg)
                      and (r['kind'] == 'audit' if stem == 'audit' else r['stem'] == stem)) for cfg in CONFIGS]
        ax.bar(range(3), values, bottom=bottom, color=color, label=stem)
        bottom += values
        handles.append(plt.Rectangle((0, 0), 1, 1, color=color, label=stem.replace('_', ' ')))
    for j, value in enumerate(bottom):
        ax.text(j, value + .3, '%.2f min' % value, ha='center')
    ax.set_xticks(range(3))
    ax.set_xticklabels(['n%d' % c for c in CONFIGS])
    ax.set_ylabel('Sum of recorded MPI launch times [minutes]')
    ax.set_ylim(0, 27)
    ax.grid(axis='y', alpha=.15)
    save(fig, output, 'native_launch_timings', 'Measured MPI execution times for the completed pilot',
         'Three audit launches and eight scan launches; every recorded exit code is zero',
         'Stacked MPI durations exclude queue waiting, shell/setup overhead and gaps between launches.\n'
         'Task start times need not coincide; these bars are not the total elapsed campaign duration.', handles)


def render(output):
    output = Path(output)
    rows = csv_read(output / 'position_projections.csv')
    projected_curves(rows, output)
    projected_curves(rows, output, own=True)
    phases(rows, output)
    differences(rows, output)
    source_sensitivity(rows, output)
    matrices(csv_read(output / 'matrix_position_diagnostics.csv'), output)
    controls(csv_read(output / 'projection_control.csv'), output)
    timings(csv_read(output / 'launch_timings.csv'), output)
    print('Rendered eight PNG/SVG figure pairs; all pilot configurations retained.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(HERE / 'results'))
    render(parser.parse_args().output)
