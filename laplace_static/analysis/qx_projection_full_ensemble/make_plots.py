#!/usr/bin/env python3
"""Combine all matched profiles into one 70-configuration ensemble."""
import argparse
import importlib.util
import json
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = [n for n in range(100, 801, 10) if n != 640]
CONFIGS = TARGET
PITERS = [4, 10, 20, 40]
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
D = Decimal
COLORS = {4: '#cd5519', 10: '#9c5f9b', 20: '#278260', 40: '#176b8e'}
COHORT = '70 matched configurations: n100-800, excluding n640'
PLOTS = [
    'profiles_all_piters_70configs',
    'su3_central_windows_all_piters_70configs',
    'profiles_p4_vs_p40_70configs',
    'profiles_preliminary_errors_70configs',
    'su3_central_windows_70configs',
]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


core = module('qx_full_core', ROOT/'laplace_static/analysis/qx_hyp_projection/make_plots.py')


def paired_estimator():
    return core.paired_estimator()


def aggregate_records(raw_rows):
    """Return the ratio of equal raw configuration means, never mean ratios."""
    groups = defaultdict(dict)
    configs, piters = set(), set()
    for row in raw_rows:
        cfg, piter = int(row['configuration']), int(row['piter'])
        if cfg not in TARGET:
            raise ValueError('Configuration outside the full ensemble')
        if piter < 1:
            raise ValueError('Projection iterations must be positive')
        key = tuple(int(row[k]) for k in KEYS)
        ident = cfg, piter
        if ident in groups[key]:
            raise ValueError('Duplicate configuration/piter/profile record')
        values = tuple(D(row[f]) for f in ('L_RE', 'LS_RE', 'S_RAW'))
        if any(not v.is_finite() for v in values) or values[0] == 0:
            raise ValueError('Non-finite measurement or zero loop denominator')
        groups[key][ident] = values
        configs.add(cfg)
        piters.add(piter)
    if not groups or not {4, 40}.issubset(piters):
        raise ValueError('Need measured piter4 and piter40 records')
    expected = {(c, p) for c in configs for p in piters}
    records = []
    with localcontext() as context:
        context.prec = 50
        for key, values in sorted(groups.items()):
            if set(values) != expected:
                raise ValueError('Incomplete matched cohort at '+str(key))
            for cfg in configs:
                if len({values[cfg, p][2] for p in piters}) != 1:
                    raise ValueError('Original-field clover changed between piters')
            result = dict(zip(KEYS, key))
            result['configurations'] = len(configs)
            for p in sorted(piters):
                totals = [sum((values[c, p][j] for c in configs), D(0)) for j in range(3)]
                means = [value/D(len(configs)) for value in totals]
                if means[0] == 0:
                    raise ValueError('Undefined ensemble ratio')
                result.update({f'mean_L_RE_p{p}': means[0], f'mean_LS_RE_p{p}': means[1],
                               f'mean_S_RAW_p{p}': means[2], f'C_p{p}': totals[1]/totals[0]-means[2]})
            for p in sorted(piters-{40}):
                result[f'delta_C_p{p}_minus_p40'] = result[f'C_p{p}']-result['C_p40']
            if len(configs) > 1:
                statistics = paired_estimator()(values, sorted(configs))
                for field in ('C_p4_jackknife_se', 'C_p40_jackknife_se', 'delta_C_paired_jackknife_se'):
                    result[field] = statistics[field]
            records.append(result)
    return records, sorted(configs), sorted(piters)

def load_data(extra_raw, extra_windows):
    workflow = module('qx_full_workflow_analysis', ROOT/'laplace_static/validation/qx_projection_full_ensemble/campaign.py')
    workflow.verify_sources()
    previous = module('qx_full_base29', ROOT/'laplace_static/analysis/qx_projection_spectrum_29configs/make_plots.py')
    _, raw, windows, _, _ = previous.load_data()
    for path in extra_raw:
        raw += core.read_csv(Path(path))
    for path in extra_windows:
        windows += core.read_csv(Path(path))
    original = core.check_sources()
    expected_keys = {tuple(int(r[k]) for k in KEYS)
                     for r in core.read_csv(ROOT/original['default_profile_tables'][0]['path'])}
    core.validate_geometry(raw, expected_keys)
    records, configs, piters = aggregate_records(raw)
    core.validate_windows(windows, configs, piters)
    if configs != TARGET or piters != PITERS or len(raw) != 278880 or len(windows) != 840 or len(records) != 996:
        raise ValueError('Incomplete full ensemble: require all 70 configurations, three windows and four prescriptions')
    return raw, windows, records


def export_results(extra_raw, extra_windows, output, make_plots=True):
    raw, windows, records = load_data(extra_raw, extra_windows)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    core.write_csv(output/'profiles_ensemble_70configs.csv', records)
    core.write_csv(output/'su3_central_windows_70configs.csv', windows)
    comparison = module('qx_full_comparison13', ROOT/'laplace_static/analysis/qx_projection_spectrum_13configs/make_plots.py')
    comparisons = comparison.comparison_summary(records, windows)
    for row in comparisons:
        row['configurations'] = 70
    core.write_csv(output/'piter_comparison_70configs.csv', comparisons)
    without100, configs69, _ = aggregate_records([r for r in raw if int(r['configuration']) != 100])
    core.write_csv(output/'profiles_ensemble_69configs_without_n100.csv', without100)
    bykey = {tuple(r[k] for k in KEYS): r for r in without100}
    influence = []
    with localcontext() as context:
        context.prec = 50
        for row in records:
            if row['axis'] != -1:
                continue
            other = bykey[tuple(row[k] for k in KEYS)]
            result = {k: row[k] for k in KEYS}
            for p in PITERS:
                field = 'C_p'+str(p)
                result[field+'_70configs'] = row[field]
                result[field+'_69configs'] = other[field]
                result[field+'_including_minus_excluding_n100'] = row[field]-other[field]
            influence.append(result)
    core.write_csv(output/'n100_inclusion_comparison.csv', influence)
    summary = dict(success=True, configurations=TARGET, configuration_count=70,
                   configuration_window_cases=210, profile_piters=PITERS, taus=[2, 4, 6],
                   raw_records=len(raw), window_records=len(windows), profile_positions=len(records),
                   combined_profile_positions=len(influence), piter_comparison=comparisons,
                   reused_configurations=29, newly_measured_configurations=41,
                   historical_without_n100_configurations=configs69,
                   estimator='mean(LS_RE)/mean(L_RE)-mean(S_RAW); equal configuration weights; precision50',
                   uncertainties='Paired delete-one-configuration jackknife for piter4/40; preliminary; no autocorrelation correction',
                   clover_convention='Raw clover; no colour-trace subtraction or additional factor 1/2',
                   scope='All 70 generated configurations form one matched ensemble. n100 was previously '
                         'used only as a regression control; it now has all four prescriptions at all '
                         'three windows. The 69-configuration export records the effect of this explicit '
                         'cohort change. n640 remains unavailable.',
                   imported_raw=[dict(name=Path(p).name, sha256=core.sha256(Path(p))) for p in extra_raw],
                   imported_windows=[dict(name=Path(p).name, sha256=core.sha256(Path(p))) for p in extra_windows])
    (output/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    if make_plots:
        plot_results(records, windows, output)
    print('PASS full ensemble: 70 configurations, 210 cases, 278880 raw records, 840 window rows')
    return summary


def plot_results(records, windows, output):
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
            ax.set_xticks(CONFIGS[::4]+[CONFIGS[-1]], rotation=45, fontsize=8)
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

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extra-raw', type=Path, action='append', required=True)
    parser.add_argument('--extra-windows', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, default=HERE/'results')
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    export_results(args.extra_raw, args.extra_windows, args.output, not args.no_plots)


if __name__ == '__main__':
    main()
