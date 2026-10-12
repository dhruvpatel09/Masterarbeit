#!/usr/bin/env python3
"""Plot qX temporal-HYP projection diagnostics.

Configuration weights are equal. Ratios are formed after averaging raw L, LS,
and S. All plotted positions and prescriptions must share a complete cohort.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import re
from collections import defaultdict
from decimal import Decimal, localcontext
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
D = Decimal


@lru_cache(maxsize=1)
def paired_estimator():
    source = json.loads((HERE/'inputs/SOURCES.json').read_text())['estimator_module']
    path = REPO/source['path']
    if sha256(path) != source['sha256']:
        raise ValueError('Estimator module digest changed')
    spec = importlib.util.spec_from_file_location('hyp_projection_estimators', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.paired_estimate


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check_sources():
    sources = json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name, expected in sources['preserved_inputs'].items():
        if sha256(HERE/'inputs'/name) != expected:
            raise ValueError('Preserved input digest changed: '+name)
    historical = sources['historical_80']
    if sha256(REPO/historical['source']) != historical['source_sha256']:
        raise ValueError('Historical piter80 source changed; review its recorded bound')
    for section in ('default_profile_tables', 'default_window_tables'):
        for source in sources[section]:
            if sha256(REPO/source['path']) != source['sha256']:
                raise ValueError('Source digest changed: '+source['path'])
    return sources


def aggregate(raw_rows):
    """Return the ratio of equal raw configuration means, never mean ratios."""
    groups = defaultdict(dict)
    configs, piters = set(), set()
    for row in raw_rows:
        cfg, piter = int(row['configuration']), int(row['piter'])
        if cfg == 100:
            raise ValueError('n100 is a regression control; exclude it from this ensemble')
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


def validate_geometry(raw_rows, expected_keys):
    actual_keys = set()
    loops = {}
    for row in raw_rows:
        key = tuple(int(row[k]) for k in KEYS)
        tau, axis, separation, along, _, _ = key
        if (int(row['Nv']), int(row['tmid']), int(row['t0']), int(row['t1'])) != (100, 98, 98-tau//2, 98+tau//2):
            raise ValueError('Changed mode count or time geometry')
        if int(row['Nsrc']) != (786432 if axis == -1 else 262144):
            raise ValueError('Changed spatial sample count')
        if D(row['long_from_midpoint']) != D(along)-D(separation)/2:
            raise ValueError('Changed longitudinal coordinate')
        actual_keys.add(key)
        loop_key = (int(row['configuration']), int(row['piter']), tau, axis, separation)
        value = D(row['L_RE'])
        if loops.setdefault(loop_key, value) != value:
            raise ValueError('Loop changed across probe positions')
    if actual_keys != expected_keys:
        raise ValueError('Incomplete or changed probe geometry; expected the preserved 996 keys')


def validate_windows(rows, configs, piters):
    expected = {(c, t, p) for c in configs for t in (2, 4, 6) for p in piters}
    seen = set()
    for row in rows:
        key = int(row['configuration']), int(row['tau']), int(row['piter'])
        if key in seen:
            raise ValueError('Duplicate central-window record')
        seen.add(key)
        _, tau, _ = key
        if (int(row['t0']), int(row['t1']), int(row['count'])) != (98-tau//2, 98+tau//2, tau*262144):
            raise ValueError('Changed SU(3) window geometry/count')
        for name in ('unit_max', 'unit_rms', 'det_max', 'det_rms'):
            value = D(row[name])
            if not value.is_finite() or value < 0:
                raise ValueError('Invalid SU(3) defect')
    if seen != expected:
        raise ValueError('Window and profile cohorts differ; supply all matched windows')


def historical_measurements(sources, times):
    # Use the SAME maximum entry metric as the historical p10/p20/p40 C tests.
    exact = []
    for p, stage in ((4, 'PITER4'), (40, 'PITER40')):
        selected = [r for r in times if r['stage'] == stage]
        exact.append({'piter': p, 'unit_max_entry': max(D(r['UNIT_MAX_ENTRY']) for r in selected),
                      'det_max': max(D(r['DET_MAX']) for r in selected), 'kind': 'measured', 'job': 5916659})
    for name, iterations, job in (('job5902645_diagnostics.txt', (10, 20), 5902645),
                                  ('job5903127_diagnostics.txt', (20, 40), 5903127)):
        body = (HERE/'inputs'/name).read_text()
        unit = next(line for line in body.splitlines() if line.startswith('UNITARITY '))
        det = next(line for line in body.splitlines() if line.startswith('DETERMINANT '))
        for p in iterations:
            u = D(re.search(r'piter'+str(p)+r'=([^\s]+)', unit).group(1))
            d = D(re.search(r'piter'+str(p)+r'=([^\s]+)', det).group(1))
            existing = next((r for r in exact if r['piter'] == p), None)
            if existing:
                # The time CSV serializes binary64 with fewer trailing digits
                # than the original %.16e log. Allow only that rounding loss.
                if any(abs(a-b) > max(abs(a), abs(b))*D('5e-16') for a, b in
                       ((existing['unit_max_entry'], u), (existing['det_max'], d))):
                    raise ValueError('Historical repeats disagree on final-link maxima')
            else:
                exact.append({'piter': p, 'unit_max_entry': u, 'det_max': d, 'kind': 'measured', 'job': job})
    bound = sources['historical_80']
    exact.append({'piter': 80, 'unit_max_entry': D(str(bound['unit_max_entry_upper_bound'])),
                  'det_max': D(str(bound['det_max_upper_bound'])), 'kind': 'strict_upper_bound', 'job': bound['job']})
    # Pair distances are component differences, NOT SU(3) defects.
    pairs = []
    for name, a, b, job in (('job5902645_diagnostics.txt', 10, 20, 5902645),
                            ('job5903127_diagnostics.txt', 20, 40, 5903127)):
        body = (HERE/'inputs'/name).read_text()
        match = re.search(r'PITER_\d+_VS_\d+ temporal=([^\s]+)', body)
        pairs.append({'piter_a': a, 'piter_b': b, 'component_max': D(match.group(1)), 'job': job})
    pairs.append({'piter_a': 40, 'piter_b': 80, 'component_max': D(bound['component_difference_40_vs_80']), 'job': bound['job']})
    return sorted(exact, key=lambda r: r['piter']), pairs


def plot_all(records, windows, times, history, pairs, configs, piters, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'qx-hyp-projection-summary'})
    colors = {2: '#176b8e', 4: '#cd5519', 6: '#437c45'}
    n = len(configs)
    cohort = f'n{configs[0]}–n{configs[-1]}' if configs == list(range(configs[0], configs[-1]+1, 10)) else ', '.join('n'+str(c) for c in configs)
    def save(fig, name):
        fig.savefig(output/(name+'.png'), dpi=180)
        fig.savefig(output/(name+'.svg'), metadata={'Date': None})
        svg_path = output/(name+'.svg')
        svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
        plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.3), sharex=True, gridspec_kw={'height_ratios': [1.5, 1]}, layout='constrained')
    for tau in (2, 4, 6):
        rows = sorted((r for r in records if (r['tau'], r['axis'], r['R'], r['long'], r['tr2']) == (tau, -1, 4, 2, 0)), key=lambda r: r['tr1'])
        x = [r['tr1'] for r in rows]
        axes[0].plot(x, [float(r['C_p40']) for r in rows], 'o-', color=colors[tau], label=rf'$\tau/a={tau}$, piter=40')
        axes[0].plot(x, [float(r['C_p4']) for r in rows], 'x--', color=colors[tau], label=rf'$\tau/a={tau}$, piter=4')
        axes[1].plot(x, [float(r['delta_C_p4_minus_p40']) for r in rows], 'o-', color=colors[tau], label=rf'$\tau/a={tau}$')
    for ax in axes:
        ax.grid(alpha=.18)
        ax.axhline(0, color='#87919b', lw=.8)
    axes[0].set_ylabel(r'$C=\overline{LS}_{\mathrm{Re}}/\overline{L}_{\mathrm{Re}}-\overline{S}_{\mathrm{RAW}}$')
    axes[0].legend(ncol=3, fontsize=10, loc='center')
    axes[1].set_ylabel(r'$C_4-C_{40}$')
    axes[1].legend(ncol=3, fontsize=10, loc='upper right')
    axes[1].set_xlabel(r'Transverse displacement $r_\perp/a$')
    axes[1].set_xticks([-2, -1, 0, 1, 2])
    fig.suptitle(f'Raw paired profiles: {n} configurations ({cohort})\n'+r'$R/a=4$, longitudinal midpoint $r_\parallel/a=2$; combined axes; no final statistical error bars', fontsize=14)
    save(fig, f'profiles_p4_vs_p40_{n}configs')

    if n > 1:
        fig, axes = plt.subplots(2, 1, figsize=(12.5, 7.3), sharex=True, layout='constrained')
        for tau in (2, 4, 6):
            rows = sorted((r for r in records if (r['tau'], r['axis'], r['R'], r['long'], r['tr2']) == (tau, -1, 4, 2, 0)), key=lambda r: r['tr1'])
            x = [r['tr1'] for r in rows]
            axes[0].errorbar(x, [float(r['C_p40']) for r in rows],
                             yerr=[float(r['C_p40_jackknife_se']) for r in rows], fmt='o-',
                             color=colors[tau], capsize=4, label=rf'$\tau/a={tau}$')
            axes[1].errorbar(x, [float(r['delta_C_p4_minus_p40']) for r in rows],
                             yerr=[float(r['delta_C_paired_jackknife_se']) for r in rows], fmt='o-',
                             color=colors[tau], capsize=4, label=rf'$\tau/a={tau}$')
        for ax in axes:
            ax.grid(alpha=.18)
            ax.axhline(0, color='#87919b', lw=.8)
            ax.legend(ncol=3, fontsize=10)
        axes[0].set_ylabel(r'$C_{40}$ ± preliminary SE')
        axes[1].set_ylabel(r'$C_4-C_{40}$ ± paired SE')
        axes[1].set_xlabel(r'Transverse displacement $r_\perp/a$')
        axes[1].set_xticks([-2, -1, 0, 1, 2])
        fig.suptitle(f'Raw profiles and projection differences: {n} configurations ({cohort})\n'+
                     r'$R/a=4$, longitudinal midpoint; matched delete-one-configuration jackknife', fontsize=14)
        fig.supxlabel('Preliminary uncertainties; autocorrelation not assessed. Raw clover convention.', fontsize=11)
        save(fig, f'profiles_preliminary_errors_{n}configs')

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.6), layout='constrained')
    for col, field, label in ((0, 'UNIT_MAX_FROB', r'$\max_{\mathbf{x}}\|U^\dagger U-I\|_F$'),
                              (1, 'DET_MAX', r'$\max_{\mathbf{x}}|\det U-1|$')):
        for p, style, color in ((4, '-', colors[4]), (40, '--', colors[2])):
            rows = [r for r in times if r['stage'] == f'PITER{p}']
            for index in (0, 1):
                chosen = rows if index == 0 else [r for r in rows if 90 <= int(r['t']) <= 105]
                axes[index, col].plot([int(r['t']) for r in chosen], [float(r[field]) for r in chosen],
                                      style, marker='o' if index else None, ms=4, color=color, label=f'piter={p}')
        for index in (0, 1):
            ax = axes[index, col]
            ax.set_yscale('log')
            ax.set_ylim(3e-16, 4 if index == 0 else 3e-4)
            ax.set_yticks([10.**e for e in range(-14, 1 if index == 0 else -3, 2)])
            ax.grid(axis='y', alpha=.2)
            ax.axvspan(94.5, 100.5, color=colors[6], alpha=.12)
            ax.set_ylabel(label)
        axes[0, col].set_title('Full physical temporal-link range')
        axes[0, col].legend(fontsize=10, loc='lower right')
        axes[1, col].set_title(r'Central region; shaded: $\tau/a=6$ links at $t/a=95\ldots100$', fontsize=11)
        axes[1, col].set_xlabel(r'Temporal-link origin $t/a$')
        p4 = [r for r in times if r['stage'] == 'PITER4']
        for t, align, xt in ((0, 'left', .09), (190, 'right', .86)):
            value = float(next(r for r in p4 if int(r['t']) == t)[field])
            axes[0, col].annotate(f't={t}: {value:.2f}', xy=(t, value), xytext=(xt, .95),
                                   textcoords='axes fraction', ha=align, va='top', fontsize=11,
                                   arrowprops={'arrowstyle': '-', 'color': '#78828b'})
    fig.suptitle('n100: piter=4 defects are largest next to the open temporal boundaries\nAbsent outgoing link at t=191 excluded', fontsize=14)
    save(fig, 'su3_boundary_vs_centre_n100')

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), layout='constrained')
    for ax, maximum, rms, label in zip(axes, ('unit_max', 'det_max'), ('unit_rms', 'det_rms'),
                                       (r'Unitarity defect $\|U^\dagger U-I\|_F$', r'Determinant defect $|\det U-1|$')):
        for p, field, style, marker, color in ((4, maximum, '-', 'o', colors[4]),
                                               (4, rms, '--', 's', colors[6]),
                                               (40, maximum, ':', '^', colors[2])):
            y = [max(float(r[field]) for r in windows if int(r['configuration']) == c and int(r['piter']) == p) for c in configs]
            ax.plot(configs, y, style, marker=marker, ms=4, color=color,
                    label=f'piter={p}: largest window '+('RMS' if field == rms else 'maximum'))
        ax.axhline(.01, color='#78828b', ls='--', lw=1)
        ax.annotate('0.01 counter threshold (not an acceptance limit)', xy=(.08, .01),
                    xycoords=('axes fraction', 'data'), xytext=(0, 5), textcoords='offset points',
                    color='#78828b', va='bottom', fontsize=10,
                    bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .5})
        ax.set_yscale('log')
        ax.set_ylim(1e-16, .05)
        ax.set_xticks(configs)
        ax.tick_params(axis='x', rotation=45)
        ax.set_xlabel('Configuration n')
        ax.grid(axis='y', alpha=.2)
        ax.set_title(label, fontsize=13)
        ax.legend(loc='lower right', fontsize=9)
    fig.suptitle(f'Central transporter windows: {n} configurations; '+r'$\tau/a=2,4,6$; midpoint $t/a=98$', fontsize=14)
    save(fig, f'su3_central_windows_{n}configs')

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), layout='constrained')
    ticks = list(range(len(history)))
    for ax, field, title in ((axes[0], 'unit_max_entry', r'Maximum entry defect $\max_{x,i,j}|(U^\dagger U-I)_{ij}|$'),
                              (axes[1], 'det_max', r'Maximum determinant defect $\max_x|\det U-1|$')):
        for i, row in enumerate(history):
            bound = row['kind'] == 'strict_upper_bound'
            ax.scatter(i, float(row[field]), marker='v' if bound else 'o', s=65,
                       facecolors='none' if bound else colors[2], edgecolors=colors[2])
            ax.annotate(('<' if bound else '')+f'{float(row[field]):.2e}', (i, float(row[field])),
                        xytext=(0, 12 if i%2 == 0 else -20), textcoords='offset points', ha='center', fontsize=9)
        ax.set_xticks(ticks, [str(r['piter']) for r in history])
        ax.set_xlim(-.5, len(history)-.5)
        ax.set_xlabel('Projection iterations (piter)')
        ax.set_title(title, fontsize=11)
        ax.set_yscale('log')
        ax.set_ylim(8e-18, 100)
        ax.grid(axis='y', alpha=.2)
    labels = [f"{r['piter_a']} vs {r['piter_b']}" for r in pairs]
    values = [float(r['component_max']) for r in pairs]
    axes[2].scatter(range(len(pairs)), values, color=colors[4], s=65)
    for i, value in enumerate(values):
        axes[2].annotate(f'{value:.2e}', (i, value), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=10)
    axes[2].set_xticks(range(len(pairs)), labels)
    axes[2].set_xlim(-.5, len(pairs)-.5)
    axes[2].set_xlabel('Measured prescription pair')
    axes[2].set_title('Field stability: maximum real/imaginary\ncomponent difference between final links', fontsize=11)
    axes[2].set_yscale('log')
    axes[2].set_ylim(8e-17, 100)
    axes[2].grid(axis='y', alpha=.2)
    fig.suptitle('n100: projection accuracy and iteration stability\nOne temporal HYP step, α=(1,1,0.5); full physical temporal-link range', fontsize=14)
    fig.supxlabel('piter=80: open triangles are documented strict upper bounds (<8×10⁻¹⁶), not exact values.\nEntry norm here differs from the Frobenius norm in the time/window plots; these are separate measurements.', fontsize=10)
    save(fig, 'piter_spectrum_n100')

    # Intermediate-piter profile curves appear only when actual complete cohorts exist.
    if len(piters) > 2:
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, layout='constrained')
        for col, tau in enumerate((2, 4, 6)):
            rows = sorted((r for r in records if (r['tau'], r['axis'], r['R'], r['long'], r['tr2']) == (tau, -1, 4, 2, 0)), key=lambda r: r['tr1'])
            for p in piters:
                axes[0, col].plot([r['tr1'] for r in rows], [float(r[f'C_p{p}']) for r in rows], 'o-', ms=4, label=f'piter={p}')
                if p != 40:
                    axes[1, col].plot([r['tr1'] for r in rows], [float(r[f'delta_C_p{p}_minus_p40']) for r in rows], 'o-', ms=4, label=f'piter={p}')
            axes[0, col].set_title(rf'$\tau/a={tau}$')
            axes[0, col].legend(fontsize=9)
            axes[1, col].set_xlabel(r'$r_\perp/a$')
            for ax in axes[:, col]:
                ax.grid(alpha=.2)
        axes[0, 0].set_ylabel(r'$C_p$')
        axes[1, 0].set_ylabel(r'$C_p-C_{40}$')
        fig.suptitle(f'Measured intermediate-piter profiles: {n} matched configurations; R/a=4\nRaw clover diagnostic; no final statistical error bars', fontsize=14)
        save(fig, f'profiles_all_piters_{n}configs')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, action='append', help='Replace default raw sources; repeat for each table')
    parser.add_argument('--windows', type=Path, action='append', help='Replace default window sources; repeat for each table')
    parser.add_argument('--extra-raw', type=Path, action='append', default=[], help='Append a completed campaign raw table')
    parser.add_argument('--extra-windows', type=Path, action='append', default=[], help='Append its matching central-window table')
    parser.add_argument('--output', type=Path, default=HERE/'results')
    args = parser.parse_args()
    sources = check_sources()
    raw_paths = (args.raw or [REPO/s['path'] for s in sources['default_profile_tables']])+args.extra_raw
    window_paths = (args.windows or [REPO/s['path'] for s in sources['default_window_tables']])+args.extra_windows
    raw = [r for path in raw_paths for r in read_csv(path)]
    template = read_csv(REPO/sources['default_profile_tables'][0]['path'])
    expected_keys = {tuple(int(r[k]) for k in KEYS) for r in template}
    validate_geometry(raw, expected_keys)
    records, configs, piters = aggregate(raw)
    windows = [r for path in window_paths for r in read_csv(path)]
    validate_windows(windows, configs, piters)
    times = read_csv(HERE/'inputs/n100_su3_by_time.csv')
    if len(times) != 573 or {(r['stage'], int(r['t'])) for r in times} != {(s, t) for s in ('ORIGINAL', 'PITER4', 'PITER40') for t in range(191)}:
        raise ValueError('Incomplete n100 time scan or absent upper-boundary link included')
    if any(int(r['N']) != 262144 for r in times):
        raise ValueError('Changed n100 per-time sample count')
    history, pairs = historical_measurements(sources, times)
    combined = [r for r in records if r['axis'] == -1]
    with localcontext() as context:
        context.prec = 50
        maximum_abs = max(abs(r['delta_C_p4_minus_p40']) for r in combined)
        relative = [abs(r['delta_C_p4_minus_p40']/r['C_p40'])*100 for r in combined if r['C_p40']]
    summary = {'configurations': configs, 'cases': len(configs)*3, 'raw_records': len(raw),
               'combined_profile_positions': len(combined), 'profile_piters': piters,
               'estimator': 'mean(LS_RE)/mean(L_RE)-mean(S_RAW); equal configuration weights; precision50',
               'uncertainties': 'Matched delete-one-configuration jackknife for piter4/40; no autocorrelation correction; preliminary',
               'max_combined_abs_delta_C_4_40': str(maximum_abs),
               'max_combined_relative_delta_percent_4_40': str(max(relative)),
               'relative_denominator_zero_points': sum(r['C_p40'] == 0 for r in combined),
               'inputs': [{'path': str(p.relative_to(REPO)) if p.is_relative_to(REPO) else str(p), 'sha256': sha256(p)} for p in raw_paths+window_paths]}
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output/f'profiles_ensemble_{len(configs)}configs.csv', records)
    write_csv(args.output/f'su3_central_windows_{len(configs)}configs.csv', windows)
    write_csv(args.output/'piter_spectrum_n100.csv', history)
    write_csv(args.output/'piter_pair_stability_n100.csv', pairs)
    (args.output/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    plot_all(records, windows, times, history, pairs, configs, piters, args.output)
    print(json.dumps({'configurations': len(configs), 'profile_piters': piters, 'max_abs_delta': str(maximum_abs), 'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()
