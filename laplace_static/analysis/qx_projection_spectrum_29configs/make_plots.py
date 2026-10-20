#!/usr/bin/env python3
"""Reproduce the qX projection comparison on 29 matched configurations."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONFIGS = list(range(110, 231, 10))+list(range(650, 801, 10))
PITERS = [4, 10, 20, 40]
COLORS = {4: '#cd5519', 10: '#9c5f9b', 20: '#278260', 40: '#176b8e'}
COHORT = '29 configurations: n110-230 and n650-800'
PLOTS = [
    'profiles_all_piters_29configs',
    'su3_central_windows_all_piters_29configs',
    'profiles_p4_vs_p40_29configs',
    'profiles_preliminary_errors_29configs',
    'su3_central_windows_29configs',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def verify_inputs():
    sources = json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name, expected in sources['preserved_inputs'].items():
        if digest(HERE/'inputs'/name) != expected:
            raise ValueError('Preserved input changed: '+name)
    for name, expected in sources['upstream_files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Upstream numerical source changed: '+name)
    validation = json.loads((HERE/'inputs/job5940445_validation.json').read_text())
    if (validation['success'], validation['cases'], validation['mpi_launches'],
            validation['raw_records'], validation['window_records']) != (True, 9, 36, 11952, 36):
        raise ValueError('Incomplete n780-800 validation')
    if validation['configurations'] != [780, 790, 800] or validation['piters'] != PITERS:
        raise ValueError('Wrong n780-800 validation cohort')
    return sources


def load_data():
    sources = verify_inputs()
    previous = module('qx_projection_base26', ROOT/'laplace_static/analysis/qx_projection_spectrum_26configs/make_plots.py')
    core, raw, windows, _, _ = previous.load_data()
    raw += core.read_csv(HERE/'inputs/job5940445_profiles_raw.csv')
    windows += core.read_csv(HERE/'inputs/job5940445_su3_windows.csv')
    original = core.check_sources()
    expected_keys = {tuple(int(r[k]) for k in core.KEYS)
                     for r in core.read_csv(ROOT/original['default_profile_tables'][0]['path'])}
    core.validate_geometry(raw, expected_keys)
    records, configs, piters = core.aggregate(raw)
    core.validate_windows(windows, configs, piters)
    if configs != CONFIGS or piters != PITERS or len(raw) != 115536 or len(windows) != 348:
        raise ValueError('Changed matched 29-configuration cohort')
    return core, raw, windows, records, sources


def comparisons(core, records, windows):
    previous = module('qx_comparison13_for29', ROOT/'laplace_static/analysis/qx_projection_spectrum_13configs/make_plots.py')
    result = previous.comparison_summary(records, windows)
    for row in result:
        row['configurations'] = len(CONFIGS)
    return result


def summary_data(raw, windows, records, comparisons, sources):
    return dict(
        success=True, configurations=CONFIGS, configuration_ranges=[[110, 230], [650, 800]],
        cases=87, raw_records=len(raw), total_window_records=len(windows),
        combined_profile_positions=sum(int(r['axis']) == -1 for r in records),
        profile_positions=len(records), profile_piters=PITERS,
        estimator='mean(LS_RE)/mean(L_RE)-mean(S_RAW); equal configuration weights; precision50',
        clover_convention='Raw clover; no colour-trace subtraction or additional factor 1/2',
        uncertainties='Matched delete-one-configuration jackknife for piter4/40; no autocorrelation correction; preliminary',
        piter_comparison=comparisons, completed_expansion_jobs=[5927227, 5940445],
        latest_expansion_cases=9, latest_expansion_launches=36,
        preserved_inputs=sources['preserved_inputs'],
        source_review_archive_sha256=sources['review_archive_sha256'],
        source_review_archive_incomplete=False,
        scope='The latest 36 logs and combined exports were revalidated. Earlier 26-configuration inputs '
              'are preserved unchanged. Only n110-230 and n650-800 enter these matched profiles; '
              'n100 is a regression control and n640 remains unresolved.')


def plot_results(records, windows, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
        'svg.hashsalt': 'qx-projection-spectrum26'})

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
            ax.set_xticks(list(range(len(CONFIGS))), [str(c) for c in CONFIGS], rotation=45, fontsize=8)
            ax.axvline(12.5, color='#87919b', ls=':', lw=1)
            ax.set_xlabel('Configuration n (two sampled ranges)')
        return fig, axes

    def window_curve(ax, values, label, color, style='-', marker='o'):
        for start, stop in ((0, 13), (13, len(CONFIGS))):
            ax.plot(list(range(start, stop)), values[start:stop], style, marker=marker,
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
    parser.add_argument('--output', type=Path, default=HERE/'results')
    args = parser.parse_args()
    core, raw, windows, records, sources = load_data()
    args.output.mkdir(parents=True, exist_ok=True)
    core.write_csv(args.output/'profiles_ensemble_29configs.csv', records)
    core.write_csv(args.output/'su3_central_windows_29configs.csv', windows)
    comparison = comparisons(core, records, windows)
    core.write_csv(args.output/'piter_comparison_29configs.csv', comparison)
    summary = summary_data(raw, windows, records, comparison, sources)
    (args.output/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    plot_results(records, windows, args.output)
    print('PASS: complete 29-configuration projection spectrum:', args.output)


if __name__ == '__main__':
    main()
