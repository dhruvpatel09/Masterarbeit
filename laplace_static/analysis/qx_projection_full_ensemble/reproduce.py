#!/usr/bin/env python3
"""Reproduce the complete qX projection study from preserved measurements."""
import argparse
import csv
import hashlib
import importlib.util
import json
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUTS = HERE/'inputs'
RESULTS = HERE/'results'
JOB = INPUTS/'job5941300'
EXTRA_PLOTS = ['profiles_intermediate_precision_70configs',
               'n100_inclusion_70configs']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


analysis = module('qx_complete70_analysis', HERE/'make_plots.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_inputs():
    sources = json.loads((INPUTS/'SOURCES.json').read_text())
    for name, expected in sources['files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Preserved source changed: '+name)
    validation = json.loads((INPUTS/'job5941300_validation.json').read_text())
    if (validation['success'], validation['cases'], validation['mpi_launches'],
            validation['raw_records'], validation['window_records']) != (True, 123, 492, 163344, 492):
        raise ValueError('Incomplete full70 measurement validation')
    if validation['configurations'] != [100]+list(range(240, 631, 10)):
        raise ValueError('Wrong newly measured configuration cohort')
    if validation['piters'] != analysis.PITERS or validation['taus'] != [2, 4, 6]:
        raise ValueError('Wrong projection prescriptions or temporal windows')
    return sources


def joined_raw(sources):
    rows = [row for name in sources['raw_profile_tables']
            for row in analysis.core.read_csv(ROOT/name)]
    rows.sort(key=lambda r: tuple(int(r[k]) for k in
                                 ('configuration', 'tau', 'piter', 'axis', 'R', 'long', 'tr1', 'tr2')))
    if len(rows) != 278880:
        raise ValueError('Incomplete assembled raw profile table')
    return rows


def diagnostics(records, windows, without100):
    D = Decimal
    keys = analysis.KEYS
    combined = [r for r in records if int(r['axis']) == -1]
    bykey = {tuple(r[k] for k in keys): r for r in without100}

    def position(row):
        return {k: int(row[k]) for k in keys}

    with localcontext() as context:
        context.prec = 50
        def p4_ratio(row, field):
            return abs(D(row['delta_C_p4_minus_p40']))/D(row[field])

        def inclusion_ratio(row):
            other = bykey[tuple(row[k] for k in keys)]
            return abs(D(row['C_p40'])-D(other['C_p40']))/D(row['C_p40_jackknife_se'])

        profile = max(combined, key=lambda r: p4_ratio(r, 'C_p40_jackknife_se'))
        paired = max(combined, key=lambda r: p4_ratio(r, 'delta_C_paired_jackknife_se'))
        inclusion = max(combined, key=inclusion_ratio)
        rms = max((r for r in windows if int(r['piter']) == 4), key=lambda r: D(r['unit_rms']))
        return dict(
            configurations=70, combined_profile_positions=len(combined),
            max_abs_p4_change_over_p40_preliminary_se=dict(
                value=str(p4_ratio(profile, 'C_p40_jackknife_se')), position=position(profile)),
            max_abs_p4_change_over_paired_preliminary_se=dict(
                value=str(p4_ratio(paired, 'delta_C_paired_jackknife_se')), position=position(paired)),
            max_abs_n100_inclusion_change_over_p40_preliminary_se=dict(
                value=str(inclusion_ratio(inclusion)), position=position(inclusion)),
            largest_p4_window_unit_rms=dict(
                value=rms['unit_rms'], configuration=int(rms['configuration']), tau=int(rms['tau'])),
            central_window_unit_above_1em2=sum(int(r['unit_above_1em2']) for r in windows),
            central_window_det_above_1em2=sum(int(r['det_above_1em2']) for r in windows),
            interpretation='The 70- and 69-configuration estimates are correlated. The n100 '
                           'inclusion change divided by the full70 standard error is a scale '
                           'comparison, not an independent-sample significance test. All '
                           'standard errors are preliminary and lack autocorrelation correction. '
                           'Maxima across probe positions are descriptive, without a multiple-test correction.')


def plot_extras(records, without100, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.hashsalt': 'qx-full70-additional'})
    D = Decimal
    keys = analysis.KEYS
    bykey = {tuple(r[k] for k in keys): r for r in without100}

    def cut(tau):
        return sorted([r for r in records if tuple(int(r[k]) for k in
                      ('tau', 'axis', 'R', 'long', 'tr2')) == (tau, -1, 4, 2, 0)],
                      key=lambda r: int(r['tr1']))

    def save(fig, name):
        fig.savefig(str(output/(name+'.png')), dpi=180)
        fig.savefig(str(output/(name+'.svg')), metadata={'Date': None})
        path = output/(name+'.svg')
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, layout='constrained')
    for column, tau in enumerate((2, 4, 6)):
        rows = cut(tau)
        x = [int(r['tr1']) for r in rows]
        for line, p in enumerate((10, 20)):
            ax = axes[line, column]
            ax.plot(x, [float(r['delta_C_p%d_minus_p40' % p]) for r in rows],
                    'o-', color=analysis.COLORS[p], ms=5)
            ax.axhline(0, color='#87919b', lw=.8)
            ax.grid(alpha=.18)
            ax.set_xticks([-2, -1, 0, 1, 2])
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
            if column == 0:
                ax.set_ylabel(r'$C_{%d}-C_{40}$' % p)
            if line == 0:
                ax.set_title(r'$\tau/a=%d$' % tau)
            else:
                ax.set_xlabel(r'$r_\perp/a$')
    fig.suptitle(analysis.COHORT+'\nIntermediate prescriptions: differences at numerical precision', fontsize=14)
    fig.supxlabel(r'$R/a=4$, longitudinal midpoint; raw clover. '
                  'The displayed differences are at the binary64 rounding scale.', fontsize=10)
    save(fig, EXTRA_PLOTS[0])

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.1), sharex=True, layout='constrained')
    for column, tau in enumerate((2, 4, 6)):
        rows = cut(tau)
        old = [bykey[tuple(r[k] for k in keys)] for r in rows]
        x = [int(r['tr1']) for r in rows]
        for selected, label, color, marker, style in [
                (old, '69: excluding n100', '#87919b', 's', '--'),
                (rows, '70: including n100', analysis.COLORS[40], 'o', '-')]:
            axes[0, column].errorbar(x, [float(r['C_p40']) for r in selected],
                                    yerr=[float(r['C_p40_jackknife_se']) for r in selected],
                                    color=color, marker=marker, ls=style, capsize=3, label=label)
        with localcontext() as context:
            context.prec = 50
            change = [float((D(r['C_p40'])-D(o['C_p40']))/D(r['C_p40_jackknife_se']))
                      for r, o in zip(rows, old)]
        axes[1, column].plot(x, change, 'o-', color=analysis.COLORS[40])
        axes[1, column].axhline(0, color='#87919b', lw=.8)
        axes[1, column].set_ylim(-.35, .35)
        axes[1, column].set_xlabel(r'$r_\perp/a$')
        axes[0, column].set_title(r'$\tau/a=%d$' % tau)
        axes[0, column].legend(fontsize=9)
        for ax in axes[:, column]:
            ax.grid(alpha=.18)
            ax.set_xticks([-2, -1, 0, 1, 2])
    axes[0, 0].set_ylabel(r'$C_{40}$ with preliminary SE')
    axes[1, 0].set_ylabel(r'$(C_{70}-C_{69})/\mathrm{SE}_{70}$')
    fig.suptitle('Effect of including the validated n100 configuration\n'+
                 r'piter=40; $R/a=4$, longitudinal midpoint', fontsize=14)
    fig.supxlabel('Correlated cohorts; normalized change is descriptive, not a significance test. '
                  'Jackknife errors lack autocorrelation correction.', fontsize=9)
    save(fig, EXTRA_PLOTS[1])


def reproduce(output=RESULTS, make_plots=True):
    sources = verify_inputs()
    output = Path(output)
    summary = analysis.export_results([JOB/'profiles_raw.csv'], [JOB/'su3_windows.csv'],
                                      output, make_plots=False)
    raw = joined_raw(sources)
    analysis.core.write_csv(output/'profiles_raw_70configs.csv', raw)
    records = analysis.core.read_csv(output/'profiles_ensemble_70configs.csv')
    windows = analysis.core.read_csv(output/'su3_central_windows_70configs.csv')
    without100 = analysis.core.read_csv(output/'profiles_ensemble_69configs_without_n100.csv')
    report = diagnostics(records, windows, without100)
    (output/'cohort_diagnostics_70configs.json').write_text(json.dumps(report, indent=2)+'\n')
    if make_plots:
        renderer = module('qx_complete70_renderer', HERE/'render_plots.py')
        typed = [{k: (int(v) if k in analysis.KEYS+('configurations',) else Decimal(v))
                  for k, v in row.items()} for row in records]
        renderer.plot_results(typed, windows, output, analysis)
        plot_extras(records, without100, output)
    print('PASS assembled raw table and n100 inclusion diagnostics; '+
          ('seven figure pairs' if make_plots else 'numerical exports only'))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=RESULTS)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    reproduce(args.output, not args.no_plots)


if __name__ == '__main__':
    main()
