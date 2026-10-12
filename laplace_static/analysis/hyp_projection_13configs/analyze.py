#!/usr/bin/env python3
"""Paired projection comparison, ratio of ensemble means, decimal arithmetic."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

CONFIGS = tuple(range(110, 231, 10))
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
FIELDS = ('L_RE', 'LS_RE', 'S_RAW')
D = Decimal


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_inputs(base):
    sources = json.loads((base / 'SOURCES.json').read_text())
    groups = defaultdict(dict)
    windows = []
    for source in sources['profile_sources']:
        for name, expected in source['tables'].items():
            if digest(base / name) != expected:
                raise ValueError('Input digest changed: ' + name)
        with (base / source['raw_table']).open(newline='') as handle:
            for row in csv.DictReader(handle):
                cfg, piter = int(row['configuration']), int(row['piter'])
                key = tuple(int(row[k]) for k in KEYS)
                tau, axis, separation, along, tr1, tr2 = key
                if cfg not in CONFIGS or piter not in (4, 40):
                    raise ValueError('Unexpected configuration or projection count')
                if (int(row['Nv']), int(row['tmid']), int(row['t0']), int(row['t1'])) != (100, 98, 98-tau//2, 98+tau//2):
                    raise ValueError('Changed mode count or time window')
                if tau not in (2, 4, 6) or axis not in (-1, 0, 1, 2) or separation not in (1, 2, 4):
                    raise ValueError('Unexpected observable geometry')
                if int(row['Nsrc']) != (786432 if axis == -1 else 262144):
                    raise ValueError('Unexpected spatial sample count')
                if D(row['long_from_midpoint']) != D(along)-D(separation)/2:
                    raise ValueError('Unexpected longitudinal coordinate')
                if (cfg, piter) in groups[key]:
                    raise ValueError('Duplicate raw profile record')
                values = tuple(D(row[f]) for f in FIELDS)
                if any(not v.is_finite() for v in values) or values[0] == 0:
                    raise ValueError('Non-finite values or zero loop denominator')
                groups[key][cfg, piter] = values
        with (base / source['window_table']).open(newline='') as handle:
            windows.extend(csv.DictReader(handle))
    expected = {(cfg, p) for cfg in CONFIGS for p in (4, 40)}
    if len(groups) != 996:
        raise ValueError('Expected 996 profile keys')
    for key, rows in groups.items():
        if set(rows) != expected:
            raise ValueError('Missing paired configuration at ' + str(key))
        for cfg in CONFIGS:
            if rows[cfg, 4][2] != rows[cfg, 40][2]:
                raise ValueError('Original-field clover differs between prescriptions')
    # The loop denominator must be the same at every probe position.
    loops = {}
    for key, rows in groups.items():
        tau, axis, separation = key[:3]
        for (cfg, piter), values in rows.items():
            loop_key = (cfg, piter, tau, axis, separation)
            if loops.setdefault(loop_key, values[0]) != values[0]:
                raise ValueError('Loop denominator changed across probe positions')
    if len(windows) != 78:
        raise ValueError('Expected 78 SU(3) window records')
    return sources, groups, windows


def estimator(values):
    count = D(len(values))
    sums = tuple(sum((v[j] for v in values), D(0)) for j in range(3))
    if sums[0] == 0:
        raise ValueError('Undefined ensemble ratio')
    means = tuple(s / count for s in sums)
    return sums[1] / sums[0] - means[2], means, sums


def jk_error(replicates):
    n = D(len(replicates))
    mean = sum(replicates, D(0)) / n
    return (((n-1)/n) * sum(((r-mean)**2 for r in replicates), D(0))).sqrt()


def paired_estimate(rows, configs):
    n = len(configs)
    estimates, means, replicates = {}, {}, {}
    for piter in (4, 40):
        values = [rows[cfg, piter] for cfg in configs]
        estimate, average, totals = estimator(values)
        estimates[piter], means[piter] = estimate, average
        replicates[piter] = []
        for value in values:
            leave = tuple(totals[j]-value[j] for j in range(3))
            if leave[0] == 0:
                raise ValueError('Undefined leave-one-out ratio')
            replicates[piter].append(leave[1]/leave[0]-leave[2]/D(n-1))
    delta = estimates[4]-estimates[40]
    error40 = jk_error(replicates[40])
    paired_error = jk_error([a-b for a, b in zip(replicates[4], replicates[40])])
    return {
        'configurations': n,
        'mean_L_RE_p4': means[4][0], 'mean_L_RE_p40': means[40][0],
        'mean_LS_RE_p4': means[4][1], 'mean_LS_RE_p40': means[40][1],
        'mean_S_RAW': means[40][2],
        'C_p4': estimates[4], 'C_p40': estimates[40],
        'delta_C_p4_minus_p40': delta,
        'C_p4_jackknife_se': jk_error(replicates[4]),
        'C_p40_jackknife_se': error40,
        'delta_C_paired_jackknife_se': paired_error,
        'absolute_relative_shift_percent': abs(delta/estimates[40])*100 if estimates[40] else None,
        'absolute_shift_over_control_se': abs(delta/error40) if error40 else None,
        'control_signal_over_se': abs(estimates[40]/error40) if error40 else None,
    }


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def extremes(rows):
    def pick(field):
        valid = [r for r in rows if r[field] is not None]
        if not valid:
            return None
        best = max(valid, key=lambda r: r[field])
        return {k: float(v) if isinstance(v, D) else v for k, v in best.items()}
    return {
        'records': len(rows),
        'max_relative_shift': pick('absolute_relative_shift_percent'),
        'max_shift_over_control_se': pick('absolute_shift_over_control_se'),
        'max_absolute_shift': {
            k: float(v) if isinstance(v, D) else v
            for k, v in max(rows, key=lambda r: abs(r['delta_C_p4_minus_p40'])).items()
        },
    }


def analyze(base, output):
    sources, groups, windows = read_inputs(base)
    output.mkdir(parents=True, exist_ok=True)
    summary = {'success': True, 'configurations': list(CONFIGS), 'raw_records': 25896,
               'profile_keys': 996, 'combined_profile_keys': 249,
               'estimator': 'mean(LS_RE)/mean(L_RE)-mean(S_RAW), equal configuration weights',
               'uncertainty': 'Matched delete-one-configuration jackknife; no autocorrelation correction',
               'scope': 'Central midpoint98, tau2/4/6, R1/2/4, Nv100; no piter4 acceptance gate',
               'subsets': {}}
    all13 = []
    with localcontext() as context:
        context.prec = 50
        for count in (5, 9, 13):
            configs = CONFIGS[:count]
            records = []
            for key, values in sorted(groups.items()):
                record = dict(zip(KEYS, key))
                record.update(paired_estimate(values, configs))
                if count == 13 or key[1] == -1:
                    records.append(record)
            combined = [r for r in records if r['axis'] == -1]
            write_csv(output / f'ensemble_{count}_combined.csv', combined)
            if count == 13:
                write_csv(output / 'ensemble_13_all_axes.csv', records)
                all13 = records
            summary['subsets'][str(count)] = {
                'configurations': list(configs), 'combined': extremes(combined),
                'combined_by_tau': {str(t): extremes([r for r in combined if r['tau'] == t]) for t in (2, 4, 6)},
            }
            if count == 13:
                summary['subsets']['13']['all_axes'] = extremes(records)
    p4_windows = [r for r in windows if int(r['piter']) == 4]
    p40_windows = [r for r in windows if int(r['piter']) == 40]
    def max_window(items, field):
        row = max(items, key=lambda r: D(r[field]))
        return {'configuration': int(row['configuration']), 'tau': int(row['tau']), 'value': float(row[field])}
    summary['SU3'] = {
        'metric': 'unit: Frobenius norm of Udagger*U-I; det: absolute det(U)-1',
        'p4_max_unit': max_window(p4_windows, 'unit_max'),
        'p4_max_det': max_window(p4_windows, 'det_max'),
        'p40_max_unit': max_window(p40_windows, 'unit_max'),
        'p40_max_det': max_window(p40_windows, 'det_max'),
        'p4_windows_above_1e_minus2': sum(int(r['unit_above_1em2']) > 0 or int(r['det_above_1em2']) > 0 for r in p4_windows),
        'threshold_scope': '1e-2 is a diagnostic count threshold, not an acceptance criterion',
    }
    write_csv(output / 'su3_windows_13.csv', windows)
    (output / 'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'success': True, 'configurations': len(CONFIGS), 'records': len(all13),
                      'combined_by_tau': summary['subsets']['13']['combined_by_tau'],
                      'SU3': summary['SU3']}, indent=2))
    return summary, all13


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    here = Path(__file__).resolve().parent
    parser.add_argument('--inputs', type=Path, default=here/'inputs')
    parser.add_argument('--output', type=Path, default=here/'results')
    args = parser.parse_args()
    analyze(args.inputs, args.output)


if __name__ == '__main__':
    main()
