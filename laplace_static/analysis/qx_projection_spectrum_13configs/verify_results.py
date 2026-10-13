#!/usr/bin/env python3
"""Check the complete scan inventory and reproduce the archived p4/p40 baseline."""
import csv
import json
from pathlib import Path
import make_plots as analysis


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    analysis.verify_inputs()
    here = analysis.HERE
    raw = rows(here/'inputs/job5926668_profiles_raw.csv')
    expected_cases = {(c, t, p) for c in range(110, 231, 10)
                      for t in (2, 4, 6) for p in (10, 20)}
    counts = {}
    seen = set()
    for row in raw:
        case = tuple(int(row[k]) for k in ('configuration', 'tau', 'piter'))
        key = case+tuple(int(row[k]) for k in ('axis', 'R', 'long', 'tr1', 'tr2'))
        if key in seen:
            raise ValueError('Duplicate imported measurement')
        seen.add(key)
        counts[case] = counts.get(case, 0)+1
    if set(counts) != expected_cases or set(counts.values()) != {332} or len(raw) != 25896:
        raise ValueError('Incomplete imported scan')
    imported_windows = rows(here/'inputs/job5926668_su3_windows.csv')
    if len(imported_windows) != 78 or {
        tuple(int(row[k]) for k in ('configuration', 'tau', 'piter'))
        for row in imported_windows} != expected_cases:
        raise ValueError('Incomplete imported central-window scan')
    baseline = rows(analysis.ROOT/'laplace_static/analysis/qx_hyp_projection/results/profiles_ensemble_13configs.csv')
    current = rows(here/'results/profiles_ensemble_13configs.csv')
    if len(current) != 996 or len(baseline) != len(current):
        raise ValueError('Changed profile-position inventory')
    for before, after in zip(baseline, current):
        if any(after[k] != before[k] for k in before):
            raise ValueError('Archived p4/p40 result or paired jackknife changed')
        if any('C_p'+str(p) not in after for p in analysis.PITERS):
            raise ValueError('Missing measured prescription')
    summary = json.loads((here/'results/SUMMARY.json').read_text())
    if summary['configurations'] != list(range(110, 231, 10)) or summary['profile_piters'] != list(analysis.PITERS):
        raise ValueError('Changed summary cohort')
    windows = rows(here/'results/su3_central_windows_13configs.csv')
    expected_all = {(c, t, p) for c in range(110, 231, 10)
                    for t in (2, 4, 6) for p in analysis.PITERS}
    if len(windows) != 156 or {
        tuple(int(r[k]) for k in ('configuration', 'tau', 'piter'))
        for r in windows} != expected_all:
        raise ValueError('Window/profile prescription inventory differs')
    for name in analysis.KEEP:
        if not (here/'results'/name).is_file():
            raise ValueError('Missing numerical export or plot: '+name)
    print('PASS: 39 cases, 78 launches, complete four-piter cohort; '
          'archived p4/p40 profiles and paired errors reproduced exactly.')


if __name__ == '__main__':
    main()

