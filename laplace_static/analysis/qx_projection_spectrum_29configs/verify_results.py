#!/usr/bin/env python3
"""Reconstruct the matched ensemble and verify both preserved baselines."""
import json
import make_plots as analysis


def same_records(expected, stored, label):
    if len(expected) != len(stored) or len(expected) != 996:
        raise ValueError('Changed profile inventory: '+label)
    for before, after in zip(expected, stored):
        if set(before) != set(after) or any(str(value) != after[k] for k, value in before.items()):
            raise ValueError('Ensemble value or paired uncertainty differs: '+label)


def main():
    core, raw, windows, records, sources = analysis.load_data()
    here = analysis.HERE
    stored = core.read_csv(here/'results/profiles_ensemble_29configs.csv')
    same_records(records, stored, '29 configurations')
    for size, selected in (
            (13, [r for r in raw if int(r['configuration']) <= 230]),
            (26, [r for r in raw if int(r['configuration']) not in (780, 790, 800)])):
        reconstructed, _, _ = core.aggregate(selected)
        baseline = core.read_csv(analysis.ROOT/(
            'laplace_static/analysis/qx_projection_spectrum_%dconfigs/results/profiles_ensemble_%dconfigs.csv'
            % (size, size)))
        same_records(reconstructed, baseline, str(size)+'-configuration baseline')
    exported_windows = core.read_csv(here/'results/su3_central_windows_29configs.csv')
    if exported_windows != windows or len(windows) != 348:
        raise ValueError('Changed central-window exports')
    comparisons = analysis.comparisons(core, records, windows)
    stored_comparisons = core.read_csv(here/'results/piter_comparison_29configs.csv')
    if len(comparisons) != len(stored_comparisons) or any(
            {k: str(value) for k, value in expected.items()} != actual
            for expected, actual in zip(comparisons, stored_comparisons)):
        raise ValueError('Changed projection comparison summary')
    expected_summary = analysis.summary_data(raw, windows, records, comparisons, sources)
    if json.loads((here/'results/SUMMARY.json').read_text()) != expected_summary:
        raise ValueError('Summary differs from reconstructed results')
    for name in analysis.PLOTS:
        for extension in ('.png', '.svg'):
            path = here/'results'/(name+extension)
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError('Missing plot: '+name+extension)
    print('PASS: 29 matched configurations, 115536 raw records, 348 window rows; '
          'ensemble values and paired errors reproduced; 13- and 26-configuration baselines unchanged.')


if __name__ == '__main__':
    main()
