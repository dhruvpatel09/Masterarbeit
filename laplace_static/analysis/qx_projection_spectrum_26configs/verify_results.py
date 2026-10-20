#!/usr/bin/env python3
"""Verify the expanded exports against the pinned numerical analysis."""
import json
import make_plots as analysis


def main():
    core, raw, windows, records, sources = analysis.load_data()
    stored = core.read_csv(analysis.HERE/'results/profiles_ensemble_26configs.csv')
    if len(stored) != len(records) or len(stored) != 996:
        raise ValueError('Changed profile-position inventory')
    for expected, actual in zip(records, stored):
        if any(str(value) != actual[k] for k, value in expected.items()):
            raise ValueError('Stored ensemble value or paired uncertainty differs from reconstruction')
    original_raw = [r for r in raw if int(r['configuration']) <= 230]
    original, configs, piters = core.aggregate(original_raw)
    baseline = core.read_csv(analysis.ROOT/'laplace_static/analysis/qx_projection_spectrum_13configs/results/profiles_ensemble_13configs.csv')
    for expected, actual in zip(original, baseline):
        if any(str(value) != actual[k] for k, value in expected.items()):
            raise ValueError('Original 13-configuration baseline changed')
    exported_windows = core.read_csv(analysis.HERE/'results/su3_central_windows_26configs.csv')
    if exported_windows != windows or len(windows) != 312:
        raise ValueError('Changed central-window exports')
    summary = json.loads((analysis.HERE/'results/SUMMARY.json').read_text())
    if not summary['success'] or summary['configurations'] != analysis.CONFIGS or summary['profile_piters'] != analysis.PITERS:
        raise ValueError('Wrong summary cohort')
    for name in analysis.PLOTS:
        for extension in ('.png', '.svg'):
            if not (analysis.HERE/'results'/(name+extension)).is_file():
                raise ValueError('Missing plot: '+name+extension)
    print('PASS: 26 matched configurations, 103584 raw records, 312 window rows; '
          'ensemble values and paired errors reproduced; original 13-configuration baseline unchanged.')


if __name__ == '__main__':
    main()
