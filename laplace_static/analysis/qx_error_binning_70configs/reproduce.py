#!/usr/bin/env python3
"""Conditional error-binning sensitivity in qX configuration-label order."""
import argparse
import hashlib
import importlib.util
import json
import tempfile
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
TEMPORAL = HERE.parent / 'qx_temporal_separation_70configs'
D = Decimal
FIELDS = ('C_tau2', 'C_tau4', 'C_tau6', 'delta_tau4_minus_tau2',
          'delta_tau6_minus_tau4', 'delta_tau6_minus_tau2')
EXPORTS = ('binning_comparison_70configs.csv', 'R4_midpoint_binning.csv',
           'binning_plans.json', 'SUMMARY.json')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def load_dependencies():
    provenance = json.loads((HERE / 'inputs/SOURCES.json').read_text())
    for name, digest in provenance['repository_references'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Preserved temporal analysis changed: ' + name)
    spec = importlib.util.spec_from_file_location('qx_temporal_binning_source', TEMPORAL / 'reproduce.py')
    temporal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(temporal)
    source, _ = temporal.verified_upstream()
    return temporal, source, provenance


def make_plan(configs, width, label_increment=10):
    if isinstance(width, bool) or not isinstance(width, int) or width < 1 or configs != sorted(set(configs)) or not configs:
        raise ValueError('Need sorted unique configuration labels and positive integer width')
    runs = [[configs[0]]]
    for c in configs[1:]:
        if c - runs[-1][-1] != label_increment:
            runs.append([])
        runs[-1].append(c)
    blocks, dropped = [], []
    for run in runs:
        complete = len(run) // width * width
        blocks.extend([run[i:i + width] for i in range(0, complete, width)])
        dropped.extend(run[complete:])
    if len(blocks) < 2:
        raise ValueError('Need at least two complete equal-size bins')
    selected = [c for block in blocks for c in block]
    return dict(bin_width=width, blocks=blocks, bins=len(blocks), configurations=len(selected),
                selected_configurations=selected, dropped_incomplete_bin_configurations=dropped,
                contiguous_label_runs=runs,
                gaps=[dict(left=runs[i][-1], right=runs[i + 1][0]) for i in range(len(runs) - 1)])


def block_estimates(values, blocks, source):
    """Ratio-of-means point and synchronized delete-one-equal-bin replicates."""
    if len(blocks) < 2 or not blocks[0] or any(len(b) != len(blocks[0]) for b in blocks):
        raise ValueError('Need at least two nonempty equal-size bins')
    selected = [c for block in blocks for c in block]
    if len(set(selected)) != len(selected):
        raise ValueError('Configuration appears in more than one bin')
    if set(values) != {2, 4, 6} or any((c, p) not in values[t]
            for t in (2, 4, 6) for c in selected for p in (4, 40)):
        raise ValueError('Incomplete matched temporal or projection cohort')
    n, width = len(selected), len(blocks[0])
    with localcontext() as context:
        context.prec = 50
        totals = {t: {p: tuple(sum((values[t][c, p][j] for c in selected), D(0)) for j in range(5))
                      for p in (4, 40)} for t in (2, 4, 6)}
        point = {t: source.observables(totals[t], n) for t in (2, 4, 6)}
        loo = []
        for block in blocks:
            sample = {}
            for t in (2, 4, 6):
                remaining = {p: tuple(totals[t][p][j] - sum((values[t][c, p][j] for c in block), D(0))
                                      for j in range(5)) for p in (4, 40)}
                sample[t] = source.observables(remaining, n - width)
            loo.append(sample)
        return point, loo


def numerical_exports(output):
    temporal, source, provenance = load_dependencies()
    groups = source.group_primaries(source.read_csv(source.HERE / 'inputs/profiles_raw_and_trace_70configs.csv'))
    archived = {tuple(int(r[k]) for k in temporal.GEOMETRY) + (r['definition'], int(r['piter'])): r
                for r in source.read_csv(TEMPORAL / 'results/temporal_comparison_70configs.csv')}
    plans = [make_plan(source.TARGET, b, provenance['label_increment']) for b in provenance['bin_widths']]
    records = []
    with localcontext() as context:
        context.prec = 50
        for geometry in sorted({k[1:] for k in groups}):
            values = {t: groups[(t,) + geometry] for t in (2, 4, 6)}
            for plan in plans:
                point, loo = block_estimates(values, plan['blocks'], source)
                # Use identical selected samples for the blocked/unblocked error comparison.
                unblocked_point, unblocked_loo = block_estimates(
                    values, [[c] for c in plan['selected_configurations']], source)
                if point != unblocked_point:
                    raise ValueError('Blocking changed the ratio-of-means point on the same cohort')
                for definition in ('RAW', 'TL'):
                    blocked, _ = temporal.channel_statistics(point, loo, definition, 40)
                    unblocked, _ = temporal.channel_statistics(unblocked_point, unblocked_loo, definition, 40)
                    full = archived[geometry + (definition, 40)]
                    row = dict(zip(temporal.GEOMETRY, geometry))
                    row.update(definition=definition, piter=40, bin_width=plan['bin_width'],
                               configurations=plan['configurations'], bins=plan['bins'])
                    for field in FIELDS:
                        point_value = blocked[field]
                        se, se0 = blocked[field + '_jackknife_se'], unblocked[field + '_jackknife_se']
                        if point_value != unblocked[field]:
                            raise ValueError('Blocked/unblocked points disagree')
                        if plan['configurations'] == 70 and (point_value != D(full[field]) or
                                se0 != D(full[field + '_jackknife_se'])):
                            raise ValueError('Full70 baseline disagrees with archived paired errors')
                        row[field] = point_value
                        row[field + '_blocked_se'] = se
                        row[field + '_unblocked_same_cohort_se'] = se0
                        row[field + '_se_ratio'] = se / se0 if se0 else ''
                        row[field + '_shift_from_full70'] = point_value - D(full[field])
                    records.append(row)
        midpoint = [r for r in records if tuple(r[k] for k in temporal.GEOMETRY) == (-1, 4, 2, 0, 0)]
        summaries = []
        for plan in plans:
            for definition in ('RAW', 'TL'):
                candidates = [r for r in records if r['axis'] == -1 and r['definition'] == definition
                              and r['bin_width'] == plan['bin_width']]
                for field in FIELDS:
                    ratios = [r[field + '_se_ratio'] for r in candidates if r[field + '_se_ratio'] != '']
                    summaries.append(dict(bin_width=plan['bin_width'], definition=definition, observable=field,
                        combined_positions=len(candidates), minimum_error_ratio=min(ratios),
                        maximum_error_ratio=max(ratios),
                        maximum_absolute_point_shift_from_full70=max(abs(r[field + '_shift_from_full70']) for r in candidates)))
        summary = dict(success=True, base_commit=provenance['base_commit'],
            purpose='Conditional binning sensitivity in increasing configuration-label order',
            primary_cohort=source.TARGET, configurations=70, excluded_configuration=640,
            piters=[40], definitions=['RAW', 'TL'], taus=[2, 4, 6],
            spatial_geometries=332, combined_spatial_geometries=83, comparison_rows=len(records),
            bins=[{k: p[k] for k in ('bin_width', 'bins', 'configurations', 'dropped_incomplete_bin_configurations', 'gaps')}
                  for p in plans],
            archived_full70_values_and_unblocked_errors_reproduced_exactly=True,
            method=dict(point='Ratio of raw configuration means; same selected configuration weights in blocked and unblocked calculations',
                blocked_error='Delete one equal-size bin, recompute all ratios, then take temporal contrasts within the same replicate. Jackknife factor (B-1)/B for B bins.',
                error_ratio='Blocked SE / unblocked SE on exactly the same selected configurations',
                bin_width_unit='Number of consecutive available configurations with nominal label increment 10',
                gap_handling=provenance['bin_policy'],
                point_shifts='Diagnostic larger-bin subsets can change the point estimate; each shift from the full70 result is exported.',
                uncertainty_status='Conditional sensitivity study. No autocorrelation time, final corrected errors or significance tests are inferred.'),
            ordering_assumption=provenance['ordering_assumption'],
            unresolved_metadata={k: provenance[k] for k in ('MDU_per_label_increment', 'replica_or_stream_assignment', 'thermalization_cut', 'p2_meaning')},
            reference_analysis_MCconversionFactor=200,
            reference_analysis_MCconversionFactor_applied_to_profile_labels=False,
            R4_midpoint=midpoint, combined_axis_error_ratio_ranges=summaries,
            limitations=[
                'Chronological single-stream ordering is a working assumption. Numerical label order is used without an MDU conversion.',
                'Binning does not prove that blocks or configurations are independent; a stable-looking curve on 70 samples is not a reliable upper bound on long autocorrelation tails.',
                'The source r0/r1 autocorrelation parameters are observable-dependent and are not reused as flux-profile tail estimates.',
                'Bins of width 1 and 2 retain the full70 cohort. Widths 4 and 8 discard recorded incomplete end bins and remain subset diagnostics.',
                'A gap starts a new bin-construction segment, not a statistically independent replica.',
                'Only piter40 controls are needed here; prior piter comparisons remain in the projection and temporal-separation studies.',
                'All three temporal differences are correlated, and the nonadjacent difference is the sum of the adjacent differences.',
                'Spatial probes and axes are correlated and do not increase the configuration sample count.'])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    source.write_csv(output / EXPORTS[0], records)
    source.write_csv(output / EXPORTS[1], midpoint)
    (output / EXPORTS[2]).write_text(json.dumps(plans, indent=2) + '\n')
    (output / EXPORTS[3]).write_text(json.dumps(summary, indent=2, default=str, allow_nan=False) + '\n')
    return records, plans, summary


def verify():
    with tempfile.TemporaryDirectory() as directory:
        numerical_exports(directory)
        for name in EXPORTS:
            if (Path(directory) / name).read_bytes() != (RESULTS / name).read_bytes():
                raise ValueError('Numerical export differs: ' + name)
    print('PASS conditional binning: all full70 baseline values and paired errors exact; widths 1/2 retain 70 configs; widths 4/8 retain 68/64; 2656 geometry/definition/bin comparisons reproduced.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    if args.verify:
        if args.plots:
            parser.error('--verify does not rewrite plots')
        verify()
    else:
        records, plans, _ = numerical_exports(RESULTS)
        if args.plots:
            import render_plots
            render_plots.render(records, plans, RESULTS)
        print('PASS conditional qX binning exports')


if __name__ == '__main__':
    main()
