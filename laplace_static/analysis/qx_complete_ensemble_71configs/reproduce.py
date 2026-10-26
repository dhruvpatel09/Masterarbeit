#!/usr/bin/env python3
"""Reproduce the matched 71-configuration qX projection, clover and temporal study."""
import argparse
import csv
import hashlib
import importlib.util
import json
import tempfile
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
TARGET = list(range(100, 801, 10))
BASE = [c for c in TARGET if c != 640]
PITERS = (4, 10, 20, 40)
TRACE_PITERS = (4, 40)
TAUS = (2, 4, 6)
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
GEOMETRY = KEYS[1:]
D = Decimal
ANALYSIS = ROOT / 'laplace_static/analysis'
RAW70 = ANALYSIS / 'qx_projection_full_ensemble'
TRACE70 = ANALYSIS / 'qx_clover_trace_full_ensemble'
TIME70 = ANALYSIS / 'qx_temporal_separation_70configs'
BIN70 = ANALYSIS / 'qx_error_binning_70configs'
EXTENSION = HERE / 'inputs/job5946008/campaign'
EXPORTS = ('profiles_raw_ensemble_71configs.csv', 'profiles_raw_and_traceless_71configs.csv',
           'su3_central_windows_71configs.csv', 'piter_comparison_71configs.csv',
           'trace_comparison_71configs.csv', 'temporal_comparison_71configs.csv',
           'temporal_covariance_71configs.csv', 'projection_control_71configs.csv',
           'R4_midpoint_71configs.csv', 'binning_comparison_71configs.csv',
           'R4_midpoint_binning_71configs.csv', 'binning_plans_71configs.json',
           'n640_inclusion_comparison.csv', 'configuration_coverage_71configs.csv',
           'VALIDATION.json', 'SUMMARY.json')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path):
    with Path(path).open(newline='') as stream:
        yield from csv.DictReader(stream)


def write_rows(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + '\n')


def dependencies():
    provenance = json.loads((HERE / 'inputs/SOURCES.json').read_text())
    for name, expected in provenance['repository_references'].items():
        if sha(ROOT / name) != expected:
            raise ValueError('Preserved repository reference changed: ' + name)
    for name, expected in provenance['preserved_inputs'].items():
        if sha(HERE / name) != expected:
            raise ValueError('Preserved n640 input changed: ' + name)
    raw = module('qx71_raw_model', RAW70 / 'make_plots.py')
    trace = module('qx71_clover_model', TRACE70 / 'reproduce.py')
    temporal = module('qx71_temporal_model', TIME70 / 'reproduce.py')
    binning = module('qx71_binning_model', BIN70 / 'reproduce.py')
    return provenance, raw, trace, temporal, binning


def validate_primary(row, configs, piters, geometry, trace=False):
    cfg, p = int(row['configuration']), int(row['piter'])
    key = tuple(int(row[k]) for k in KEYS)
    tau, axis = key[:2]
    if cfg not in configs or p not in piters or key not in geometry:
        raise ValueError('Unexpected configuration, piter or profile geometry')
    if tuple(int(row[k]) for k in ('Nv', 't0', 't1', 'tmid', 'Nsrc')) != (
            100, 98 - tau // 2, 98 + tau // 2, 98, 64 ** 3 * (3 if axis == -1 else 1)):
        raise ValueError('Incorrect mode count, time window or source count')
    if D(row['long_from_midpoint']) != D(row['long']) - D(row['R']) / 2:
        raise ValueError('Incorrect longitudinal coordinate')
    fields = ('L_RE', 'LS_RE', 'S_RAW') + (('T_TRACE', 'LT_RE') if trace else ())
    values = tuple(D(row[k]) for k in fields)
    all_fields = fields + ('L_IM', 'LS_IM') + (('LT_IM',) if trace else ())
    if not all(D(row[k]).is_finite() for k in all_fields) or not values[0]:
        raise ValueError('Nonfinite primary or zero loop denominator')
    if trace and not 0 <= values[3] <= values[2]:
        raise ValueError('Invalid trace component')
    return key, (cfg, p), values


def grouped(paths, configs, piters, geometry, trace=False):
    groups = defaultdict(dict)
    for path, permitted in paths:
        for row in read_rows(path):
            if int(row['piter']) not in permitted:
                continue
            key, ident, values = validate_primary(row, configs, piters, geometry, trace)
            if ident in groups[key]:
                raise ValueError('Duplicate configuration/piter/profile')
            groups[key][ident] = values
    expected = {(c, p) for c in configs for p in piters}
    if set(groups) != geometry:
        raise ValueError('Incomplete profile geometry')
    for key, values in groups.items():
        if set(values) != expected:
            raise ValueError('Incomplete matched configuration cohort')
        for c in configs:
            for j in ((2, 3) if trace else (2,)):
                if len({values[c, p][j] for p in piters}) != 1:
                    raise ValueError('Original-field clover changed across piters')
    return groups


def raw_statistics(values, configs, raw_model):
    with localcontext() as context:
        context.prec = 50
        result = {}
        for p in PITERS:
            total = [sum((values[c, p][j] for c in configs), D(0)) for j in range(3)]
            if not total[0]:
                raise ValueError('Zero ensemble loop denominator')
            result.update({'mean_%s_p%d' % (field, p): total[j] / len(configs)
                           for j, field in enumerate(('L_RE', 'LS_RE', 'S_RAW'))})
            result['C_p%d' % p] = total[1] / total[0] - total[2] / len(configs)
        for p in (4, 10, 20):
            result['delta_C_p%d_minus_p40' % p] = result['C_p%d' % p] - result['C_p40']
        estimates = raw_model.paired_estimator()(values, configs)
        for name in ('C_p4_jackknife_se', 'C_p40_jackknife_se', 'delta_C_paired_jackknife_se'):
            result[name] = estimates[name]
        return result


def check_reference(actual, archived, key_fields, excluded=('configurations',)):
    if len(actual) != len(archived):
        raise ValueError('Archived result row count changed')
    lookup = {tuple(str(r[k]) for k in key_fields): r for r in actual}
    if len(lookup) != len(actual):
        raise ValueError('Duplicate result geometry')
    compared = 0
    for ref in archived:
        key = tuple(str(ref[k]) for k in key_fields)
        row = lookup[key]
        for field, value in ref.items():
            if field in key_fields or field in excluded:
                continue
            if value == '':
                if row[field] != '':
                    raise ValueError('Undefined archived value changed')
            elif D(str(row[field])) != D(value):
                raise ValueError('Archived baseline changed: %s %s' % (key, field))
            compared += 1
    return compared


def compare_piters(records, windows):
    combined = [r for r in records if r['axis'] == -1]
    comparisons = []
    with localcontext() as context:
        context.prec = 50
        for p in PITERS:
            delta = lambda r: r['C_p%d' % p] - r['C_p40']
            absolute = max(combined, key=lambda r: abs(delta(r)))
            relative = max(combined, key=lambda r: abs(delta(r) / r['C_p40']))
            w = [r for r in windows if int(r['piter']) == p]
            maximum = max(w, key=lambda r: D(r['unit_max']))
            determinant = max(w, key=lambda r: D(r['det_max']))
            comparisons.append(dict(piter=p, configurations=len(TARGET),
                combined_positions=len(combined), maximum_absolute_delta=abs(delta(absolute)),
                maximum_relative_delta_percent=100 * abs(delta(relative) / relative['C_p40']),
                maximum_unit_defect=maximum['unit_max'], unit_configuration=maximum['configuration'],
                unit_tau=maximum['tau'], maximum_det_defect=determinant['det_max']))
    return comparisons


def temporal_records(groups, configs, source, temporal, archived_profiles):
    records, covariance, controls = [], [], []
    spatial = sorted({k[1:] for k in groups})
    for geometry in spatial:
        selected = {t: {ident: values for ident, values in groups[(t,) + geometry].items()
                        if ident[0] in configs} for t in TAUS}
        point, loo = temporal.paired_estimates(selected, configs, source)
        for definition in ('RAW', 'TL'):
            for p in TRACE_PITERS:
                stats, cov = temporal.channel_statistics(point, loo, definition, p)
                row = dict(zip(GEOMETRY, geometry))
                row.update(definition=definition, piter=p, configurations=len(configs))
                row.update(stats)
                records.append(row)
                field = 'C_%s_p%d' % (definition, p)
                for t in TAUS:
                    ref = archived_profiles[(t,) + geometry]
                    if stats['C_tau%d' % t] != ref[field] or stats['C_tau%d_jackknife_se' % t] != ref[field + '_jackknife_se']:
                        raise ValueError('Temporal marginal differs from profile estimate')
                with localcontext() as context:
                    context.prec = 50
                    for i, t in enumerate(TAUS):
                        for j, u in enumerate(TAUS):
                            product = cov[i][i] * cov[j][j]
                            rho = cov[i][j] / product.sqrt() if product > 0 else ''
                            out = dict(zip(GEOMETRY, geometry))
                            out.update(definition=definition, piter=p, tau_i=t, tau_j=u,
                                       covariance=cov[i][j], correlation=rho)
                            covariance.append(out)
            for control in temporal.projection_controls(point, loo, definition):
                row = dict(zip(GEOMETRY, geometry))
                row.update(definition=definition, configurations=len(configs))
                row.update(control)
                controls.append(row)
    return records, covariance, controls


def binned_records(groups, source, temporal, binning, full):
    plans = [binning.make_plan(TARGET, b, 10) for b in (1, 2, 4, 8)]
    reference = {tuple(r[k] for k in GEOMETRY) + (r['definition'], r['piter']): r for r in full}
    records = []
    for geometry in sorted({k[1:] for k in groups}):
        values = {t: groups[(t,) + geometry] for t in TAUS}
        for plan in plans:
            point, loo = binning.block_estimates(values, plan['blocks'], source)
            unblocked_point, unblocked_loo = binning.block_estimates(
                values, [[c] for c in plan['selected_configurations']], source)
            if point != unblocked_point:
                raise ValueError('Blocking changed the selected-cohort point estimate')
            for definition in ('RAW', 'TL'):
                blocked, _ = temporal.channel_statistics(point, loo, definition, 40)
                unblocked, _ = temporal.channel_statistics(unblocked_point, unblocked_loo, definition, 40)
                ref = reference[geometry + (definition, 40)]
                row = dict(zip(GEOMETRY, geometry))
                row.update(definition=definition, piter=40, bin_width=plan['bin_width'],
                           configurations=plan['configurations'], bins=plan['bins'])
                with localcontext() as context:
                    context.prec = 50
                    for field in binning.FIELDS:
                        if blocked[field] != unblocked[field]:
                            raise ValueError('Blocked and unblocked points differ')
                        se, se0 = blocked[field + '_jackknife_se'], unblocked[field + '_jackknife_se']
                        if plan['bin_width'] == 1 and (blocked[field] != ref[field] or se0 != ref[field + '_jackknife_se']):
                            raise ValueError('Binning baseline differs from full71 temporal result')
                        row[field] = blocked[field]
                        row[field + '_blocked_se'] = se
                        row[field + '_unblocked_same_cohort_se'] = se0
                        row[field + '_se_ratio'] = se / se0 if se0 else ''
                        row[field + '_shift_from_full71'] = blocked[field] - ref[field]
                records.append(row)
    return records, plans


def numerical_exports(output):
    provenance, raw_model, source, temporal, binning = dependencies()
    review = module('qx71_review_replay', HERE / 'verify_review.py').verify()
    geometry = source.expected_geometry()
    groups = grouped([(RAW70 / 'results/profiles_raw_70configs.csv', PITERS),
                      (EXTENSION / 'profiles_raw_and_trace.csv', PITERS)], TARGET, PITERS, geometry)
    traces = grouped([(TRACE70 / 'inputs/profiles_raw_and_trace_70configs.csv', TRACE_PITERS),
                      (EXTENSION / 'profiles_raw_and_trace.csv', TRACE_PITERS)], TARGET, TRACE_PITERS, geometry, trace=True)
    for key in geometry:
        for ident, values in traces[key].items():
            if values[:3] != groups[key][ident]:
                raise ValueError('Raw and trace campaigns disagree on original primaries')
    windows = list(read_rows(RAW70 / 'results/su3_central_windows_70configs.csv')) + list(read_rows(EXTENSION / 'su3_windows.csv'))
    windows.sort(key=lambda r: (int(r['configuration']), int(r['tau']), int(r['piter'])))
    raw_model.core.validate_windows(windows, TARGET, PITERS)
    if len(windows) != 852:
        raise ValueError('Incomplete full71 central windows')
    records, trace_records, raw_base, trace_base = [], [], [], []
    for key in sorted(geometry):
        row = dict(zip(KEYS, key)); row['configurations'] = len(TARGET)
        row.update(raw_statistics(groups[key], TARGET, raw_model)); records.append(row)
        old = dict(zip(KEYS, key)); old['configurations'] = len(BASE)
        old.update(raw_statistics(groups[key], BASE, raw_model)); raw_base.append(old)
        row = dict(zip(KEYS, key)); row['configurations'] = len(TARGET)
        row.update(source.paired_statistics(traces[key], TARGET)); trace_records.append(row)
        old = dict(zip(KEYS, key)); old['configurations'] = len(BASE)
        old.update(source.paired_statistics(traces[key], BASE)); trace_base.append(old)
        for p in TRACE_PITERS:
            if row['C_RAW_p%d' % p] != records[-1]['C_p%d' % p]:
                raise ValueError('RAW estimator differs between matched analyses')
    raw_fields = check_reference(raw_base, list(read_rows(RAW70 / 'results/profiles_ensemble_70configs.csv')), KEYS)
    trace_fields = check_reference(trace_base, list(read_rows(TRACE70 / 'results/profiles_ensemble_70configs.csv')), KEYS)
    full, covariance, controls = temporal_records(traces, TARGET, source, temporal,
        {tuple(r[k] for k in KEYS): r for r in trace_records})
    # Replay the earlier temporal result with the original cohort as a separate regression.
    previous_time, _, _ = temporal_records(traces, BASE, source, temporal,
        {tuple(r[k] for k in KEYS): r for r in trace_base})
    time_fields = check_reference(previous_time, list(read_rows(TIME70 / 'results/temporal_comparison_70configs.csv')),
                                  GEOMETRY + ('definition', 'piter'))
    binned, plans = binned_records(traces, source, temporal, binning, full)
    piter_comparison = compare_piters(records, windows)
    combined = [r for r in trace_records if r['axis'] == -1]
    trace_comparison = []
    with localcontext() as context:
        context.prec = 50
        for p in TRACE_PITERS:
            fraction = 'trace_fraction_percent_p%d' % p
            trace_comparison.append(dict(piter=p, combined_positions=len(combined),
                minimum_trace_fraction_percent=min(r[fraction] for r in combined),
                maximum_trace_fraction_percent=max(r[fraction] for r in combined)))
    influence = []
    old_lookup = {tuple(r[k] for k in KEYS): r for r in trace_base}
    with localcontext() as context:
        context.prec = 50
        for row in combined:
            old = old_lookup[tuple(row[k] for k in KEYS)]
            result = {k: row[k] for k in KEYS}
            for definition in ('RAW', 'TL'):
                field = 'C_%s_p40' % definition
                se = row[field + '_jackknife_se']
                result[field + '_70configs'] = old[field]
                result[field + '_71configs'] = row[field]
                result[field + '_71configs_jackknife_se'] = se
                result[field + '_including_minus_excluding_n640'] = row[field] - old[field]
                result[field + '_change_over_full71_se'] = (row[field] - old[field]) / se if se else ''
            influence.append(result)
    midpoint = [r for r in full if tuple(r[k] for k in GEOMETRY) == (-1, 4, 2, 0, 0)]
    midpoint_bins = [r for r in binned if tuple(r[k] for k in GEOMETRY) == (-1, 4, 2, 0, 0)]
    coverage = [dict(configuration=c, taus='2;4;6', raw_piters='4;10;20;40',
                     traceless_matched_piters='4;40', raw_rows=3984, trace_rows=1992,
                     window_rows=12, recovered_configuration=(c == 640)) for c in TARGET]
    validation = dict(success=True, n640_review=review, configurations=TARGET, configuration_count=len(TARGET), cases=213,
        assembled_raw_primary_rows=sum(len(g) for g in groups.values()),
        matched_trace_primary_rows=sum(len(g) for g in traces.values()), window_rows=len(windows),
        baseline70=dict(raw_fields_exact=raw_fields, trace_fields_exact=trace_fields,
                        temporal_fields_exact=time_fields),
        temporal_marginals_reproduced_exactly=True, raw_trace_regression_exact=True)
    summary = dict(success=True, base_commit=provenance['base_commit'], configurations=TARGET,
        configuration_count=71, excluded_configurations=[], cases=213,
        raw_piters=list(PITERS), traceless_matched_piters=list(TRACE_PITERS), taus=list(TAUS),
        raw_primary_rows=validation['assembled_raw_primary_rows'],
        trace_primary_rows=validation['matched_trace_primary_rows'], central_window_rows=len(windows),
        profile_geometry_rows=len(records), combined_geometry_rows=len(combined),
        temporal_rows=len(full), covariance_entries=len(covariance), projection_controls=len(controls),
        binning_rows=len(binned), piter_comparison=piter_comparison, trace_comparison=trace_comparison,
        R4_midpoint=midpoint, R4_midpoint_binning=midpoint_bins,
        maximum_abs_n640_inclusion_change_over_full71_se={definition:max(abs(r['C_%s_p40_change_over_full71_se' % definition]) for r in influence)
            for definition in ('RAW', 'TL')},
        binning=[{k:p[k] for k in ('bin_width', 'configurations', 'bins', 'dropped_incomplete_bin_configurations', 'gaps')} for p in plans],
        estimator='Equal-configuration ratio of means; 50 decimal digits applied to preserved primaries; no bias correction.',
        convention='Original-field raw six-plane clover; traceless side observable removes sum |Tr F|^2/3; no additional factor 1/2.',
        statistics='Paired delete-one-configuration jackknife; temporal differences and trace fractions recomputed within each synchronized omission.',
        monte_carlo_metadata={k:provenance['monte_carlo_metadata'][k] for k in provenance['monte_carlo_metadata']},
        limitations=[
            'All main full71 profile and temporal estimates include n640. Older 70-configuration exports are unchanged.',
            'Traceless piter10/20 at n640 alone do not supply matched 71-configuration traceless piter10/20 data.',
            'The selected four-configuration piter6/8/80 study remains a separate diagnostic sample.',
            'Errors are preliminary without a verified Monte Carlo-time assignment or autocorrelation correction.',
            'Binning uses conditional label order. Only width1 retains all71; complete equal-size widths2/4/8 retain70/68/64. Dropped labels and same-cohort error denominators are explicit.',
            'With n640 restored, bins after the former gap are regrouped. Changes from older binning plots cannot be attributed solely to one added sample.',
            'No temporal plateau or ground-state saturation is established; no renormalized physical-density interpretation is made.',
            'The correlated 70/71 inclusion change divided by SE71 is descriptive, not an independent-sample significance test.',
            'Windows, probes and spatial axes are correlated measurements, not additional independent configurations.'])
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    tables = (records, trace_records, windows, piter_comparison, trace_comparison, full,
              covariance, controls, midpoint, binned, midpoint_bins)
    for name, table in zip(EXPORTS[:11], tables):
        write_rows(output / name, table)
    write_json(output / EXPORTS[11], plans)
    write_rows(output / EXPORTS[12], influence)
    write_rows(output / EXPORTS[13], coverage)
    write_json(output / EXPORTS[14], validation)
    write_json(output / EXPORTS[15], summary)
    return records, trace_records, windows, full, covariance, controls, binned, plans, influence, summary


def verify():
    with tempfile.TemporaryDirectory() as directory:
        result = numerical_exports(directory)
        for name in EXPORTS:
            if (Path(directory) / name).read_bytes() != (RESULTS / name).read_bytes():
                raise ValueError('Numerical export differs: ' + name)
    manifest = json.loads((RESULTS / 'FIGURES.json').read_text())
    for name, expected in manifest.items():
        if sha(RESULTS / name) != expected:
            raise ValueError('Preserved figure changed: ' + name)
    print('PASS full71: 71 configurations; 213 cases; 282864 raw and 141432 matched trace primaries; 852 windows; 12 n640 stdout logs replayed; all 70-configuration profile and temporal baselines exact; all exports reproduced.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--plots', action='store_true')
    parser.add_argument('--output', type=Path, default=RESULTS)
    args = parser.parse_args()
    if args.verify:
        if args.plots or args.output != RESULTS:
            parser.error('--verify reads the preserved outputs without rewriting them')
        verify()
    else:
        result = numerical_exports(args.output)
        if args.plots:
            module('qx71_renderer', HERE / 'render_plots.py').render(*result[:-1], args.output)
            write_json(args.output / 'FIGURES.json', {p.name:sha(p) for p in sorted(args.output.iterdir())
                if p.suffix in ('.png', '.svg')})
        print('PASS complete 71-configuration numerical exports')


if __name__ == '__main__':
    main()
