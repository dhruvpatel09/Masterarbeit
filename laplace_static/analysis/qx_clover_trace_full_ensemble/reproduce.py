#!/usr/bin/env python3
"""Reproduce the matched 70-configuration raw and traceless clover comparison."""
import argparse
import csv
import hashlib
import json
import tempfile
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE/'results'
TARGET = [c for c in range(100, 801, 10) if c != 640]
PITERS = (4, 40)
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
IDENT = ('configuration', 'tau', 'piter', 'axis', 'R', 'long', 'tr1', 'tr2')
FIELDS = ('L_RE', 'LS_RE', 'S_RAW', 'T_TRACE', 'LT_RE')
D = Decimal


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def expected_geometry():
    offsets = {(r, l, x, z) for r in (1, 2, 4)
               for l, x, z in ({(0, 1, 2)} | {(l, x, 0) for l in range(-1, r+2) for x in range(-2, 3)})}
    return {(tau, axis, r, l, x, z) for tau in (2, 4, 6) for axis in (-1, 0, 1, 2)
            for r, l, x, z in offsets}


def group_primaries(rows):
    groups = defaultdict(dict)
    geometry = expected_geometry()
    for row in rows:
        cfg, p = int(row['configuration']), int(row['piter'])
        key = tuple(int(row[k]) for k in KEYS)
        tau, axis = key[:2]
        if cfg not in TARGET or p not in PITERS or key not in geometry:
            raise ValueError('Unexpected configuration, prescription or geometry')
        if (int(row['Nv']), int(row['t0']), int(row['t1']), int(row['tmid']), int(row['Nsrc'])) != (
                100, 98-tau//2, 98+tau//2, 98, 64**3*(3 if axis == -1 else 1)):
            raise ValueError('Incorrect measurement geometry or source count')
        if D(row['long_from_midpoint']) != D(row['long'])-D(row['R'])/2:
            raise ValueError('Incorrect longitudinal coordinate')
        if (cfg, p) in groups[key]:
            raise ValueError('Duplicate configuration/prescription/geometry')
        all_fields = [D(row[k]) for k in FIELDS+('L_IM', 'LS_IM', 'LT_IM')]
        values = tuple(D(row[k]) for k in FIELDS)
        if not all(v.is_finite() for v in all_fields) or values[0] == 0:
            raise ValueError('Nonfinite primary or zero loop denominator')
        if not 0 <= values[3] <= values[2]:
            raise ValueError('Invalid nonnegative trace component')
        groups[key][cfg, p] = values
    if set(groups) != geometry or len(rows) != 139440:
        raise ValueError('Incomplete 70-configuration profile geometry')
    expected = {(c, p) for c in TARGET for p in PITERS}
    for key, values in groups.items():
        if set(values) != expected:
            raise ValueError('Incomplete matched cohort at '+str(key))
        for c in TARGET:
            if any(values[c, 4][j] != values[c, 40][j] for j in (2, 3)):
                raise ValueError('Original clover or trace changed between prescriptions')
    return groups


def observables(sums, count):
    result = {}
    for p in PITERS:
        L, LS, S, T, LT = sums[p]
        if L == 0:
            raise ValueError('Undefined full or jackknife loop ratio')
        raw, trace = LS/L-S/D(count), LT/L-T/D(count)
        result['C_RAW_p%d' % p] = raw
        result['C_TRACE_p%d' % p] = trace
        result['C_TL_p%d' % p] = (LS-LT)/L-(S-T)/D(count)
        result['trace_fraction_percent_p%d' % p] = 100*trace/raw if raw else None
    for field in ('RAW', 'TRACE', 'TL'):
        result['delta_C_%s_p4_minus_p40' % field] = result['C_%s_p4' % field]-result['C_%s_p40' % field]
    return result


def paired_statistics(values, configs):
    """Form ratios of means, then omit the same configuration in all observables."""
    n = len(configs)
    if n < 2:
        raise ValueError('Need at least two configurations for jackknife errors')
    with localcontext() as context:
        context.prec = 50
        totals = {p: tuple(sum((values[c, p][j] for c in configs), D(0)) for j in range(5))
                  for p in PITERS}
        point = observables(totals, n)
        loo = []
        for c in configs:
            omitted = {p: tuple(totals[p][j]-values[c, p][j] for j in range(5)) for p in PITERS}
            loo.append(observables(omitted, n-1))
        result = {k: '' if v is None else v for k, v in point.items()}
        for p in PITERS:
            for j, field in enumerate(FIELDS):
                result['mean_%s_p%d' % (field, p)] = totals[p][j]/D(n)
        for field in point:
            samples = [r[field] for r in loo]
            if point[field] is None or any(v is None for v in samples):
                result[field+'_jackknife_se'] = ''
                continue
            center = sum(samples, D(0))/D(n)
            variance = D(n-1)/D(n)*sum(((v-center)**2 for v in samples), D(0))
            result[field+'_jackknife_se'] = variance.sqrt()
        return result


def verify_inputs():
    sources = json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name, digest in sources['files'].items():
        if sha(HERE/name) != digest:
            raise ValueError('Preserved input changed: '+name)
    for name, digest in sources['repository_references'].items():
        if sha(ROOT/name) != digest:
            raise ValueError('Preserved reference changed: '+name)
    report = json.loads((HERE/'inputs/VALIDATION.json').read_text())
    if (report['success'], report['configurations'], report['cases'], report['new_cases'],
            report['mpi_launches'], report['new_mpi_launches'], report['raw_records'], report['window_records'],
            report['raw_regression_exact'], D(report['raw_regression_max_absolute_difference']),
            report['central_windows_exactly_reproduced']) != (
                True, TARGET, 210, 209, 420, 418, 139440, 420, True, D(0), True):
        raise ValueError('Incomplete or inconsistent archived validation')
    if {(r['configuration'], r['tau']) for r in report['checks']} != {(c,t) for c in TARGET for t in (2,4,6)}:
        raise ValueError('Incomplete validation case inventory')
    if not all(r['success'] and r['raw_records'] == 664 and r['window_records'] == 2 for r in report['checks']):
        raise ValueError('Failed case validation')
    return sources, report


def verify_raw_regression(rows, sources):
    current = {tuple(int(r[k]) for k in IDENT): r for r in rows}
    name = sources['raw_reference_table']
    count = 0
    with (ROOT/name).open(newline='') as stream:
        for ref in csv.DictReader(stream):
            if int(ref['piter']) not in PITERS:
                continue
            key = tuple(int(ref[k]) for k in IDENT)
            if key not in current:
                raise ValueError('Missing preserved raw geometry')
            row = current[key]
            for field in ('L_RE', 'L_IM', 'S_RAW', 'LS_RE', 'LS_IM'):
                if D(row[field]) != D(ref[field]):
                    raise ValueError('Raw regression changed: '+field+' '+str(key))
            count += 1
    if count != 139440:
        raise ValueError('Incomplete raw regression')


def validate_windows(windows):
    keys = [(int(r['configuration']),int(r['tau']),int(r['piter'])) for r in windows]
    if len(keys) != 420 or set(keys) != {(c,t,p) for c in TARGET for t in (2,4,6) for p in PITERS}:
        raise ValueError('Incomplete or duplicate central windows')
    for row in windows:
        tau, p = int(row['tau']), int(row['piter'])
        if (int(row['t0']),int(row['t1']),int(row['count'])) != (98-tau//2,98+tau//2,tau*64**3):
            raise ValueError('Incorrect central window metadata')
        for prefix in ('unit','det'):
            maximum, rms = D(row[prefix+'_max']), D(row[prefix+'_rms'])
            if not maximum.is_finite() or not rms.is_finite() or not 0 <= rms <= maximum:
                raise ValueError('Invalid central-window defect')
            if p == 40 and maximum > D('1e-12'):
                raise ValueError('High-iteration SU(3) control failed')


def position(row):
    return {k:int(row[k]) for k in KEYS}


def summarize(records, report):
    combined = [r for r in records if r['axis'] == -1]
    comparisons = []
    for p in PITERS:
        fraction = 'trace_fraction_percent_p%d' % p
        trace = 'C_TRACE_p%d' % p
        most = max(combined, key=lambda r: abs(r[fraction]))
        maximum = max(combined, key=lambda r: abs(r[trace]))
        comparisons.append(dict(piter=p,combined_positions=len(combined),
             minimum_trace_fraction_percent=str(min(r[fraction] for r in combined)),
             maximum_trace_fraction_percent=str(max(r[fraction] for r in combined)),
             largest_absolute_trace_fraction_position=position(most),
             maximum_absolute_connected_trace_component=str(abs(maximum[trace])),
             largest_absolute_trace_component_position=position(maximum),
             all_raw_trace_and_traceless_connected_profiles_negative=all(
                 r['C_%s_p%d' % (kind,p)] < 0 for kind in ('RAW','TRACE','TL') for r in combined)))
    projection = []
    with localcontext() as context:
        context.prec = 50
        for field in ('RAW','TL'):
            delta = 'delta_C_%s_p4_minus_p40' % field
            absolute = max(combined, key=lambda r: abs(r[delta]))
            relative = max(combined, key=lambda r: abs(r[delta]/r['C_%s_p40' % field]))
            projection.append(dict(definition=field,
                maximum_absolute_p4_minus_p40=str(abs(absolute[delta])),
                maximum_absolute_position=position(absolute),
                maximum_relative_p4_minus_p40_percent=str(100*abs(relative[delta]/relative['C_%s_p40' % field])),
                maximum_relative_position=position(relative)))
    midpoint = []
    for r in combined:
        if tuple(r[k] for k in ('R','long','tr1','tr2')) == (4,2,0,0):
            midpoint.append({k:str(v) if isinstance(v,D) else v for k,v in r.items()})
    return dict(success=True,configurations=TARGET,excluded_configuration=640,configurations_count=70,
         cases=210,taus=[2,4,6],piters=list(PITERS),raw_records=139440,window_records=420,
         ensemble_geometry_rows=len(records),combined_geometry_rows=len(combined),
         raw_regression_exact=True,raw_regression_max_absolute_difference='0',
         maximum_direct_clover_identity_residual=report['maximum_direct_clover_identity_residual'],
         convention=dict(field_strength='F=(Q-Q^dagger)/8 on the original gauge field',
            plane_sum='Six unordered planes mu<nu; no additional factor 1/2',
            S_RAW='sum Tr(F^dagger F)',T_TRACE='sum |Tr(F)|^2/3',S_TL='S_RAW-T_TRACE',
            C_RAW='sum(LS_RE)/sum(L_RE)-mean(S_RAW)',
            C_TL='(sum(LS_RE)-sum(LT_RE))/sum(L_RE)-mean(S_RAW-T_TRACE)',
            C_TRACE='C_RAW-C_TL=sum(LT_RE)/sum(L_RE)-mean(T_TRACE)',
            averaging='Equal configuration weights; form ensemble ratios after averaging primaries'),
         statistics=dict(method='Paired delete-one-configuration jackknife; omit the same configuration across definitions, piters and geometries',
            standard_error='sqrt((N-1)/N*sum((leave_one_out_estimate-mean_leave_one_out_estimate)^2))',
            point_estimate='Full-sample ratio of means; no jackknife bias correction',
            limitations='Preliminary errors without autocorrelation correction. Windows, axes and probe positions are correlated. Configuration label spacing is not a verified Monte Carlo time interval. Geometry maxima are descriptive.'),
         trace_comparison=comparisons,projection_comparison=projection,R4_midpoint=midpoint,
         scope='One matched 70-configuration comparison of raw and traceless clover at piter4/40. Existing raw piter10/20 curves remain in the projection study; traceless piter10/20 primaries were not measured. Bare lattice observables; no renormalized physical-density claim.')


def numerical_exports(output):
    sources, report = verify_inputs()
    rows = read_csv(HERE/'inputs/profiles_raw_and_trace_70configs.csv')
    groups = group_primaries(rows)
    verify_raw_regression(rows, sources)
    windows = read_csv(HERE/'inputs/su3_windows_70configs.csv')
    validate_windows(windows)
    records = []
    for key, values in sorted(groups.items()):
        row = dict(zip(KEYS,key))
        row['configurations'] = 70
        row.update(paired_statistics(values,TARGET))
        records.append(row)
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    write_csv(output/'profiles_ensemble_70configs.csv',records)
    summary = summarize(records,report)
    write_csv(output/'R4_midpoint_70configs.csv',summary['R4_midpoint'])
    write_csv(output/'projection_comparison_70configs.csv',summary['projection_comparison'])
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return records, summary


def verify():
    with tempfile.TemporaryDirectory() as directory:
        numerical_exports(Path(directory))
        for name in ('profiles_ensemble_70configs.csv','R4_midpoint_70configs.csv',
                     'projection_comparison_70configs.csv','SUMMARY.json'):
            if (Path(directory)/name).read_bytes() != (RESULTS/name).read_bytes():
                raise ValueError('Numerical export differs: '+name)
    print('PASS full70 clover: 70 configurations; 210 cases; 139440 primary rows; exact raw regression; all ensemble ratios and paired jackknife errors reproduced.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--plots',action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify()
    else:
        records, _ = numerical_exports(RESULTS)
        if args.plots:
            import render_plots
            render_plots.render(records,RESULTS)
        print('PASS full70 clover numerical exports')


if __name__ == '__main__':
    main()
