#!/usr/bin/env python3
"""Reproduce paired qX temporal-separation contrasts on the matched 70-config cohort."""
import argparse
import csv
import hashlib
import importlib.util
import json
import tempfile
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / 'results'
UPSTREAM = HERE.parent / 'qx_clover_trace_full_ensemble'
D = Decimal
TAUS = (2, 4, 6)
PITERS = (4, 40)
DEFINITIONS = ('RAW', 'TL')
GEOMETRY = ('axis', 'R', 'long', 'tr1', 'tr2')
CONTRASTS = ((4, 2), (6, 4), (6, 2))
EXPORTS = ('temporal_comparison_70configs.csv', 'temporal_covariance_70configs.csv',
           'projection_control_70configs.csv', 'R4_midpoint_70configs.csv', 'SUMMARY.json')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
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


def verified_upstream():
    sources = json.loads((HERE / 'inputs/SOURCES.json').read_text())
    for name, digest in sources['repository_references'].items():
        if sha(ROOT / name) != digest:
            raise ValueError('Upstream input changed: ' + name)
    spec = importlib.util.spec_from_file_location('qx_clover_temporal_source', UPSTREAM / 'reproduce.py')
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    report = json.loads((UPSTREAM / 'inputs/VALIDATION.json').read_text())
    if (report['success'], report['configurations'], report['cases'], report['raw_records'],
            report['raw_regression_exact']) != (True, source.TARGET, 210, 139440, True):
        raise ValueError('Upstream cohort is not the validated matched full ensemble')
    return source, sources


def jackknife_covariance(samples):
    """Covariance of the estimator from synchronized leave-one-out vectors."""
    n = len(samples)
    if n < 2 or not samples[0] or any(len(v) != len(samples[0]) for v in samples):
        raise ValueError('Need at least two equal-length jackknife vectors')
    if any(not x.is_finite() for row in samples for x in row):
        raise ValueError('Nonfinite jackknife estimate')
    with localcontext() as context:
        context.prec = 50
        m = len(samples[0])
        means = [sum((v[i] for v in samples), D(0)) / D(n) for i in range(m)]
        cov = []
        for i in range(m):
            row = []
            for j in range(m):
                # Keep the upstream variance expression identical on the diagonal.
                products = ((v[i] - means[i]) ** 2 if i == j else
                            (v[i] - means[i]) * (v[j] - means[j]) for v in samples)
                row.append(D(n - 1) / D(n) * sum(products, D(0)))
            cov.append(row)
        return cov


def paired_estimates(values_by_tau, configs, source):
    """Recompute each ratio after dropping the same configuration at every tau/piter."""
    n = len(configs)
    if n < 2 or len(set(configs)) != n or set(values_by_tau) != set(TAUS):
        raise ValueError('Invalid matched configuration or temporal inventory')
    required = {(c, p) for c in configs for p in PITERS}
    if any(set(values_by_tau[t]) != required for t in TAUS):
        raise ValueError('Temporal separations must use an identical matched cohort')
    with localcontext() as context:
        context.prec = 50
        totals = {t: {p: tuple(sum((values_by_tau[t][c, p][j] for c in configs), D(0))
                             for j in range(5)) for p in PITERS} for t in TAUS}
        point = {t: source.observables(totals[t], n) for t in TAUS}
        leave_one_out = []
        for c in configs:
            sample = {}
            for t in TAUS:
                sums = {p: tuple(totals[t][p][j] - values_by_tau[t][c, p][j]
                                 for j in range(5)) for p in PITERS}
                sample[t] = source.observables(sums, n - 1)
            leave_one_out.append(sample)
        return point, leave_one_out


def channel_statistics(point, leave_one_out, definition, piter):
    field = 'C_%s_p%d' % (definition, piter)
    with localcontext() as context:
        context.prec = 50
        samples = [[s[t][field] for t in TAUS] for s in leave_one_out]
        covariance = jackknife_covariance(samples)
        result = {}
        for i, t in enumerate(TAUS):
            result['C_tau%d' % t] = point[t][field]
            result['C_tau%d_jackknife_se' % t] = covariance[i][i].sqrt()
        for high, low in CONTRASTS:
            name = 'delta_tau%d_minus_tau%d' % (high, low)
            differences = [[s[high][field] - s[low][field]] for s in leave_one_out]
            variance = jackknife_covariance(differences)[0][0]
            i, j = TAUS.index(high), TAUS.index(low)
            matrix_variance = covariance[i][i] + covariance[j][j] - 2 * covariance[i][j]
            if abs(variance - matrix_variance) > max(abs(variance), D('1e-90')) * D('1e-44'):
                raise ValueError('Direct paired variance disagrees with covariance contraction')
            result[name] = point[high][field] - point[low][field]
            result[name + '_jackknife_se'] = variance.sqrt()
            result[name + '_se_ignoring_tau_covariance'] = (covariance[i][i] + covariance[j][j]).sqrt()
        return result, covariance


def projection_controls(point, leave_one_out, definition):
    """Compare temporal contrasts between p4 and p40 within each jackknife sample."""
    with localcontext() as context:
        context.prec = 50
        f4, f40 = 'C_%s_p4' % definition, 'C_%s_p40' % definition
        rows = []
        for high, low in CONTRASTS:
            delta = lambda r: (r[high][f4] - r[low][f4]) - (r[high][f40] - r[low][f40])
            error = jackknife_covariance([[delta(s)] for s in leave_one_out])[0][0].sqrt()
            rows.append(dict(tau_high=high, tau_low=low, contrast_p4_minus_p40=delta(point),
                             contrast_p4_minus_p40_jackknife_se=error))
        return rows


def numerical_exports(output):
    source, sources = verified_upstream()
    rows = read_csv(UPSTREAM / 'inputs/profiles_raw_and_trace_70configs.csv')
    groups = source.group_primaries(rows)
    archived = {tuple(int(r[k]) for k in source.KEYS): r
                for r in read_csv(UPSTREAM / 'results/profiles_ensemble_70configs.csv')}
    spatial = sorted({key[1:] for key in groups})
    records, covariance_rows, controls = [], [], []
    with localcontext() as context:
        context.prec = 50
        for geometry in spatial:
            point, loo = paired_estimates({t: groups[(t,) + geometry] for t in TAUS},
                                          source.TARGET, source)
            for definition in DEFINITIONS:
                for p in PITERS:
                    stats, covariance = channel_statistics(point, loo, definition, p)
                    row = dict(zip(GEOMETRY, geometry))
                    row.update(definition=definition, piter=p, configurations=len(source.TARGET))
                    row.update(stats)
                    records.append(row)
                    field = 'C_%s_p%d' % (definition, p)
                    for t in TAUS:
                        reference = archived[(t,) + geometry]
                        if (stats['C_tau%d' % t] != D(reference[field]) or
                                stats['C_tau%d_jackknife_se' % t] != D(reference[field + '_jackknife_se'])):
                            raise ValueError('Temporal marginal disagrees with archived profile')
                    for i, t in enumerate(TAUS):
                        for j, u in enumerate(TAUS):
                            product = covariance[i][i] * covariance[j][j]
                            rho = covariance[i][j] / product.sqrt() if product > 0 else ''
                            covrow = dict(zip(GEOMETRY, geometry))
                            covrow.update(definition=definition, piter=p, tau_i=t, tau_j=u,
                                          covariance=covariance[i][j], correlation=rho)
                            covariance_rows.append(covrow)
                for control in projection_controls(point, loo, definition):
                    row = dict(zip(GEOMETRY, geometry))
                    row.update(definition=definition, configurations=len(source.TARGET))
                    row.update(control)
                    controls.append(row)
        midpoint = [r for r in records if tuple(r[k] for k in GEOMETRY) == (-1, 4, 2, 0, 0)]
        combined = [r for r in records if r['axis'] == -1 and r['piter'] == 40]
        summary = dict(success=True, base_commit=sources['base_commit'],
            configurations=source.TARGET, configurations_count=70, excluded_configuration=640,
            taus=list(TAUS), piters=list(PITERS), definitions=list(DEFINITIONS),
            input_primary_rows=len(rows), original_profile_geometries=len(groups),
            spatial_geometries=len(spatial), combined_spatial_geometries=83,
            comparison_rows=len(records), covariance_rows=len(covariance_rows),
            projection_control_rows=len(controls), archived_marginals_reproduced_exactly=True,
            statistics=dict(point_estimate='Ratio of configuration means, without bias correction',
                method='Paired delete-one-configuration jackknife. All taus, definitions and piters omit the same configuration.',
                covariance='(N-1)/N * sum((theta_leave_one_out-mean_leave_one_out)_i * (theta_leave_one_out-mean_leave_one_out)_j)',
                contrast='Recompute each ratio at each tau within each jackknife sample, then take the temporal difference.',
                precision='50 decimal digits; saved primaries set the input precision',
                uncertainty_status='Preliminary; no autocorrelation correction',
                sampling='70 configuration samples. Windows, spatial probes and axes are correlated; their count is not an independent sample size.'),
            monte_carlo_metadata=sources['monte_carlo_metadata'],
            R4_midpoint=midpoint,
            combined_axis_p40_maximum_absolute_temporal_contrasts=[
                dict(definition=definition, tau_high=high, tau_low=low,
                     maximum_absolute_difference=max(abs(r['delta_tau%d_minus_tau%d' % (high, low)])
                        for r in combined if r['definition'] == definition))
                for definition in DEFINITIONS for high, low in CONTRASTS],
            limitations=[
                'No temporal plateau fit or excited-state model is performed. Agreement at one point or within current errors does not establish ground-state saturation.',
                'The reported errors do not account for Markov-chain autocorrelation. Chronological order, saving interval, streams and thermalization remain unverified.',
                'Configuration n640 is absent. Configuration identifiers are preserved without treating their numerical spacing as Monte Carlo time.',
                'The two adjacent contrasts share tau4 and are correlated. The nonadjacent contrast is their sum; the three are not independent.',
                'Raw clover remains the main convention; traceless clover is the side comparison. Both are bare lattice observables.',
                'The fixed midpoint t/a=98, tau/a=2,4,6, R/a=1,2,4 and Nv=100 define this study. No new projection prescriptions were measured.'])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / EXPORTS[0], records)
    write_csv(output / EXPORTS[1], covariance_rows)
    write_csv(output / EXPORTS[2], controls)
    write_csv(output / EXPORTS[3], midpoint)
    (output / EXPORTS[4]).write_text(json.dumps(summary, indent=2, allow_nan=False,
                                               default=lambda v: str(v) if isinstance(v, D) else v) + '\n')
    return records, covariance_rows, controls, summary


def verify():
    with tempfile.TemporaryDirectory() as directory:
        numerical_exports(directory)
        for name in EXPORTS:
            if (Path(directory) / name).read_bytes() != (RESULTS / name).read_bytes():
                raise ValueError('Numerical export differs: ' + name)
    print('PASS paired temporal analysis: 70 configurations; 210 cases; all 996 archived profile marginals reproduced exactly; 1328 comparison rows; 11952 covariance entries; 1992 projection controls.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    if args.verify:
        if args.plots:
            parser.error('--verify compares numerical exports without rewriting plots')
        verify()
    else:
        records, covariance, controls, _ = numerical_exports(RESULTS)
        if args.plots:
            import render_plots
            render_plots.render(records, covariance, controls, RESULTS)
        print('PASS paired temporal numerical exports')


if __name__ == '__main__':
    main()
