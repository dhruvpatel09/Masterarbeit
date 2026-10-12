#!/usr/bin/env python3
"""Recompute the n100 finite-projection comparison from its uploaded CSVs.

Only Python's standard library is required. Values derived with Decimal are
diagnostics computed from the printed input precision, not new measurements.
"""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from pathlib import Path
from zipfile import ZipFile

getcontext().prec = 50
KEYS = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
FIELDS = ('L_RE', 'L_IM', 'S_RAW', 'LS_RE', 'LS_IM')


def read_csv(bundle, fragment):
    names = [n for n in bundle.namelist() if fragment in n and n.endswith('.csv')]
    if len(names) != 1:
        raise ValueError(f'Expected exactly one {fragment} CSV')
    return list(csv.DictReader(io.StringIO(bundle.read(names[0]).decode())))


def offsets(r):
    return {(0, 1, 2)} | {(l, t, 0) for l in range(-1, r + 2) for t in range(-2, 3)}


def key(row):
    return tuple(int(row[f]) for f in KEYS)


def complex_ratio(p):
    lr, li = p['L_RE'], p['L_IM']
    sr, si = p['LS_RE'], p['LS_IM']
    denom = lr * lr + li * li
    if denom == 0:
        raise ValueError('Zero complex denominator')
    return ((sr * lr + si * li) / denom - p['S_RAW'],
            (si * lr - sr * li) / denom)


def maximum(records, field):
    row = max(records, key=lambda r: abs(r[field]))
    return {'absolute_value': str(abs(row[field])),
            'signed_value': str(row[field]), 'key': list(row['key'])}


def summary(records):
    columns = ['C_difference', 'C_fractional_change',
               'complex_C_RE_difference', 'complex_C_IM_difference']
    columns += [f'{f}_{s}' for f in FIELDS if f != 'S_RAW'
                for s in ('difference', 'fractional_change')]
    result = {'pairs': len(records), 'maxima': {f: maximum(records, f) for f in columns}}
    result['C_p40_min'] = str(min(r['C_p40'] for r in records))
    result['C_p40_max'] = str(max(r['C_p40'] for r in records))
    result['C_difference_rms'] = str((sum(r['C_difference']**2 for r in records)
                                        / len(records)).sqrt())
    return result


def recompute(archive):
    with ZipFile(archive) as bundle:
        if bundle.testzip() is not None:
            raise ValueError('Input archive integrity failure')
        raw = read_csv(bundle, '/profile_raw_')
        exported = read_csv(bundle, '/profile_piter4_minus40_')
        errors = [bundle.read(n).decode() for n in bundle.namelist() if n.endswith('.err')]
        for text in errors:
            if 'MPI_Abort' in text or 'Abort(' in text or 'FAIL:' in text:
                raise ValueError('Application error in supplied logs')
        member_hashes = {n: hashlib.sha256(bundle.read(n)).hexdigest()
                         for n in bundle.namelist()}

    expected = {(tau, axis, r, *offset) for tau in (2, 4, 6)
                for axis in (-1, 0, 1, 2) for r in (1, 2, 4)
                for offset in offsets(r)}
    pairs = defaultdict(dict)
    for row in raw:
        k = key(row)
        piter = int(row['piter'])
        if row['configuration'] != 'qX_p2n100' or int(row['Nv']) != 100:
            raise ValueError('Configuration or Nv mismatch')
        if (int(row['tmid']), int(row['t0']), int(row['t1'])) != (98, 98-k[0]//2, 98+k[0]//2):
            raise ValueError('Temporal interval mismatch')
        if Decimal(row['long_from_midpoint']) != Decimal(k[3])-Decimal(k[2])/2:
            raise ValueError('Probe midpoint coordinate mismatch')
        nsrc = (3 if k[1] == -1 else 1) * 64**3
        if int(row['Nsrc']) != nsrc:
            raise ValueError('Source count mismatch')
        if (piter, row['mode']) not in ((4, 'hyp4'), (40, 'hyp40_reference')):
            raise ValueError('Mode or projection-iteration mismatch')
        if piter in pairs[k]:
            raise ValueError('Duplicate input record')
        values = {f: Decimal(row[f]) for f in FIELDS}
        if any(not v.is_finite() for v in values.values()):
            raise ValueError('Nonfinite input')
        pairs[k][piter] = values
    if set(pairs) != expected or len(raw) != 1992:
        raise ValueError('Unexpected profile key set or raw record count')
    if any(set(p) != {4, 40} for p in pairs.values()):
        raise ValueError('Missing projection partner')

    exp = {}
    for row in exported:
        k = key(row)
        if k in exp:
            raise ValueError('Duplicate exported difference')
        exp[k] = row
    if set(exp) != expected or len(exported) != 996:
        raise ValueError('Difference export key/count mismatch')

    max_combination_residual = {f: Decimal(0) for f in FIELDS}
    for k in expected:
        if k[1] != -1:
            continue
        for piter in (4, 40):
            for f in FIELDS:
                mean = sum(pairs[(k[0], ax, *k[2:])][piter][f] for ax in (0, 1, 2))/3
                residual = abs(pairs[k][piter][f]-mean)
                max_combination_residual[f] = max(max_combination_residual[f], residual)
                # The combined record was produced from a separate floating-
                # point reduction; allow rounding at the printed precision.
                scale = max(abs(mean), Decimal('1e-20'))
                if residual > Decimal('2e-14')*scale:
                    raise ValueError('Combined-axis reduction inconsistency')

    records = []
    for k in sorted(pairs):
        a, b = pairs[k][4], pairs[k][40]
        if a['S_RAW'] != b['S_RAW']:
            raise ValueError('Original clover vacuum average changed')
        out = {'key': k}
        for f in FIELDS:
            # Independently reproduce the original exporter's binary-float
            # difference and signed/absolute-reference fractional convention.
            af, bf = float(a[f]), float(b[f])
            for suffix, value in (('p4', af), ('p40', bf), ('difference', af-bf)):
                if float(exp[k][f+'_'+suffix]) != value:
                    raise ValueError(f'Export mismatch: {k} {f} {suffix}')
            fraction = (af-bf)/abs(bf) if bf else None
            printed = exp[k][f+'_fractional_change']
            if ((fraction is None and printed != '') or
                    (fraction is not None and float(printed) != fraction)):
                raise ValueError('Export fractional-change mismatch')
            out[f+'_difference'] = a[f]-b[f]
            out[f+'_fractional_change'] = (a[f]-b[f])/abs(b[f]) if b[f] else Decimal(0)
        if a['L_RE'] == 0 or b['L_RE'] == 0:
            raise ValueError('Zero denominator for real-parts ratio')
        out['C_p4'] = a['LS_RE']/a['L_RE']-a['S_RAW']
        out['C_p40'] = b['LS_RE']/b['L_RE']-b['S_RAW']
        if out['C_p40'] == 0:
            raise ValueError('Undefined fractional connected diagnostic')
        out['C_difference'] = out['C_p4']-out['C_p40']
        out['C_fractional_change'] = out['C_difference']/abs(out['C_p40'])
        car, cai = complex_ratio(a)
        cbr, cbi = complex_ratio(b)
        out.update(complex_C_RE_p4=car, complex_C_RE_p40=cbr,
                   complex_C_IM_p4=cai, complex_C_IM_p40=cbi,
                   complex_C_RE_difference=car-cbr, complex_C_IM_difference=cai-cbi)
        records.append(out)
    summary_data = {
        'input_archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
        'input_members_sha256': member_hashes,
        'precision': 'Decimal arithmetic at precision 50 on printed inputs',
        'scope': 'Single configuration qX_p2n100; raw clover convention; no ensemble uncertainty',
        'checks': {'raw_rows': len(raw), 'pairs': len(records),
                   'pairs_per_tau': dict(Counter(k[0] for k in pairs)),
                   'export_matches_float_recomputation': True,
                   'original_S_RAW_identical_in_all_pairs': True,
                   'combined_axis_max_absolute_residual':
                       {f: str(v) for f, v in max_combination_residual.items()}},
        'all_records': summary(records),
        'combined_axes': summary([r for r in records if r['key'][1] == -1]),
        'by_tau': {str(t): summary([r for r in records if r['key'][0] == t]) for t in (2,4,6)},
        'combined_by_tau': {str(t): summary([r for r in records
                                            if r['key'][0] == t and r['key'][1] == -1])
                            for t in (2,4,6)},
        'by_tau_axis_R': {f'{t}/{ax}/{r}': summary([v for v in records
                                                  if v['key'][:3] == (t, ax, r)])
                         for t in (2,4,6) for ax in (-1,0,1,2) for r in (1,2,4)},
    }
    return records, summary_data


def write_outputs(records, summaries, destination):
    destination.mkdir(parents=True, exist_ok=False)
    columns = [f for f in records[0] if f != 'key']
    with (destination/'vacuum_subtracted_diagnostics.csv').open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[*KEYS, *columns])
        writer.writeheader()
        for row in records:
            writer.writerow(dict(zip(KEYS, row['key'])) | {c: str(row[c]) for c in columns})
    (destination/'summary.json').write_text(json.dumps(summaries, indent=2)+'\n')
    print('PASS input inventory, source counts, pairing, finite values, export recomputation, '
          'unchanged clover averages and combined-axis reductions')
    print(json.dumps({'all_records': summaries['all_records'],
                      'combined_axes': summaries['combined_axes']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    records, summaries = recompute(args.archive)
    write_outputs(records, summaries, args.destination)
