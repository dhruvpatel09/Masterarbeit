#!/usr/bin/env python3
"""Validate raw and traceless clover primaries in a matched n100 pilot."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import time
import zipfile
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EXE = HERE/'build/test_profiles_trace_qx'
TRACE_FIELDS = ('T_TRACE', 'LT_RE', 'LT_IM')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def verify_sources():
    metadata = json.loads((HERE/'INPUTS.json').read_text())
    for name, expected in metadata['files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Source changed: '+name)
    legacy = module('qx_trace_legacy', ROOT/'laplace_static/validation/qx_projection_spectrum/campaign.py')
    legacy.verify_sources()
    return legacy


def build_manifest():
    path = HERE/'build/manifest.sha256'
    for line in path.read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        target = Path(name)
        if not target.is_absolute():
            target = ROOT/target
        if digest(target) != expected:
            raise ValueError('Build changed: '+name)
    if not EXE.is_file():
        raise ValueError('Build the trace pilot executable first')
    return digest(path)


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}', name):
        raise ValueError('Invalid pilot name')
    return HERE/'campaigns'/name


def case_plan():
    return [dict(configuration=100, tau=4, piter=p) for p in (4, 40)]


def prepare(args):
    verify_sources()
    manifest = build_manifest()
    full = module('qx_trace_full', ROOT/'laplace_static/validation/qx_projection_full_ensemble/campaign.py')
    _, known = full.inventory()
    if 100 not in known:
        raise ValueError('The selected-central-time n100 saved-mode audit is required')
    run = known[100]
    marker = Path(run['run'])/'COMPLETED.json'
    if digest(marker) != run['completion_marker_sha256']:
        raise ValueError('Audited n100 completion marker changed')
    completion = json.loads(marker.read_text())
    if (completion['configuration'], completion['timeslices'], completion['eigenvector_files']) != (100, 192, 192):
        raise ValueError('Incomplete n100 generation inventory')
    legacy = verify_sources()
    inputs = {}
    paths = [legacy.GAUGES/'qX_p2n100', marker]
    paths += [Path(run['run'])/'eigenvectors'/('qX_p2n100_evec_t%d.h5' % t) for t in (96, 100)]
    for path in paths:
        stat = path.stat()
        if stat.st_size <= 0:
            raise ValueError('Empty measurement input: '+str(path))
        inputs[str(path)] = dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    data = dict(name=args.name, cases=case_plan(), run=run, inputs=inputs,
                executable_sha256=digest(EXE), build_manifest_sha256=manifest,
                source_manifest_sha256=digest(HERE/'INPUTS.json'),
                settings=dict(Nv=100, midpoint=98, alpha=[1, 1, .5], hyp_steps=1,
                              grid=[1, 4, 4, 4], ranks=64, nodes=8),
                scope='Two-launch n100 tau4 validation pilot; raw and trace primaries from the same field and loops')
    target = folder(args.name)
    target.mkdir(parents=True, exist_ok=False)
    write_json(target/'PREPARED.json', data)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared n100 tau4: piter4 and piter40; two MPI launches; both clover definitions')


def prepared(name):
    verify_sources()
    target = folder(name)
    if digest(target/'PREPARED.json') != (target/'prepared.sha256').read_text().strip():
        raise ValueError('Preparation changed')
    data = json.loads((target/'PREPARED.json').read_text())
    if (data['cases'] != case_plan() or data['source_manifest_sha256'] != digest(HERE/'INPUTS.json') or
            data['executable_sha256'] != digest(EXE) or data['build_manifest_sha256'] != build_manifest()):
        raise ValueError('Prepared implementation changed')
    for name, expected in data['inputs'].items():
        stat = Path(name).stat()
        if (stat.st_size, stat.st_mtime_ns) != (expected['bytes'], expected['mtime_ns']):
            raise ValueError('Prepared input changed: '+name)
    return target, data


def trace_records(text, raw, avg):
    checks = avg.rows(text, 'CLOVER_TRACE_CHECK')
    if (len(checks) != 1 or checks[0]['status'] != 'PASS' or
            int(checks[0]['count']) != 64**3):
        raise ValueError('Missing direct traceless-clover identity check')
    residual = float(checks[0]['identity_max'])
    if not math.isfinite(residual) or not 0 <= residual <= 1e-11:
        raise ValueError('Clover trace identity failed')
    probes = avg.rows(text, 'CLOVER_TRACE_PROBE')
    if len(probes) != 3 or {r['kind'] for r in probes} != {'local', 'cross_boundary', 'spatial_average'}:
        raise ValueError('Missing clover trace probes')
    raw_probes = {r['kind']: r for r in avg.rows(text, 'PROBE')}
    for row in probes:
        values = [float(row[k]) for k in ('T_TRACE', 'S_TL')]
        if int(row['t']) != 98 or not all(math.isfinite(v) and v >= -1e-13 for v in values):
            raise ValueError('Invalid trace probe')
        avg.near(sum(values), float(raw_probes[row['kind']]['S_RAW']), 'local trace identity')
    vacuum = float(next(r['T_TRACE'] for r in probes if r['kind'] == 'spatial_average'))
    records = {}
    for prefix in ('TRACE_AXIS', 'TRACE_COMBINED'):
        for row in avg.rows(text, prefix):
            axis = int(row['axis']) if prefix == 'TRACE_AXIS' else -1
            key = (axis, *(int(row[k]) for k in ('R', 'long', 'tr1', 'tr2')))
            if key in records or key not in raw:
                raise ValueError('Duplicate or unexpected trace profile')
            if int(row['Nsrc']) != int(raw[key]['Nsrc']):
                raise ValueError('Trace/raw source-count mismatch')
            if not all(math.isfinite(float(row[k])) for k in TRACE_FIELDS):
                raise ValueError('Nonfinite trace primary')
            T, S = float(row['T_TRACE']), float(raw[key]['S_RAW'])
            if not -1e-13 <= T <= S+1e-12:
                raise ValueError('Invalid nonnegative trace component')
            avg.near(T, vacuum, 'translated trace-vacuum average')
            records[key] = row
    if records.keys() != raw.keys():
        raise ValueError('Incomplete matched trace profile')
    for key, row in records.items():
        axis, r, l, x, y = key
        if axis == -1:
            for field in TRACE_FIELDS:
                values = [float(records[(a, r, l, x, y)][field]) for a in range(3)]
                avg.near(float(row[field]), math.fsum(values)/3, 'combined trace '+field,
                         math.fsum(map(abs, values))/3)
    return records


def connected_from_primaries(row):
    """row contains configuration means, or already formed ensemble means."""
    with localcontext() as context:
        context.prec = 50
        L, LS, S, T, LT = (Decimal(row[k]) for k in ('L_RE', 'LS_RE', 'S_RAW', 'T_TRACE', 'LT_RE'))
        if not all(v.is_finite() for v in (L, LS, S, T, LT)) or L == 0:
            raise ValueError('Invalid primaries or zero loop mean')
        raw, correction = LS/L-S, LT/L-T
        return dict(C_RAW=raw, C_TRACE_COMPONENT=correction,
                    C_TL=(LS-LT)/L-(S-T))


def validate(target):
    legacy = verify_sources()
    check = legacy.checker()
    scans, traces, rows = {}, {}, []
    for entry in case_plan():
        p = entry['piter']
        path = target/'outputs'/('p%d.out' % p)
        scans[p] = check.read_measurement(path, 100, 4, p)[1]
        traces[p] = trace_records(path.read_text(), scans[p], check.avg)
        for row in check.raw_rows(100, 4, p, scans[p]):
            key = tuple(int(row[k]) for k in ('axis', 'R', 'long', 'tr1', 'tr2'))
            row.update({k: traces[p][key][k] for k in TRACE_FIELDS})
            connected_from_primaries(row)
            rows.append(row)
    check.regression_check(scans)
    for key in scans[4]:
        if Decimal(scans[4][key]['S_RAW']) != Decimal(scans[40][key]['S_RAW']):
            raise ValueError('Original raw clover changed across piters')
        if Decimal(traces[4][key]['T_TRACE']) != Decimal(traces[40][key]['T_TRACE']):
            raise ValueError('Original trace component changed across piters')
    check.write_csv(target/'profiles_raw_and_trace.csv', check.RAW_COLUMNS+list(TRACE_FIELDS), rows)
    return dict(success=True, configuration=100, tau=4, piters=[4, 40], records=len(rows),
                raw_regression_job=5915129, geometry_records_per_launch=332,
                scope='Single-configuration validation only; no ensemble or autocorrelation claim',
                convention=json.loads((HERE/'INPUTS.json').read_text())['convention'])


def run(args):
    target, data = prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS', '0')) != 64 or int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 8:
        raise ValueError('Run in the submitted 8-node, 64-rank allocation')
    out = target/'outputs'
    out.mkdir(exist_ok=False)
    legacy = verify_sources()
    launches = []
    try:
        for entry in data['cases']:
            p = entry['piter']
            env = os.environ.copy()
            env.update(QX_CFG='100', QX_ROOT=str(legacy.GAUGES),
                       EVROOT=str(Path(data['run']['run'])/'eigenvectors'), T0='96', TAU_LENGTH='4',
                       HYP_MODE='2', HYP_PITER=str(p), PRINT_ALL='1', GX='0', GY='0', GZ='0')
            command = ['mpiexec', '-n', '64', str(EXE)]
            start = time.monotonic()
            with (out/('p%d.out' % p)).open('x') as stdout, (out/('p%d.err' % p)).open('x') as stderr:
                code = subprocess.run(command, env=env, stdout=stdout, stderr=stderr, check=False).returncode
            launches.append(dict(entry, exit_code=code, seconds=time.monotonic()-start,
                                 job_id=os.environ['SLURM_JOB_ID']))
            write_json(out/('launch_p%d.json' % p), launches[-1])
            if code:
                raise ValueError('MPI trace measurement failed at piter%d' % p)
        write_json(target/'LAUNCHES.json', launches)
        write_json(target/'SUMMARY.json', validate(target))
    except Exception as exc:
        write_json(target/'FAILED.json', dict(success=False, error=str(exc), launches=launches))
        raise
    print('PASS raw/traceless clover pilot')


def submit(args):
    target, _ = prepared(args.name)
    if (target/'SUBMISSION.json').exists() or (target/'outputs').exists():
        raise ValueError('Pilot has a submission/output record; inspect before submitting again')
    command = ['sbatch', '--parsable', str(HERE/'submit_micro.sh'), args.name]
    write_json(target/'SUBMISSION.json', dict(state='submitting', command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=False)
    (target/'SUBMISSION.json').write_text(json.dumps(dict(command=command, returncode=result.returncode,
                                                       stdout=result.stdout, stderr=result.stderr), indent=2)+'\n')
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?', result.stdout):
        raise ValueError('Submission did not return a clear job ID; inspect SUBMISSION.json')
    print('SUBMITTED', result.stdout.strip(), 'raw/traceless clover pilot')


def collect(args):
    target, _ = prepared(args.name)
    if (target/'FAILED.json').exists():
        raise ValueError('Pilot has a FAILED.json; inspect before collection')
    launches = json.loads((target/'LAUNCHES.json').read_text())
    if ([{k: row[k] for k in ('configuration', 'tau', 'piter')} for row in launches] != case_plan() or
            any(row['exit_code'] for row in launches)):
        raise ValueError('Incomplete pilot launch inventory')
    summary = validate(target)
    if json.loads((target/'SUMMARY.json').read_text()) != summary:
        raise ValueError('Saved pilot validation differs')
    output = ROOT/('qx_clover_trace_%s_review.zip' % args.name)
    paths = sorted(path for path in target.rglob('*') if path.is_file())
    paths += [HERE/'INPUTS.json', HERE/'build/manifest.sha256']
    with zipfile.ZipFile(str(output), 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            archive.write(str(path), path.relative_to(ROOT).as_posix())
    print('PASS raw regression, direct clover identity, trace translations and axis reductions')
    print('REVIEW_ARCHIVE', output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    sub.add_parser('verify')
    for action in ('prepare', 'submit', 'run', 'collect'):
        command = sub.add_parser(action)
        command.add_argument('--name', default='pilot')
    args = parser.parse_args()
    if not args.action:
        parser.error('Choose verify, prepare, submit, run or collect')
    if args.action == 'verify':
        verify_sources()
        print('PASS trace pilot source pins')
    else:
        globals()[args.action](args)


if __name__ == '__main__':
    main()
