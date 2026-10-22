#!/usr/bin/env python3
"""Measure a matched 70-configuration raw/traceless pair and a central piter pilot."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import time
import zipfile
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = [n for n in range(100, 801, 10) if n != 640]
TAUS = [2, 4, 6]
SPECTRUM_CONFIGS = [100, 170, 670, 740]
GEOMETRY = ('axis', 'R', 'long', 'tr1', 'tr2')
REFERENCE = ROOT/'laplace_static/analysis/qx_projection_full_ensemble/results/profiles_raw_70configs.csv'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


pilot = module('qx_trace_ensemble_pilot', ROOT/'laplace_static/validation/qx_clover_trace/campaign.py')
full = module('qx_trace_ensemble_full', ROOT/'laplace_static/validation/qx_projection_full_ensemble/campaign.py')
analysis = module('qx_trace_ensemble_reference', ROOT/'laplace_static/analysis/qx_clover_trace_pilot/reproduce.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def verify_sources():
    pilot.verify_sources()
    full.verify_sources()
    for name, expected in json.loads((HERE/'INPUTS.json').read_text())['files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Ensemble workflow source changed: '+name)


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}', name):
        raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def case_plan(stage):
    if stage not in ('full70', 'spectrum'):
        raise ValueError('Choose full70 or spectrum')
    configs = TARGET if stage == 'full70' else SPECTRUM_CONFIGS
    piters = [4, 40] if stage == 'full70' else [6, 8, 80]
    pairs = [(c, t) for c in configs for t in TAUS if stage != 'full70' or (c, t) != (100, 4)]
    return [dict(index=i, configuration=c, tau=t, piters=piters[:]) for i, (c, t) in enumerate(pairs)]


def slot_plan(stage):
    cases = case_plan(stage)
    size = 6 if stage == 'full70' else 1
    return [list(range(i, min(i+size, len(cases)))) for i in range(0, len(cases), size)]


def prepare(args):
    verify_sources()
    analysis.verify()
    manifest = pilot.build_manifest()
    _, known = full.inventory()
    configs = TARGET if args.stage == 'full70' else SPECTRUM_CONFIGS
    missing = [c for c in configs if c not in known]
    if missing:
        raise ValueError('Missing selected-central-time audit: '+','.join(map(str, missing)))
    inputs, runs = {}, []
    for cfg in configs:
        run = known[cfg]
        marker = Path(run['run'])/'COMPLETED.json'
        if digest(marker) != run['completion_marker_sha256']:
            raise ValueError('Audited completion marker changed for n%d' % cfg)
        done = json.loads(marker.read_text())
        if (done['configuration'], done['timeslices'], done['eigenvector_files']) != (cfg, 192, 192):
            raise ValueError('Incomplete saved-mode inventory for n%d' % cfg)
        runs.append(run)
        paths = [full.legacy.GAUGES/('qX_p2n%d' % cfg), marker]
        paths += [Path(run['run'])/'eigenvectors'/('qX_p2n%d_evec_t%d.h5' % (cfg, t))
                  for t in (95, 96, 97, 99, 100, 101)]
        for path in paths:
            stat = path.stat()
            if stat.st_size <= 0:
                raise ValueError('Empty input: '+str(path))
            inputs[str(path)] = dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    data = dict(name=args.name, stage=args.stage, configurations=configs, cases=case_plan(args.stage),
                slots=slot_plan(args.stage), runs=runs, inputs=inputs,
                reused_cases=[dict(configuration=100, tau=4, piters=[4, 40], job_id=5942805)]
                if args.stage == 'full70' else [],
                executable_sha256=digest(pilot.EXE), build_manifest_sha256=manifest,
                source_manifest_sha256=digest(HERE/'INPUTS.json'),
                pilot_manifest_sha256=digest(pilot.HERE/'INPUTS.json'),
                settings=dict(Nv=100, midpoint=98, alpha=[1, 1, .5], hyp_steps=1,
                              grid=[1, 4, 4, 4], ranks=64, nodes=8),
                scope=('One matched 70-configuration cohort; raw and traceless at piter4/40 and tau2/4/6. '
                       'The existing raw piter10/20 results remain in the projection study.')
                if args.stage == 'full70' else
                'Selected central-window diagnostic at piter6/8/80; no full-lattice field-stability claim')
    target = folder(args.name)
    target.mkdir(parents=True, exist_ok=False)
    write_json(target/'PREPARED.json', data)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared', args.stage, len(data['cases']), 'new cases;',
          sum(len(c['piters']) for c in data['cases']), 'MPI launches;', len(data['slots']), 'array tasks')


def prepared(name):
    verify_sources()
    target = folder(name)
    if digest(target/'PREPARED.json') != (target/'prepared.sha256').read_text().strip():
        raise ValueError('Preparation changed')
    data = json.loads((target/'PREPARED.json').read_text())
    if data['cases'] != case_plan(data['stage']) or data['slots'] != slot_plan(data['stage']):
        raise ValueError('Case or slot plan changed')
    if (digest(HERE/'INPUTS.json') != data['source_manifest_sha256'] or
            digest(pilot.HERE/'INPUTS.json') != data['pilot_manifest_sha256'] or
            digest(pilot.EXE) != data['executable_sha256'] or
            pilot.build_manifest() != data['build_manifest_sha256']):
        raise ValueError('Measurement implementation changed')
    for name, expected in data['inputs'].items():
        stat = Path(name).stat()
        if (stat.st_size, stat.st_mtime_ns) != (expected['bytes'], expected['mtime_ns']):
            raise ValueError('Prepared input changed: '+name)
    return target, data


def references(entries):
    wanted = {(e['configuration'], e['tau']) for e in entries}
    records = {}
    with REFERENCE.open(newline='') as stream:
        for row in csv.DictReader(stream):
            cfg, tau, p = (int(row[k]) for k in ('configuration', 'tau', 'piter'))
            if (cfg, tau) in wanted and p in (4, 40):
                records[(cfg, tau, p, *(int(row[k]) for k in GEOMETRY))] = row
    if len(records) != len(wanted)*664:
        raise ValueError('Incomplete preserved raw reference')
    return records


def read_case(out, entry, refs):
    cfg, tau = entry['configuration'], entry['tau']
    check = full.legacy.checker()
    raw, windows, traces = [], [], {}
    for p in entry['piters']:
        path = out/('p%d.out' % p)
        _, scan, _, window = check.read_measurement(path, cfg, tau, p)
        if p >= 40 and max(float(window[k]) for k in ('unit_max', 'det_max')) > 1e-12:
            raise ValueError('High-iteration central SU(3) control failed')
        traces[p] = pilot.trace_records(path.read_text(), scan, check.avg)
        for key, row in scan.items():
            ref = refs[(cfg, tau, p if p in (4, 40) else 40, *key)]
            fields = check.FIELDS if p in (4, 40) else ('S_RAW',)
            for field in fields:
                check.avg.near(float(row[field]), float(ref[field]), 'preserved raw regression '+field)
        for row in check.raw_rows(cfg, tau, p, scan):
            key = tuple(int(row[k]) for k in GEOMETRY)
            row.update({k: traces[p][key][k] for k in pilot.TRACE_FIELDS})
            raw.append({k: str(v) for k, v in row.items()})
        windows.append({k: str(v) for k, v in window.items()})
    for key in traces[entry['piters'][0]]:
        if len({Decimal(traces[p][key]['T_TRACE']) for p in entry['piters']}) != 1:
            raise ValueError('Original trace component changed across piters')
    return raw, windows


def current_attempt(target, index):
    case = target/('case_%03d' % index)
    marker = case/'LATEST.json'
    if not marker.is_file():
        return None
    record = json.loads(marker.read_text())
    name = record['attempt']
    if not re.fullmatch(r'attempt_\d{4}', name):
        raise ValueError('Invalid attempt name')
    attempt = case/name
    if digest(attempt/'SUMMARY.json') != record['summary_sha256']:
        raise ValueError('Completed attempt summary changed')
    if json.loads((attempt/'SUMMARY.json').read_text()).get('success') is not True:
        raise ValueError('Completed attempt lacks success')
    return attempt


def verify_launches(out, entry):
    launches = json.loads((out/'LAUNCHES.json').read_text())
    actual = [(r['configuration'], r['tau'], r['piter'], r['exit_code']) for r in launches]
    if actual != [(entry['configuration'], entry['tau'], p, 0) for p in entry['piters']]:
        raise ValueError('Incomplete case launch inventory')


def run(args):
    target, data = prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS', '0')) != 64 or int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 8:
        raise ValueError('Use the submitted 8-node, 64-rank allocation')
    slot = int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0 <= slot < len(data['slots']):
        raise ValueError('Array task outside slot plan')
    entries = [data['cases'][i] for i in data['slots'][slot]]
    refs = references(entries)
    for entry in entries:
        index, cfg, tau = entry['index'], entry['configuration'], entry['tau']
        previous = current_attempt(target, index)
        if previous:
            verify_launches(previous, entry)
            read_case(previous, entry, refs)
            print('REUSED completed case', index)
            continue
        case = target/('case_%03d' % index)
        case.mkdir(exist_ok=True)
        numbers = [int(p.name.split('_')[1]) for p in case.glob('attempt_*') if re.fullmatch(r'attempt_\d{4}', p.name)]
        out = case/('attempt_%04d' % (max(numbers, default=0)+1))
        out.mkdir(exist_ok=False)
        write_json(out/'CASE.json', dict(entry, prepared_sha256=digest(target/'PREPARED.json'),
                                       job_id=os.environ['SLURM_JOB_ID']))
        run_path = next(r['run'] for r in data['runs'] if r['configuration'] == cfg)
        launches = []
        try:
            for p in entry['piters']:
                env = os.environ.copy()
                env.update(QX_CFG=str(cfg), QX_ROOT=str(full.legacy.GAUGES),
                           EVROOT=str(Path(run_path)/'eigenvectors'), T0=str(98-tau//2), TAU_LENGTH=str(tau),
                           HYP_MODE='2', HYP_PITER=str(p), PRINT_ALL='1', GX='0', GY='0', GZ='0')
                start = time.monotonic()
                with (out/('p%d.out' % p)).open('x') as stdout, (out/('p%d.err' % p)).open('x') as stderr:
                    code = subprocess.run(['mpiexec', '-n', '64', str(pilot.EXE)], env=env,
                                          stdout=stdout, stderr=stderr, check=False).returncode
                launches.append(dict(configuration=cfg, tau=tau, piter=p, exit_code=code,
                                     seconds=time.monotonic()-start, job_id=os.environ['SLURM_JOB_ID']))
                write_json(out/('launch_p%d.json' % p), launches[-1])
                if code:
                    raise ValueError('MPI launch failed for piter%d' % p)
            write_json(out/'LAUNCHES.json', launches)
            raw, windows = read_case(out, entry, refs)
            check = full.legacy.checker()
            check.write_csv(out/'profiles_raw_and_trace.csv', check.RAW_COLUMNS+list(pilot.TRACE_FIELDS), raw)
            check.write_csv(out/'su3_windows.csv', check.SU3_COLUMNS, windows)
            write_json(out/'SUMMARY.json', dict(success=True, index=index, configuration=cfg, tau=tau,
                                              piters=entry['piters'], raw_records=len(raw), window_records=len(windows)))
            marker = case/('LATEST_'+out.name+'.tmp')
            write_json(marker, dict(attempt=out.name, summary_sha256=digest(out/'SUMMARY.json')))
            marker.replace(case/'LATEST.json')
        except Exception as exc:
            write_json(out/'FAILED.json', dict(success=False, error=str(exc), launches=launches))
            raise
        print('PASS case', index, 'n%d tau%d' % (cfg, tau))


def submission_records(target):
    return [json.loads(p.read_text()) for p in sorted(target.glob('SUBMISSION_*.json'))]


def submit(args):
    target, data = prepared(args.name)
    analysis.verify()
    if not 1 <= args.concurrent <= 4:
        raise ValueError('Use one to four simultaneous array tasks')
    previous = submission_records(target)
    if previous:
        if not args.resume:
            raise ValueError('Campaign has submission records; use --resume after its jobs finish')
        ids = set()
        for record in previous:
            match = re.fullmatch(r'(\d+)(?:;\S+)?\n?', record.get('stdout', ''))
            if record.get('returncode') != 0 or not match:
                raise ValueError('Ambiguous previous submission; inspect records before resubmitting')
            ids.add(match.group(1))
        live = subprocess.run(['squeue', '-h', '-u', os.environ['USER'], '-o', '%A %F'],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False)
        if live.returncode or ids.intersection(re.findall(r'\d+', live.stdout)):
            raise ValueError('Previous campaign may still be active; no duplicate submission made')
    slots = [i for i, indices in enumerate(data['slots'])
             if any(current_attempt(target, index) is None for index in indices)]
    if not slots:
        print('All cases already completed; collect the campaign')
        return
    command = ['sbatch', '--parsable', '--array='+','.join(map(str, slots))+'%'+str(args.concurrent),
               str(HERE/'submit_micro.sh'), args.name]
    record = target/('SUBMISSION_%03d.json' % (len(previous)+1))
    write_json(record, dict(state='submitting', command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=False)
    record.write_text(json.dumps(dict(command=command, returncode=result.returncode,
                                     stdout=result.stdout, stderr=result.stderr), indent=2)+'\n')
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?', result.stdout):
        raise ValueError('Submission did not return a clear job ID; inspect '+str(record))
    print('SUBMITTED', result.stdout.strip(), data['stage'], len(slots), 'array tasks')


def seed_data(stage):
    if stage != 'full70':
        return [], []
    rows, _ = analysis.validated_rows()
    raw = [{k: v for k, v in row.items() if k not in analysis.DERIVED} for row in rows]
    check = full.legacy.checker()
    windows = [check.read_measurement(analysis.INPUT/'outputs'/('p%d.out' % p), 100, 4, p)[3] for p in (4, 40)]
    return raw, [{k: str(v) for k, v in row.items()} for row in windows]


def assert_cohort(raw, windows, stage):
    configs = TARGET if stage == 'full70' else SPECTRUM_CONFIGS
    piters = [4, 40] if stage == 'full70' else [6, 8, 80]
    groups = {(c, t, p) for c in configs for t in TAUS for p in piters}
    observed = {}
    for row in raw:
        group = tuple(int(row[k]) for k in ('configuration', 'tau', 'piter'))
        key = tuple(int(row[k]) for k in GEOMETRY)
        if key in observed.setdefault(group, set()):
            raise ValueError('Duplicate raw profile key')
        observed[group].add(key)
    check = full.legacy.checker()
    expected = {(a, r, *offset) for a in (-1, 0, 1, 2) for r in (1, 2, 4)
                for offset in check.profile.offsets(r)}
    if set(observed) != groups or any(keys != expected for keys in observed.values()):
        raise ValueError('Incomplete matched raw/trace cohort')
    window_keys = [tuple(int(r[k]) for k in ('configuration', 'tau', 'piter')) for r in windows]
    if len(window_keys) != len(groups) or set(window_keys) != groups:
        raise ValueError('Incomplete or duplicate SU(3) window cohort')


def collect(args):
    target, data = prepared(args.name)
    raw, windows = seed_data(data['stage'])
    missing = []
    for indices in data['slots']:
        entries = [data['cases'][i] for i in indices]
        refs = references(entries)
        for entry in entries:
            out = current_attempt(target, entry['index'])
            if out is None:
                missing.append(entry['index'])
                continue
            verify_launches(out, entry)
            r, w = read_case(out, entry, refs)
            check = full.legacy.checker()
            for name, actual in [('profiles_raw_and_trace.csv', r), ('su3_windows.csv', w)]:
                with (out/name).open(newline='') as stream:
                    if list(csv.DictReader(stream)) != actual:
                        raise ValueError('Saved case table differs: '+str(out/name))
            raw += r
            windows += w
    if missing and not args.partial:
        raise ValueError('%d cases are incomplete; use collect --partial for a diagnostic export' % len(missing))
    if not missing:
        assert_cohort(raw, windows, data['stage'])
    check = full.legacy.checker()
    raw.sort(key=lambda r: tuple(int(r[k]) for k in ('configuration', 'tau', 'piter', *GEOMETRY)))
    windows.sort(key=lambda r: tuple(int(r[k]) for k in ('configuration', 'tau', 'piter')))
    suffix = '_partial' if missing else ''
    check.write_csv(target/('profiles_raw_and_trace'+suffix+'.csv'), check.RAW_COLUMNS+list(pilot.TRACE_FIELDS), raw)
    check.write_csv(target/('su3_windows'+suffix+'.csv'), check.SU3_COLUMNS, windows)
    report = dict(success=not missing, stage=data['stage'], target_configurations=data['configurations'],
                  missing_cases=missing, raw_records=len(raw), window_records=len(windows),
                  reused_cases=data['reused_cases'], scope=data['scope'])
    (target/('SUMMARY'+suffix+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    output = ROOT/('qx_clover_trace_%s_review%s.zip' % (args.name, suffix))
    with zipfile.ZipFile(str(output), 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(target.rglob('*')):
            if path.is_file():
                archive.write(str(path), path.relative_to(ROOT).as_posix())
        for path in (HERE/'INPUTS.json', pilot.HERE/'INPUTS.json', pilot.HERE/'build/manifest.sha256'):
            archive.write(str(path), path.relative_to(ROOT).as_posix())
    print('PASS complete matched cohort' if not missing else 'PARTIAL export; incomplete cases remain')
    print('REVIEW_ARCHIVE', output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    command = sub.add_parser('prepare')
    command.add_argument('--name', required=True)
    command.add_argument('--stage', choices=['full70', 'spectrum'], required=True)
    command = sub.add_parser('submit')
    command.add_argument('--name', required=True)
    command.add_argument('--concurrent', type=int, default=4)
    command.add_argument('--resume', action='store_true')
    command = sub.add_parser('run')
    command.add_argument('--name', required=True)
    command = sub.add_parser('collect')
    command.add_argument('--name', required=True)
    command.add_argument('--partial', action='store_true')
    args = parser.parse_args()
    if not args.action:
        parser.error('Choose prepare, submit, run or collect')
    globals()[args.action](args)


if __name__ == '__main__':
    main()
