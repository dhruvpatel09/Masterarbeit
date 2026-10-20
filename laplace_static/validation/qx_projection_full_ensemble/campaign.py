#!/usr/bin/env python3
"""Complete one matched qX ensemble, reusing the preserved 29 configurations."""
import argparse
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import tempfile
import time
import zipfile
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = [n for n in range(100, 801, 10) if n != 640]
BASELINE = list(range(110, 231, 10))+list(range(650, 801, 10))
MISSING = sorted(set(TARGET)-set(BASELINE))
PITERS = [4, 10, 20, 40]
TAUS = [2, 4, 6]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


legacy = module('qx_full_legacy', ROOT/'laplace_static/validation/qx_projection_spectrum/campaign.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def verify_sources():
    legacy.verify_sources()
    pins = json.loads((HERE/'INPUTS.json').read_text())['files']
    for name, expected in pins.items():
        if digest(ROOT/name) != expected:
            raise ValueError('Full-ensemble source changed: '+name)


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}', name):
        raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def case_plan():
    return [dict(index=i, configuration=cfg, tau=tau, piters=PITERS[:])
            for i, (cfg, tau) in enumerate((cfg, tau) for cfg in MISSING for tau in TAUS)]


def slot_plan(cases, per_task=4):
    if per_task not in (1, 2, 3, 4):
        raise ValueError('Use one to four cases per array task')
    slots = [list(range(start, min(start+per_task, len(cases))))
             for start in range(0, len(cases), per_task)]
    if len(slots) > 39:
        raise ValueError('Use enough cases per task to keep the array below 40 submitted tasks')
    return slots


def inventory():
    verify_sources()
    known = dict(legacy.choose_runs())
    # Older saved-mode audits may use names predating qx_saved_modes*/results_*.
    # Admit only the same selected time coverage and a successful numerical audit.
    summaries = sorted((ROOT/'mental/validation').glob('**/SUMMARY.json'))
    targeted = ROOT/'laplace_static/validation/qx_saved_modes_n100/results_n100/SUMMARY.json'
    if targeted.is_file():
        summaries.append(targeted)
    for summary in summaries:
        prep = summary.parent/'PREPARED.json'
        if not prep.is_file():
            continue
        metadata = json.loads(prep.read_text())
        result = json.loads(summary.read_text())
        if summary == targeted:
            if ((summary.parent/'FAILED.json').exists() or
                    result.get('n650_layout_control_passed') is not True or
                    digest(prep) != (summary.parent/'prepared.sha256').read_text().strip()):
                continue
        if (metadata.get('Nv') != 100 or result.get('success') is not True or
                not set(range(95, 102)).issubset(metadata.get('times', []))):
            continue
        passed = {int(r['configuration']) for r in result.get('runs', [])
                  if r.get('success') is True and r.get('mpi_exit_code') == 0 and
                  0 <= r.get('maximum_relative_residual', float('inf')) <= 1.1e-9}
        for row in metadata.get('runs', []):
            cfg = int(row['configuration'])
            if cfg in passed and cfg not in known and row.get('completion_marker_sha256'):
                known[cfg] = {k: row[k] for k in ('configuration', 'run', 'completion_marker_sha256')}
    missing_audits = [cfg for cfg in MISSING if cfg not in known]
    return dict(target_configurations=TARGET, reused_configurations=BASELINE,
                configurations_to_measure=MISSING, cases_to_measure=len(case_plan()),
                mpi_launches_to_measure=len(case_plan())*len(PITERS),
                array_tasks=len(slot_plan(case_plan())),
                admitted_pending_configurations=[cfg for cfg in MISSING if cfg in known],
                missing_selected_central_time_audits=missing_audits,
                ready=not missing_audits), known


def inspect(args):
    report, _ = inventory()
    print(json.dumps(report, indent=2))
    print('READY for full-ensemble preparation' if report['ready'] else
          'AUDIT_REQUIRED: '+','.join(map(str, report['missing_selected_central_time_audits'])))


def prepare(args):
    report, known = inventory()
    if not report['ready']:
        print(json.dumps(report, indent=2))
        raise ValueError('Selected central-time saved-mode audits are missing; no campaign was prepared')
    base = module('qx_full_baseline29', ROOT/'laplace_static/analysis/qx_projection_spectrum_29configs/make_plots.py')
    base.load_data()
    manifest = legacy.build_manifest()
    inputs, runs = {}, []
    for cfg in MISSING:
        run = known[cfg]
        marker = Path(run['run'])/'COMPLETED.json'
        if digest(marker) != run['completion_marker_sha256']:
            raise ValueError('Audited completion marker changed for n%d' % cfg)
        done = json.loads(marker.read_text())
        if (done['configuration'], done['timeslices'], done['eigenvector_files']) != (cfg, 192, 192):
            raise ValueError('Incomplete eigenpair inventory for n%d' % cfg)
        runs.append(run)
        paths = [legacy.GAUGES/('qX_p2n%d' % cfg), marker]
        paths += [Path(run['run'])/'eigenvectors'/('qX_p2n%d_evec_t%d.h5' % (cfg, t))
                  for t in (95, 96, 97, 99, 100, 101)]
        for path in paths:
            st = path.stat()
            if st.st_size <= 0:
                raise ValueError('Empty input: '+str(path))
            inputs[str(path)] = dict(bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    cases = case_plan()
    slots = slot_plan(cases)
    data = dict(name=args.name, stage='full70', target_configurations=TARGET,
                reused_configurations=BASELINE, cases=cases, slots=slots, runs=runs,
                inputs=inputs, source_manifest_sha256=digest(HERE/'INPUTS.json'),
                executable_sha256=digest(legacy.EXE), build_manifest_sha256=manifest,
                legacy_package_sha256=digest(legacy.CASE/'package.sha256'),
                settings=dict(Nv=100, midpoint=98, alpha=[1, 1, .5], hyp_steps=1,
                              grid=[1, 4, 4, 4], ranks=64, nodes=8),
                scope='All 70 generated configurations including n100; n640 is unavailable. '
                      'Reused validated profiles and new measurements form one matched ensemble.')
    target = folder(args.name)
    target.mkdir(parents=True, exist_ok=False)
    write_json(target/'PREPARED.json', data)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared full70:', len(BASELINE), 'reused configurations;', len(MISSING),
          'to measure;', len(cases), 'cases;', len(slots), 'array tasks')


def prepared(name):
    verify_sources()
    path = folder(name)
    if digest(path/'PREPARED.json') != (path/'prepared.sha256').read_text().strip():
        raise ValueError('Preparation changed')
    data = json.loads((path/'PREPARED.json').read_text())
    if data['cases'] != case_plan() or data['slots'] != slot_plan(case_plan()):
        raise ValueError('Full-ensemble plan changed')
    if (data['source_manifest_sha256'] != digest(HERE/'INPUTS.json') or
            data['legacy_package_sha256'] != digest(legacy.CASE/'package.sha256') or
            data['executable_sha256'] != digest(legacy.EXE) or
            data['build_manifest_sha256'] != legacy.build_manifest()):
        raise ValueError('Measurement implementation changed')
    for name, expected in data['inputs'].items():
        st = Path(name).stat()
        if (st.st_size, st.st_mtime_ns) != (expected['bytes'], expected['mtime_ns']):
            raise ValueError('Prepared input changed: '+name)
    return path, data


def submit(args):
    path, data = prepared(args.name)
    if not 1 <= args.concurrent <= 4:
        raise ValueError('Use one to four simultaneous array tasks')
    legacy.pilot_gate(args.pilot)
    if any(path.glob('case_*')):
        raise ValueError('Campaign already has case outputs; inspect before submitting')
    command = ['sbatch', '--parsable', '--array=0-%d%%%d' % (len(data['slots'])-1, args.concurrent),
               str(HERE/'submit_micro.sh'), args.name]
    record = path/'SUBMISSION.json'
    write_json(record, dict(state='submitting', command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=False)
    record.write_text(json.dumps(dict(command=command, returncode=result.returncode,
                                     stdout=result.stdout, stderr=result.stderr), indent=2)+'\n')
    if result.returncode:
        raise ValueError('Slurm rejected submission: '+result.stderr.strip())
    if not re.fullmatch(r'\d+(?:;\S+)?\n?', result.stdout):
        raise ValueError('Ambiguous submission response; inspect SUBMISSION.json before retrying')
    print('SUBMITTED', result.stdout.strip(), 'full70;', len(data['slots']), 'array tasks')


def read_case(out, entry, data):
    if entry['configuration'] != 100:
        return legacy.read_case(out, entry, data)
    check = legacy.checker()
    scans, raw, windows = {}, [], []
    for p in PITERS:
        _, records, _, window = check.read_measurement(out/('p%d.out' % p), 100, entry['tau'], p)
        scans[p] = records
        raw.extend(check.raw_rows(100, entry['tau'], p, records))
        windows.append(window)
    for key in scans[40]:
        if len({Decimal(scans[p][key]['S_RAW']) for p in PITERS}) != 1:
            raise ValueError('Original clover changed across n100 prescriptions')
    if entry['tau'] == 4:
        check.regression_check({p: scans[p] for p in (4, 40)})
    check.write_csv(out/'profiles_raw.csv', check.RAW_COLUMNS, raw)
    check.write_csv(out/'su3_windows.csv', check.SU3_COLUMNS, windows)
    return dict(success=True, configuration=100, tau=entry['tau'], piters=PITERS,
                raw_records=len(raw), window_records=len(windows),
                scope='Matched n100 measurement; tau4 also reproduces the preserved n100 regression.')


def validate_case(out, entry, data, prepared_hash):
    if (out/'FAILED.json').exists():
        raise ValueError('Case has a failure record: '+str(out))
    identity = json.loads((out/'CASE.json').read_text())
    if any(identity[k] != value for k, value in entry.items()) or identity['prepared_sha256'] != prepared_hash:
        raise ValueError('Wrong case identity')
    launches = json.loads((out/'LAUNCHES.json').read_text())
    wanted = [(entry['configuration'], entry['tau'], p, 'p'+str(p)) for p in PITERS]
    if [(r['configuration'], r['tau'], r['piter'], r['label']) for r in launches] != wanted or any(r['exit_code'] for r in launches):
        raise ValueError('Incomplete/failed launch inventory')
    result = read_case(out, entry, data)
    if result != json.loads((out/'SUMMARY.json').read_text()):
        raise ValueError('Case summary differs from measurements')
    return result


def run(args):
    path, data = prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS', '0')) != 64 or int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 8:
        raise ValueError('Run in the prepared 8-node, 64-rank Slurm allocation')
    slot = int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0 <= slot < len(data['slots']):
        raise ValueError('Array task is outside the prepared slots')
    phash = digest(path/'PREPARED.json')
    for index in data['slots'][slot]:
        entry = data['cases'][index]
        out = path/('case_%03d' % index)
        out.mkdir(exist_ok=False)
        write_json(out/'CASE.json', dict(entry, prepared_sha256=phash,
                                       job_id=os.environ['SLURM_JOB_ID'], array_task=slot))
        chosen = next(r['run'] for r in data['runs'] if r['configuration'] == entry['configuration'])
        launches = []
        try:
            for p in PITERS:
                env = os.environ.copy()
                env.update(QX_CFG=str(entry['configuration']), QX_ROOT=str(legacy.GAUGES),
                           EVROOT=str(Path(chosen)/'eigenvectors'), T0=str(98-entry['tau']//2),
                           TAU_LENGTH=str(entry['tau']), HYP_MODE='2', HYP_PITER=str(p),
                           PRINT_ALL='1', GX='0', GY='0', GZ='0')
                start = time.monotonic()
                with (out/('p%d.out' % p)).open('x') as stdout, (out/('p%d.err' % p)).open('x') as stderr:
                    code = subprocess.run(['mpiexec', '-n', '64', str(legacy.EXE)], env=env,
                                          stdout=stdout, stderr=stderr, check=False).returncode
                launches.append(dict(configuration=entry['configuration'], tau=entry['tau'], piter=p,
                                     label='p'+str(p), exit_code=code, seconds=time.monotonic()-start))
                write_json(out/('launch_%02d.json' % len(launches)), launches[-1])
                if code:
                    raise ValueError('MPI measurement failed')
            write_json(out/'LAUNCHES.json', launches)
            result = read_case(out, entry, data)
            write_json(out/'SUMMARY.json', result)
        except Exception as exc:
            write_json(out/'FAILED.json', dict(success=False, error=str(exc), launches=launches))
            raise
        print('PASS full70 case', index, 'n%d tau%d' % (entry['configuration'], entry['tau']))


def collect(args):
    path, data = prepared(args.name)
    phash = digest(path/'PREPARED.json')
    raw, windows, summaries = [], [], []
    for entry in data['cases']:
        out = path/('case_%03d' % entry['index'])
        summaries.append(validate_case(out, entry, data, phash))
        raw.extend(legacy.csv_rows(out/'profiles_raw.csv'))
        windows.extend(legacy.csv_rows(out/'su3_windows.csv'))
    check = legacy.checker()
    collection = path/'collection'
    collection.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as name:
        temp = Path(name)
        check.write_csv(temp/'profiles_raw.csv', check.RAW_COLUMNS, raw)
        check.write_csv(temp/'su3_windows.csv', check.SU3_COLUMNS, windows)
        write_json(temp/'SUMMARY.json', dict(success=True, stage='full70', cases=summaries,
                                            raw_records=len(raw), window_records=len(windows)))
        for source in temp.iterdir():
            target = collection/source.name
            if target.exists() and target.read_bytes() != source.read_bytes():
                raise ValueError('Existing collection differs from revalidated outputs')
            target.write_bytes(source.read_bytes())
    analysis = module('qx_full_analysis_collect', ROOT/'laplace_static/analysis/qx_projection_full_ensemble/make_plots.py')
    analysis.export_results([collection/'profiles_raw.csv'], [collection/'su3_windows.csv'],
                            path/'analysis', make_plots=False)
    archive = ROOT/('qx_projection_full_ensemble_'+args.name+'_review.zip')
    files = sorted(set([p for p in path.rglob('*') if p.is_file() and '__pycache__' not in p.parts] +
                       [p for p in HERE.iterdir() if p.is_file() and p.suffix not in ('.out', '.err')]))
    descriptor, temporary = tempfile.mkstemp(prefix='.qx_full_review_', suffix='.tmp', dir=str(ROOT))
    os.close(descriptor)
    temporary = Path(temporary)
    try:
        with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as z:
            for source in files:
                if source.is_symlink() or source.suffix in ('.h5', '.o', '.a', '.zip'):
                    raise ValueError('Unexpected review file')
                limit = 64*1024**2 if source.parent == collection and source.suffix == '.csv' else 8*1024**2
                if source.stat().st_size > limit:
                    raise ValueError('Review entry exceeds its size limit: '+str(source))
                z.write(source, source.relative_to(ROOT).as_posix())
        with zipfile.ZipFile(temporary) as z:
            if z.testzip() is not None or len(z.namelist()) != len(files):
                raise ValueError('New review ZIP failed verification')
        if archive.exists():
            with zipfile.ZipFile(archive) as old, zipfile.ZipFile(temporary) as new:
                if set(old.namelist()) != set(new.namelist()) or any(old.read(n) != new.read(n) for n in old.namelist()):
                    raise ValueError('Existing review ZIP differs; it was preserved')
        else:
            os.replace(str(temporary), str(archive))
    finally:
        if temporary.exists():
            temporary.unlink()
    print('PASS full70 collection:', len(raw), 'new records;', len(windows), 'new window rows')
    print('REVIEW_ARCHIVE', archive)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    sub.add_parser('inventory')
    for action in ('prepare', 'run', 'collect'):
        p = sub.add_parser(action)
        p.add_argument('--name', default='full70')
    p = sub.add_parser('submit')
    p.add_argument('--name', default='full70')
    p.add_argument('--pilot', default='pilot')
    p.add_argument('--concurrent', type=int, default=4)
    args = parser.parse_args()
    if args.action is None:
        parser.error('Select inventory, prepare, submit, run or collect')
    if args.action == 'inventory':
        inspect(args)
    else:
        globals()[args.action](args)


if __name__ == '__main__':
    main()
