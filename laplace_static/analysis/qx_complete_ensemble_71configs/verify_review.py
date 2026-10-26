#!/usr/bin/env python3
"""Replay the preserved n640 measurement logs without accessing HPC data files."""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def verify(base=None):
    base = Path(base) if base is not None else HERE / 'inputs/job5946008'
    campaign = base / 'campaign'
    plan = read(campaign / 'PREPARED.json')
    if digest(campaign / 'PREPARED.json') != (campaign / 'prepared.sha256').read_text().strip():
        raise ValueError('Prepared n640 plan changed')
    expected = [dict(index=i, configuration=640, tau=t, piters=[4, 10, 20, 40])
                for i, t in enumerate((2, 4, 6))]
    if plan['cases'] != expected:
        raise ValueError('Unexpected n640 measurement inventory')
    for name, sha in plan['package_files'].items():
        if digest(base / 'workflow' / name) != sha:
            raise ValueError('Archived workflow changed: ' + name)
    if digest(base / 'build_provenance/trace_manifest.sha256') != plan['build_manifest_sha256']:
        raise ValueError('Archived build manifest changed')
    recovery = read(base / 'workflow/inputs/recovery_summary.json')
    if plan['recovery'] != recovery or recovery['complete_eigenpair_set_valid'] is not True:
        raise ValueError('Measurement lacks the reviewed complete recovery receipt')
    pins = read(base / 'workflow/INPUTS.json')['source_pins']
    for name in ('laplace_static/validation/qx_clover_trace/campaign.py',
                 'laplace_static/validation/qx_projection_spectrum/campaign.py'):
        if digest(ROOT / name) != pins[name]:
            raise ValueError('Measurement parser source changed: ' + name)
    # Reuse numerical log parsers; no remote build/gauge-path checks are invoked.
    checker_dir = ROOT / 'laplace_static/validation/hyp_projection_configs_qx_200_230'
    sys.path.insert(0, str(checker_dir))
    checker = module('qx71_measurement_checker', checker_dir / 'check_task.py')
    trace = module('qx71_trace_checker', ROOT / 'laplace_static/validation/qx_clover_trace/campaign.py')
    workflow = module('qx71_archived_n640_workflow', base / 'workflow/campaign.py')
    raw, windows, launches = [], [], []
    for entry in expected:
        out = campaign / ('case_%03d' % entry['index'])
        case = read(out / 'CASE.json')
        if case != dict(entry, prepared_sha256=digest(campaign / 'PREPARED.json')):
            raise ValueError('Case does not match prepared plan')
        launches.extend(workflow.launch_records(out, entry))
        r, w = workflow.read_case(out, entry, trace, checker)
        if rows(out / 'profiles_raw_and_trace.csv') != r or rows(out / 'su3_windows.csv') != w:
            raise ValueError('Case exports differ from stdout replay')
        summary = dict(success=True, configuration=640, tau=entry['tau'],
                       piters=entry['piters'], raw_records=1328, window_records=4)
        if read(out / 'SUMMARY.json') != summary:
            raise ValueError('Case summary differs from replay')
        raw.extend(r)
        windows.extend(w)
    if (len(raw), len(windows), len(launches)) != (3984, 12, 12):
        raise ValueError('Incomplete n640 measurements')
    if rows(campaign / 'profiles_raw_and_trace.csv') != raw or rows(campaign / 'su3_windows.csv') != windows:
        raise ValueError('Joined n640 exports differ from replay')
    summary = dict(success=True, configuration=640, cases=3, piters=[4, 10, 20, 40],
                   taus=[2, 4, 6], raw_records=3984, window_records=12,
                   ensemble_71_recomputed=False, scope=plan['scope'])
    if read(campaign / 'SUMMARY.json') != summary:
        raise ValueError('Joined campaign summary differs from replay')
    return dict(success=True, job_id=5946008, configuration=640, cases=3, mpi_launches=12,
                raw_records=3984, window_records=12, output_stream_hashes=24,
                workflow_files_verified=len(plan['package_files']),
                raw_and_trace_tables_reproduced_exactly=True,
                MPI_minutes_per_tau=[dict(tau=t, minutes=sum(x['seconds'] for x in launches
                    if x['tau'] == t) / 60) for t in (2, 4, 6)],
                recovered_eigenpairs=recovery,
                scope='Offline replay of all measurement logs and archived recovery receipt. '
                      'The physical eigenvector read-back audit was completed on SuperMUC before measurement.')


if __name__ == '__main__':
    result = verify()
    print('PASS n640 replay: 12 stdout logs; 24 stream hashes; 3984 primary rows; 12 windows; trace identities and MPI checks')
