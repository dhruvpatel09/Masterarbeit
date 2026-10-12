"""Submission bookkeeping; reuse the original scientific checks unchanged."""
import hashlib
import json
import math
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORIGINAL = HERE.with_name('hyp_projection_configs_qx')
sys.path.insert(0, str(ORIGINAL))
import common as c
import check_task as checks

LEDGER = HERE / 'RUNS.json'
PROVENANCE = ('package.sha256', 'prepared.sha256', 'build/manifest.sha256')


def now():
    return datetime.now().astimezone().isoformat(timespec='seconds')


def say(event, **details):
    print(now(), event, json.dumps(details, allow_nan=False), flush=True)


def atomic_json(path, value):
    """Flush bookkeeping before/after each submission; never overwrite outputs."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix=path.name+'.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory = os.open(str(path.parent), os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def verify_environment():
    c.verify_manifest(HERE / 'package.sha256')
    c.verify_manifest(HERE / 'inputs.sha256')
    c.prepared()


def provenance():
    return {name: c.digest(c.CASE / name) for name in PROVENANCE}


def validate_ledger(data):
    if data.get('schema') != 1 or data.get('measurement_directory') != str(c.CASE):
        raise RuntimeError('Ledger format or measurement directory differs')
    c.array_id(data['first_job'])
    if data.get('provenance') != provenance():
        raise RuntimeError('Build/package/preparation changed since controller start')
    entries = data.get('cases', [])
    if len(entries) != 15:
        raise RuntimeError('Ledger must describe exactly 15 cases')
    for i, entry in enumerate(entries):
        if any(entry.get(k) != v for k, v in c.case_for(i).items()):
            raise RuntimeError('Ledger case mapping differs at index '+str(i))
        if entry.get('state') not in ('planned', 'submitting', 'submitted', 'completed'):
            raise RuntimeError('Invalid controller state at index '+str(i))
        if entry.get('array_job_id') is not None:
            c.array_id(entry['array_job_id'])
        if entry['state'] in ('planned', 'submitting') and entry.get('array_job_id'):
            raise RuntimeError('Unsubmitted case already has a job ID')
        if entry['state'] in ('submitted', 'completed') and not entry.get('array_job_id'):
            raise RuntimeError('Submitted case has no recorded job ID')
    if entries[0].get('array_job_id') != data['first_job']:
        raise RuntimeError('First-job mapping changed')
    jobs = [entry['array_job_id'] for entry in entries if entry.get('array_job_id')]
    if len(jobs) != len(set(jobs)):
        raise RuntimeError('Expected one distinct single-element submission per case')
    return data


def load_ledger():
    return validate_ledger(json.loads(LEDGER.read_text()))


def output_directory(entry):
    job = c.array_id(entry['array_job_id'])
    return c.CASE / f"results_{job}_{entry['index']}"


def validate_case(entry, recheck=True):
    """Require case identity, matching provenance and successful actual launches."""
    index, job = entry['index'], c.array_id(entry['array_job_id'])
    selected = c.case_for(index)
    folder = output_directory(entry)
    if (folder / 'FAILED.json').exists():
        raise RuntimeError(f'Measurement failed: {folder / "FAILED.json"}')
    if json.loads((folder / 'CASE.json').read_text()) != dict(selected, array_job_id=job):
        raise RuntimeError('Wrong case identity: '+str(folder))
    for name in PROVENANCE:
        if (folder / name.replace('/', '_')).read_bytes() != (c.CASE / name).read_bytes():
            raise RuntimeError('Measurement provenance changed: '+str(folder))
    launches = json.loads((folder / 'LAUNCHES.json').read_text())
    expected = [(selected['configuration'], selected['tau'], p, f'p{p}') for p in (4,40)]
    if index == 0:
        expected += [(100,4,p,f'n100_p{p}') for p in (4,40)]
    actual = [(r['configuration'],r['tau'],r['piter'],r['label']) for r in launches]
    if actual != expected or any(r['exit_code'] != 0 for r in launches):
        raise RuntimeError('MPI launch inventory or return status failed: '+str(folder))
    if any(not math.isfinite(r['seconds']) or r['seconds'] < 0 for r in launches):
        raise RuntimeError('Invalid launch timing: '+str(folder))
    saved = json.loads((folder / 'SUMMARY.json').read_text())
    if saved.get('success') is not True:
        raise RuntimeError('Successful measurement marker absent: '+str(folder))
    if (saved.get('configuration'), saved.get('tau')) != (selected['configuration'], selected['tau']):
        raise RuntimeError('Wrong summary identity: '+str(folder))
    if index == 0 and saved.get('n100_regression',{}).get('success') is not True:
        raise RuntimeError('The first case must include the successful n100 regression')
    summary = checks.check(folder, selected['configuration'], selected['tau'], regression=(index==0)) if recheck else saved
    if summary.get('success') is not True:
        raise RuntimeError('Output recheck failed: '+str(folder))
    summary.update(index=index, array_job_id=job,
                   mpi_seconds=sum(r['seconds'] for r in launches))
    return summary
