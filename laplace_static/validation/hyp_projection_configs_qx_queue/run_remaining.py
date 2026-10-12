#!/usr/bin/env python3
"""Submit indices 1..14 one at a time, then collect all distinct job IDs."""
import argparse
import fcntl
import json
import os
import re
import subprocess
import sys
import time
import queue_common as q


def parse_job_id(stdout):
    match = re.fullmatch(r'\s*(\d+)(?:;[A-Za-z0-9_.-]+)?\s*', stdout)
    if not match:
        raise RuntimeError('Unrecognized sbatch reply; inspect RUNS.json before any new submission')
    return match.group(1)


def slurm_environment():
    # Command-line array selection overrides the original 0-14 array directive.
    # Prevent ambient SBATCH_* or SQUEUE_* variables from changing this workflow.
    return {k:v for k,v in os.environ.items() if not k.startswith(('SBATCH_', 'SQUEUE_'))}


def active_job(entry):
    # Query the user's inventory: selecting a finished ID directly can produce
    # "Invalid job id". A failed inventory query never means job completion.
    result = subprocess.run(
        ['squeue','--noheader','--array','--user',str(os.getuid()),'--format=%i|%T'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        timeout=30, env=slurm_environment())
    if result.returncode:
        raise RuntimeError('squeue failed; controller stopped with job ID retained: '+result.stderr.strip())
    wanted = {entry['array_job_id'], f"{entry['array_job_id']}_{entry['index']}"}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.strip().split('|')
        if len(parts) != 2:
            raise RuntimeError('Unrecognized squeue inventory; completion not assumed')
        if parts[0].strip() in wanted:
            return parts[1].strip()
    return None


def wait_for_completion(entry, poll_seconds):
    last_state = object()
    while True:
        state = active_job(entry)
        if state is None:
            # Never infer success from queue disappearance. Recheck real outputs.
            return q.validate_case(entry)
        if state != last_state:
            q.say('WAIT', index=entry['index'], job=entry['array_job_id'], state=state)
            last_state = state
        time.sleep(poll_seconds)


def submit_one(data, entry):
    entry.pop('sbatch_reply', None)
    entry.update(state='submitting', submission_started=q.now())
    q.atomic_json(q.LEDGER, data)
    command = ['sbatch','--parsable',f"--array={entry['index']}",str(q.c.CASE / 'submit.sh')]
    # A timeout/exception leaves the intent recorded. Never blindly resubmit:
    # Slurm might have accepted the job before the reply was interrupted.
    result = subprocess.run(command, cwd=str(q.c.ROOT), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True,
                            timeout=60, env=slurm_environment())
    entry['sbatch_reply'] = dict(returncode=result.returncode,
                                stdout=result.stdout, stderr=result.stderr)
    q.atomic_json(q.LEDGER, data)
    if result.returncode:
        # Only the explicit submit-limit rejection establishes nonacceptance.
        if 'QOSMaxSubmitJobPerUserLimit' in result.stderr:
            entry['state'] = 'planned'
            q.atomic_json(q.LEDGER, data)
            raise RuntimeError('QoS submit limit rejected this case. No job ID was accepted; '
                               'resume later with --resume after other submissions finish.')
        raise RuntimeError('sbatch failed; recorded reply requires inspection before retry: '+result.stderr.strip())
    job = parse_job_id(result.stdout)
    if any(e.get('array_job_id') == job for e in data['cases']):
        raise RuntimeError('sbatch returned a previously recorded parent job ID')
    entry.update(state='submitted', array_job_id=job, submitted=q.now())
    q.atomic_json(q.LEDGER, data)
    q.say('SUBMITTED', index=entry['index'], configuration=entry['configuration'],
          tau=entry['tau'], job=job)


def initialize(first_job):
    q.c.array_id(first_job)
    data = dict(schema=1, first_job=first_job, created=q.now(),
                measurement_directory=str(q.c.CASE), provenance=q.provenance(),
                cases=[dict(q.c.case_for(i), state='planned', array_job_id=None) for i in range(15)])
    first = data['cases'][0]
    first.update(state='completed', array_job_id=first_job)
    q.validate_case(first)
    first['validated'] = q.now()
    q.atomic_json(q.LEDGER, data)
    return data


def recover_reply(data, entry):
    reply = entry.get('sbatch_reply', {})
    if reply.get('returncode') == 0:
        job = parse_job_id(reply.get('stdout',''))
        if any(e.get('array_job_id') == job for e in data['cases']):
            raise RuntimeError('Recorded sbatch job ID is already mapped to another case')
        entry.update(state='submitted', array_job_id=job)
        q.atomic_json(q.LEDGER, data)
        q.say('RECOVERED_SUBMISSION', index=entry['index'], job=job)
        return
    if reply.get('returncode') != 0 and 'QOSMaxSubmitJobPerUserLimit' in reply.get('stderr',''):
        entry['state'] = 'planned'
        q.atomic_json(q.LEDGER, data)
        return
    raise RuntimeError('Interrupted/ambiguous submission at index '+str(entry['index'])+
                       '. Inspect RUNS.json and squeue; do not start a second campaign or resubmit blindly.')


def run(data, poll_seconds):
    for entry in data['cases']:
        if entry['state'] == 'submitting':
            recover_reply(data, entry)
        if entry['state'] == 'planned':
            # Recheck inputs/build before every future allocation.
            q.verify_environment()
            submit_one(data, entry)
        if entry['state'] == 'submitted':
            summary = wait_for_completion(entry, poll_seconds)
        else:
            summary = q.validate_case(entry)
        entry.update(state='completed', validated=q.now(), mpi_seconds=summary['mpi_seconds'])
        q.atomic_json(q.LEDGER, data)
        q.say('VALIDATED', index=entry['index'], configuration=entry['configuration'],
              tau=entry['tau'], job=entry['array_job_id'])
    data['completed'] = q.now()
    q.atomic_json(q.LEDGER, data)
    import collect_mixed
    archive = collect_mixed.collect(data)
    data['review_archive'] = str(archive)
    q.atomic_json(q.LEDGER, data)
    q.say('COMPLETE', cases=15, review_archive=str(archive))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first-job', help='Already successful parent job for index 0, e.g. 5918419')
    parser.add_argument('--resume', action='store_true', help='Continue the existing RUNS.json without duplicating submissions')
    parser.add_argument('--poll-seconds', type=int, default=30)
    args = parser.parse_args()
    if args.poll_seconds < 15 or args.poll_seconds > 300:
        parser.error('poll-seconds must be 15..300')
    if args.resume == bool(args.first_job):
        parser.error('Choose --first-job JOB or --resume')
    with (q.HERE / 'controller.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another controller is already running')
        q.verify_environment()
        if args.resume:
            data = q.load_ledger()
        else:
            if q.LEDGER.exists():
                raise RuntimeError('RUNS.json already exists; use --resume')
            data = initialize(args.first_job)
        q.say('CONTROLLER_START', resume=args.resume, pending_limit=1)
        run(data, args.poll_seconds)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        q.say('STOPPED', error=str(exc))
        sys.exit(1)
