#!/usr/bin/env python3
"""Launch two prescriptions on one configuration/tau; retain failed runs."""
import json
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path
import common as c
import check_task

def launch(outdir, cfg, tau, piter, run, label):
    env = os.environ.copy()
    env.update(QX_CFG=str(cfg), QX_ROOT=str(c.GAUGE_ROOT),
               EVROOT=str(Path(run) / 'eigenvectors'), T0=str(98-tau//2),
               TAU_LENGTH=str(tau), HYP_MODE='2', HYP_PITER=str(piter),
               PRINT_ALL='1', GX='0', GY='0', GZ='0')
    # Use the launcher form from the successful n100 profile jobs.
    command = ['mpiexec', '-n', '64',
               str(c.CASE / 'build/test_profiles_qx')]
    start = time.monotonic()
    with (outdir / (label+'.out')).open('x') as stdout, (outdir / (label+'.err')).open('x') as stderr:
        try:
            result = subprocess.run(command, env=env, stdout=stdout, stderr=stderr, check=False)
            code = result.returncode
        except OSError as exc:
            stderr.write(str(exc)+'\n')
            code = 127
    record = dict(configuration=cfg, tau=tau, piter=piter, label=label,
                  exit_code=code, seconds=time.monotonic()-start, command=command)
    print('MEASUREMENT_DONE', json.dumps(record), flush=True)
    return record

def main():
    data = c.prepared()
    index = int(os.environ['SLURM_ARRAY_TASK_ID'])
    job = c.array_id(os.environ['SLURM_ARRAY_JOB_ID'])
    if int(os.environ['SLURM_NTASKS']) != 64:
        raise RuntimeError('Require the prepared 64-rank observable layout')
    selected = c.case_for(index)
    cfg, tau = selected['configuration'], selected['tau']
    run = next(x['run'] for x in data['runs'] if x['configuration']==cfg)
    outdir = c.CASE / f'results_{job}_{index}'
    outdir.mkdir()
    c.write_json(outdir / 'CASE.json', dict(selected, array_job_id=job))
    print('START', datetime.now().astimezone().isoformat(), 'CASE', json.dumps(selected), flush=True)
    # Save the exact active manifests with each measurement pair.
    for source in ('package.sha256', 'prepared.sha256', 'build/manifest.sha256'):
        name = source.replace('/', '_')
        (outdir/name).write_bytes((c.CASE/source).read_bytes())
    records = [launch(outdir,cfg,tau,p,run,f'p{p}') for p in (4,40)]
    if index == 0:
        records += [launch(outdir,100,4,p,data['n100_run'],f'n100_p{p}') for p in (4,40)]
    c.write_json(outdir / 'LAUNCHES.json', records)
    try:
        if any(r['exit_code'] != 0 for r in records):
            raise RuntimeError('An MPI measurement returned a nonzero exit code; see LAUNCHES.json')
        summary = check_task.check(outdir,cfg,tau,regression=(index==0))
        c.write_json(outdir / 'SUMMARY.json', summary)
    except Exception as exc:
        c.write_json(outdir / 'FAILED.json', dict(success=False,error=str(exc)))
        raise
    print('PASS paired profile checks; piter4 acceptance remains a physics decision', flush=True)
    print('END', datetime.now().astimezone().isoformat(), flush=True)

if __name__ == '__main__':
    main()
