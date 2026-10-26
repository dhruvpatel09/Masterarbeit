#!/usr/bin/env python3
"""Measure the recovered n640 configuration using the existing trace executable."""
import argparse
import csv
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

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
GEOMETRY=('axis','R','long','tr1','tr2')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def package_pins():
    return {str(p.relative_to(HERE)):digest(p) for p in sorted(HERE.rglob('*')) if p.is_file()
        and not any(x in ('campaigns','__pycache__') for x in p.relative_to(HERE).parts)
        and p.suffix not in ('.out','.err','.log','.zip')}


def dependencies():
    spec=read(HERE/'INPUTS.json')
    if spec['configuration']!=640 or spec['taus']!=[2,4,6] or spec['piters']!=[4,10,20,40]:
        raise ValueError('Measurement scope changed')
    for table,root in [('source_pins',ROOT),('input_pins',HERE)]:
        for name,sha in spec[table].items():
            if digest(root/name)!=sha:raise ValueError('Source/input changed: '+name)
    trace=module('qx_n640_trace',ROOT/'laplace_static/validation/qx_clover_trace/campaign.py')
    legacy=trace.verify_sources();manifest=trace.build_manifest();check=legacy.checker()
    return spec,trace,check,manifest


def recovery():
    spec,trace,check,manifest=dependencies()
    full=module('qx_n640_completed_modes',ROOT/'laplace_static/validation/qx_n640_recovery_full/campaign.py')
    folder,plan=full.prepared(spec['recovery_name'])
    result=full.final_replay(folder,plan)
    if result!=read(HERE/'inputs/recovery_summary.json') or not result['complete_eigenpair_set_valid']:
        raise ValueError('Completed recovery differs from reviewed result')
    return result


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def cases():return [dict(index=i,configuration=640,tau=tau,piters=[4,10,20,40]) for i,tau in enumerate((2,4,6))]


def prepare(args):
    spec,trace,check,manifest=dependencies();result=recovery()
    target=folder(args.name)
    if target.exists():raise ValueError('Campaign already exists; preserve it')
    plan=dict(name=args.name,cases=cases(),settings=spec['settings'],recovery=result,
        executable=str(trace.EXE),executable_sha256=digest(trace.EXE),build_manifest_sha256=manifest,
        package_files=package_pins(),scope=spec['scope'])
    target.mkdir(parents=True);write(target/'PREPARED.json',plan)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared n640: three tau cases; four piters each; 12 MPI launches; both clover definitions')


def prepared(name):
    spec,trace,check,manifest=dependencies();target=folder(name);plan=read(target/'PREPARED.json')
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared plan changed')
    if (plan['cases']!=cases() or plan['package_files']!=package_pins() or
        plan['executable_sha256']!=digest(trace.EXE) or plan['build_manifest_sha256']!=manifest):
        raise ValueError('Prepared measurement implementation changed')
    if recovery()!=plan['recovery']:raise ValueError('Recovered input set changed')
    return target,plan,trace,check


def submit(args):
    target,plan,trace,check=prepared(args.name)
    if (target/'SUBMISSION.json').exists():raise ValueError('Already submitted; do not duplicate the job')
    command=['sbatch','--parsable','--array=0-2%3',str(HERE/'submit_micro.sh'),args.name]
    write(target/'SUBMISSION.json',dict(state='submitting',command=command))
    r=subprocess.run(command,cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    write(target/'SUBMISSION.json',dict(command=command,returncode=r.returncode,stdout=r.stdout,stderr=r.stderr))
    if r.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?',r.stdout):raise ValueError('Submission failed or ambiguous; inspect SUBMISSION.json')
    print('SUBMITTED '+r.stdout.strip()+' n640 profiles: three tasks, at most three concurrent')


def read_case(out,entry,trace,check):
    cfg,tau=entry['configuration'],entry['tau'];raw=[];windows=[];scans={};traces={}
    for p in entry['piters']:
        _,scan,_,window=check.read_measurement(out/('p%d.out'%p),cfg,tau,p)
        if len(scan)!=332:raise ValueError('Expected 332 profile records per launch')
        trace_rows=trace.trace_records((out/('p%d.out'%p)).read_text(),scan,check.avg)
        scans[p]=scan;traces[p]=trace_rows
        for row in check.raw_rows(cfg,tau,p,scan):
            key=tuple(int(row[k]) for k in GEOMETRY)
            row.update({k:trace_rows[key][k] for k in trace.TRACE_FIELDS})
            raw.append({k:str(v) for k,v in row.items()})
        windows.append({k:str(v) for k,v in window.items()})
    first=entry['piters'][0]
    for p in entry['piters']:
        if scans[p].keys()!=scans[first].keys():raise ValueError('Projection profile geometries differ')
        for key in scans[first]:
            for field in ('S_RAW','Nsrc'):
                if Decimal(scans[p][key][field])!=Decimal(scans[first][key][field]):
                    raise ValueError('Original clover/source count changed across piters')
            if Decimal(traces[p][key]['T_TRACE'])!=Decimal(traces[first][key]['T_TRACE']):
                raise ValueError('Original trace component changed across piters')
    return raw,windows


def write_csv(path,columns,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)


def launch_records(out,entry):
    launches=read(out/'LAUNCHES.json')
    if [(x['configuration'],x['tau'],x['piter'],x['exit_code']) for x in launches]!=[(640,entry['tau'],p,0) for p in entry['piters']]:
        raise ValueError('Incomplete profile launch inventory')
    for row in launches:
        for ext,key in [('out','stdout_sha256'),('err','stderr_sha256')]:
            if digest(out/('p%d.%s'%(row['piter'],ext)))!=row[key]:raise ValueError('Profile output changed')
    return launches


def run(args):
    target,plan,trace,check=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:
        raise ValueError('Use an eight-node, 64-rank allocation')
    slot=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0<=slot<3:raise ValueError('Invalid array task')
    entry=plan['cases'][slot];out=target/('case_%03d'%slot);out.mkdir()
    write(out/'CASE.json',dict(entry,prepared_sha256=digest(target/'PREPARED.json')))
    launches=[]
    for p in entry['piters']:
        env=os.environ.copy();tau=entry['tau']
        env.update(QX_CFG='640',QX_ROOT='/hppfs/scratch/00/di98qes3/qX',
            EVROOT=plan['recovery']['eigenvector_directory'],T0=str(98-tau//2),TAU_LENGTH=str(tau),
            HYP_MODE='2',HYP_PITER=str(p),PRINT_ALL='1',GX='0',GY='0',GZ='0')
        command=['mpiexec','-n','64','-ppn','8',str(trace.EXE)];began=time.monotonic()
        with (out/('p%d.out'%p)).open('x') as stdout,(out/('p%d.err'%p)).open('x') as stderr:
            r=subprocess.run(command,cwd=str(ROOT),env=env,stdout=stdout,stderr=stderr)
        launches.append(dict(configuration=640,tau=tau,piter=p,exit_code=r.returncode,seconds=time.monotonic()-began,
            command=command,stdout_sha256=digest(out/('p%d.out'%p)),stderr_sha256=digest(out/('p%d.err'%p))))
        write(out/'LAUNCHES.json',launches)
        if r.returncode:raise ValueError('MPI launch failed for piter%d'%p)
    launch_records(out,entry);raw,windows=read_case(out,entry,trace,check)
    write_csv(out/'profiles_raw_and_trace.csv',check.RAW_COLUMNS+list(trace.TRACE_FIELDS),raw)
    write_csv(out/'su3_windows.csv',check.SU3_COLUMNS,windows)
    write(out/'SUMMARY.json',dict(success=True,configuration=640,tau=entry['tau'],piters=entry['piters'],raw_records=len(raw),window_records=len(windows)))
    print('PASS n640 tau%d: four piters; raw and traceless primaries'%entry['tau'])


def collect(args):
    target=folder(args.name)
    if not target.is_dir():raise ValueError('Campaign does not exist')
    try:
        target,plan,trace,check=prepared(args.name);raw=[];windows=[]
        for entry in plan['cases']:
            out=target/('case_%03d'%entry['index']);launch_records(out,entry)
            r,w=read_case(out,entry,trace,check)
            for name,actual in [('profiles_raw_and_trace.csv',r),('su3_windows.csv',w)]:
                with (out/name).open(newline='') as f:
                    if list(csv.DictReader(f))!=actual:raise ValueError('Saved table differs from log replay')
            raw.extend(r);windows.extend(w)
        if len(raw)!=3984 or len(windows)!=12:raise ValueError('Incomplete n640 matched cohort')
        write_csv(target/'profiles_raw_and_trace.csv',check.RAW_COLUMNS+list(trace.TRACE_FIELDS),raw)
        write_csv(target/'su3_windows.csv',check.SU3_COLUMNS,windows)
        summary=dict(success=True,configuration=640,cases=3,piters=[4,10,20,40],taus=[2,4,6],raw_records=len(raw),window_records=len(windows),
            ensemble_71_recomputed=False,scope=plan['scope'])
    except (OSError,ValueError,RuntimeError,KeyError,TypeError,IndexError) as error:
        summary=dict(success=False,error=str(error),ensemble_71_recomputed=False)
    write(target/'SUMMARY.json',summary)
    archive=ROOT/('qx_n640_profiles_'+args.name+'_review.zip')
    if archive.exists():raise ValueError('Review ZIP already exists; preserve or rename it first')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(target.rglob('*')):
            if path.is_file():z.write(path,'campaign/'+str(path.relative_to(target)))
        for name in package_pins():z.write(HERE/name,'workflow/'+name)
        manifest=ROOT/'laplace_static/validation/qx_clover_trace/build/manifest.sha256'
        if manifest.is_file():z.write(manifest,'build_provenance/trace_manifest.sha256')
    print('PASS complete n640 profile extension' if summary['success'] else 'INCOMPLETE: '+summary['error'])
    print('REVIEW_ARCHIVE '+str(archive))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','submit','run','collect','check-inputs'])
    p.add_argument('--name',default='matched');args=p.parse_args()
    if args.action=='check-inputs':dependencies();recovery();print('PASS existing measurement build and recovered n640 inputs')
    else:globals()[args.action](args)


if __name__=='__main__':main()
