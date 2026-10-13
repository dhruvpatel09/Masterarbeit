#!/usr/bin/env python3
"""Measure missing qX projection prescriptions with the validated 64-rank code."""
import argparse
import csv
import hashlib
import importlib
import json
import os
import re
import subprocess
import sys
import time
import zipfile
from decimal import Decimal
from pathlib import Path

CASE=Path(__file__).resolve().parent
ROOT=CASE.parents[2]
OLD=ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230'
EXE=OLD/'build/test_profiles_qx'
GAUGES=Path('/hppfs/scratch/00/di98qes3/qX')
RUNS=Path('/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs')
TAUS=(2,4,6)
PITERS=(4,10,20,40)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,data):
    with path.open('x') as handle:json.dump(data,handle,indent=2,allow_nan=False);handle.write('\n')


def csv_rows(path):
    with Path(path).open(newline='') as handle:return list(csv.DictReader(handle))


def verify_sources():
    for line in (CASE/'package.sha256').read_text().splitlines():
        sha,name=line.split(maxsplit=1)
        if digest(CASE/name)!=sha:raise ValueError('Measurement package changed: '+name)
    pins=json.loads((CASE/'INPUTS.json').read_text())['files']
    for name,sha in pins.items():
        if digest(ROOT/name)!=sha:raise ValueError('Scientific source changed: '+name)


def checker():
    verify_sources()
    sys.path.insert(0,str(OLD))
    return importlib.import_module('check_task')


def reference_rows():
    paths=json.loads((ROOT/'laplace_static/analysis/qx_hyp_projection/inputs/SOURCES.json').read_text())['default_profile_tables']
    rows=[r for p in paths for r in csv_rows(ROOT/p['path'])]
    return {(int(r['configuration']),int(r['tau']),int(r['piter']),
             *(int(r[k]) for k in ('axis','R','long','tr1','tr2'))):r for r in rows}


def choose_runs():
    chosen={}
    paths=[ROOT/'laplace_static/validation'/name/'audited_runs.json' for name in
           ('hyp_projection_configs_qx','hyp_projection_configs_qx_160_190','hyp_projection_configs_qx_200_230')]
    paths.append(CASE/'new_audited_runs.json')
    for p in paths:
        for row in json.loads(p.read_text())['runs']:
            cfg=int(row['configuration'])
            if cfg in chosen and chosen[cfg]!=row:raise ValueError('Conflicting audited run')
            chosen[cfg]=row
    # Earlier completion audits remain on SuperMUC, outside the Git repository.
    audit_root=ROOT/'mental/validation'
    for p in sorted(audit_root.glob('qx_saved_modes*/results_*/SUMMARY.json')):
        data=json.loads(p.read_text())
        prep=p.parent/'PREPARED.json'
        if data.get('success') is not True or not prep.exists():continue
        metadata=json.loads(prep.read_text())
        if not set(range(95,102)).issubset(metadata.get('times',[])):continue
        successful={int(r['configuration']) for r in data.get('runs',[])
                    if r.get('success') is True and r.get('mpi_exit_code')==0}
        for row in metadata.get('runs',[]):
            cfg=int(row['configuration'])
            if cfg in chosen or cfg not in successful:continue
            if row.get('completion_marker_sha256'):
                chosen[cfg]={k:row[k] for k in ('configuration','run','completion_marker_sha256')}
    return chosen


def build_manifest():
    path=OLD/'build/manifest.sha256'
    for line in path.read_text().splitlines():
        sha,name=line.split(maxsplit=1);p=Path(name)
        if not p.is_absolute():p=ROOT/p
        if digest(p)!=sha:raise ValueError('Existing observable build changed: '+name)
    if not EXE.is_file():raise ValueError('Validated observable executable is absent')
    return digest(path)


def campaign(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):raise ValueError('Invalid campaign name')
    return CASE/'campaigns'/name


def configs(text):
    selected=[int(s) for s in text.split(',')]
    if len(selected)!=len(set(selected)) or any(c<110 or c>800 or c%10 or c==640 for c in selected):
        raise ValueError('Use distinct n110..800 configurations, excluding n640 and n100')
    return sorted(selected)


def case_plan(selected,stage):
    if stage=='pilot' and len(selected)!=1:
        raise ValueError('Pilot requires exactly one configuration')
    if stage=='extended' and any(c<=230 for c in selected):
        raise ValueError('Extended stage adds new configurations after n230')
    taus=(4,) if stage=='pilot' else TAUS
    if stage=='intermediate13' and selected!=list(range(110,231,10)):
        raise ValueError('Intermediate stage must match the complete archived 13-configuration cohort')
    return [dict(index=i,configuration=c,tau=t,piters=list((10,20) if stage=='intermediate13' else PITERS))
            for i,(c,t) in enumerate((c,t) for c in selected for t in taus)]


def prepare(args):
    verify_sources();manifest=build_manifest()
    selected=configs(args.configs);plan=case_plan(selected,args.stage);known=choose_runs()
    inputs={};runs=[]
    for cfg in selected:
        if cfg not in known:raise ValueError(f'No selected central-time saved-mode audit for n{cfg}; inspect its audit before admission')
        run=known[cfg];marker=Path(run['run'])/'COMPLETED.json'
        if digest(marker)!=run['completion_marker_sha256']:raise ValueError('Audited completion marker changed')
        done=json.loads(marker.read_text())
        if (done['configuration'],done['timeslices'],done['eigenvector_files'])!=(cfg,192,192):raise ValueError('Incomplete generation inventory')
        runs.append(run)
        paths=[GAUGES/f'qX_p2n{cfg}',marker]
        times={98-t//2 for t in ((4,) if args.stage=='pilot' else TAUS)}|{98+t//2 for t in ((4,) if args.stage=='pilot' else TAUS)}
        paths += [Path(run['run'])/'eigenvectors'/f'qX_p2n{cfg}_evec_t{t}.h5' for t in sorted(times)]
        for p in paths:
            st=p.stat()
            if st.st_size<=0:raise ValueError('Empty measurement input: '+str(p))
            inputs[str(p)]=dict(bytes=st.st_size,mtime_ns=st.st_mtime_ns)
    n100=RUNS/'n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4'
    if args.stage=='pilot':
        for p in [GAUGES/'qX_p2n100']+[n100/'eigenvectors'/f'qX_p2n100_evec_t{t}.h5' for t in (96,100)]:
            st=p.stat()
            if st.st_size<=0:raise ValueError('Empty regression input')
            inputs[str(p)]=dict(bytes=st.st_size,mtime_ns=st.st_mtime_ns)
    folder=campaign(args.name);folder.mkdir(parents=True,exist_ok=False)
    data=dict(name=args.name,stage=args.stage,cases=plan,runs=runs,inputs=inputs,
              executable_sha256=digest(EXE),build_manifest_sha256=manifest,
              package_sha256=digest(CASE/'package.sha256'),n100_run=str(n100),
              settings=dict(Nv=100,midpoint=98,alpha=[1,1,0.5],hyp_steps=1,grid=[1,4,4,4],ranks=64))
    write_json(folder/'PREPARED.json',data)
    (folder/'prepared.sha256').write_text(digest(folder/'PREPARED.json')+'\n')
    print('Prepared',len(plan),'cases:',args.name,args.stage)
    print('Test pilot first; full measurements use the micro partition in arrays of at most 39 cases.')


def prepared(name):
    verify_sources();folder=campaign(name);path=folder/'PREPARED.json'
    if digest(path)!=(folder/'prepared.sha256').read_text().strip():raise ValueError('Preparation changed')
    data=json.loads(path.read_text())
    if digest(CASE/'package.sha256')!=data['package_sha256'] or digest(EXE)!=data['executable_sha256'] or build_manifest()!=data['build_manifest_sha256']:
        raise ValueError('Measurement implementation changed')
    for name,expected in data['inputs'].items():
        st=Path(name).stat()
        if (st.st_size,st.st_mtime_ns)!=(expected['bytes'],expected['mtime_ns']):raise ValueError('Prepared measurement input changed: '+name)
    return folder,data


def read_case(folder,entry,data):
    check=checker();cfg,tau=entry['configuration'],entry['tau'];scans={};raw=[];windows=[]
    refs=reference_rows()
    for p in entry['piters']:
        points,records,times,window=check.read_measurement(folder/f'p{p}.out',cfg,tau,p)
        scans[p]=records;windows.append(window)
        raw.extend(check.raw_rows(cfg,tau,p,records))
        if cfg<=230:
            for key,row in records.items():
                ref=refs[(cfg,tau,40,*key)]
                if Decimal(row['S_RAW'])!=Decimal(ref['S_RAW']):raise ValueError('Original clover differs from archived cohort')
                if p in (4,40):
                    for field in check.FIELDS:check.avg.near(float(row[field]),float(refs[(cfg,tau,p,*key)][field]),'archived p4/p40 reproduction')
    for key in next(iter(scans.values())):
        if len({Decimal(scans[p][key]['S_RAW']) for p in scans})!=1:raise ValueError('Original clover changed across prescriptions')
    if data['stage']=='pilot':
        control={p:check.read_measurement(folder/f'n100_p{p}.out',100,4,p)[1] for p in (4,40)}
        check.regression_check(control)
    check.write_csv(folder/'profiles_raw.csv',check.RAW_COLUMNS,raw)
    check.write_csv(folder/'su3_windows.csv',check.SU3_COLUMNS,windows)
    return dict(success=True,configuration=cfg,tau=tau,piters=entry['piters'],raw_records=len(raw),window_records=len(windows),
                scope='Validated geometry, profile reductions, unchanged original clover and measured SU(3) diagnostics; no low-piter physics acceptance gate.')


def run(args):
    folder,data=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:
        raise ValueError('Run inside the prepared 8-node, 64-rank Slurm allocation')
    index=int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))
    if not 0<=index<len(data['cases']):raise ValueError('Array index outside campaign')
    entry=data['cases'][index];cfg,tau=entry['configuration'],entry['tau']
    out=folder/f'case_{index:03d}';out.mkdir(exist_ok=False)
    write_json(out/'CASE.json',dict(entry,prepared_sha256=digest(folder/'PREPARED.json'),job_id=os.environ['SLURM_JOB_ID']))
    chosen=next(r['run'] for r in data['runs'] if r['configuration']==cfg)
    tasks=[(cfg,tau,p,chosen,f'p{p}') for p in entry['piters']]
    if data['stage']=='pilot':tasks += [(100,4,p,data['n100_run'],f'n100_p{p}') for p in (4,40)]
    launches=[]
    try:
        for current,t,p,run_path,label in tasks:
            env=os.environ.copy();env.update(QX_CFG=str(current),QX_ROOT=str(GAUGES),EVROOT=str(Path(run_path)/'eigenvectors'),
                T0=str(98-t//2),TAU_LENGTH=str(t),HYP_MODE='2',HYP_PITER=str(p),PRINT_ALL='1',GX='0',GY='0',GZ='0')
            start=time.monotonic();command=['mpiexec','-n','64',str(EXE)]
            with (out/(label+'.out')).open('x') as stdout,(out/(label+'.err')).open('x') as stderr:
                code=subprocess.run(command,env=env,stdout=stdout,stderr=stderr,check=False).returncode
            launches.append(dict(configuration=current,tau=t,piter=p,label=label,exit_code=code,seconds=time.monotonic()-start))
            write_json(out/f'launch_{len(launches):02d}.json',launches[-1])
            if code:raise ValueError('MPI measurement failed: '+label)
        write_json(out/'LAUNCHES.json',launches)
        result=read_case(out,entry,data);write_json(out/'SUMMARY.json',result)
    except Exception as exc:
        write_json(out/'FAILED.json',dict(success=False,error=str(exc),launches=launches));raise
    print('PASS',args.name,index)


def pilot_gate(name):
    folder,data=prepared(name)
    if data['stage']!='pilot':raise ValueError('Require a completed pilot campaign')
    out=folder/'case_000'
    if (out/'FAILED.json').exists():raise ValueError('Pilot has a failure record')
    summary=json.loads((out/'SUMMARY.json').read_text())
    if summary.get('success') is not True:raise ValueError('Pilot failed')
    launches=json.loads((out/'LAUNCHES.json').read_text())
    wanted=[(data['cases'][0]['configuration'],4,p) for p in PITERS]+[(100,4,p) for p in (4,40)]
    if [(r['configuration'],r['tau'],r['piter']) for r in launches]!=wanted or any(r['exit_code'] for r in launches):raise ValueError('Pilot launch inventory failed')
    read_case(out,data['cases'][0],data)


def submit(args):
    folder,data=prepared(args.name)
    if data['stage']=='pilot':
        if len(data['cases'])!=1:raise ValueError('Pilot submission must contain exactly one configuration')
        command=['sbatch','--parsable',str(CASE/'submit_test.sh'),args.name]
    else:
        pilot_gate(args.pilot)
        stop=min(args.start+args.count,len(data['cases']))
        if not 0<=args.start<stop or stop-args.start>39 or not 1<=args.concurrent<=4:raise ValueError('Use at most 39 submitted cases and 1..4 concurrent tasks')
        if any((folder/f'case_{i:03d}').exists() for i in range(args.start,stop)):raise ValueError('A selected case already has output; inspect before resubmitting')
        command=['sbatch','--parsable',f'--array={args.start}-{stop-1}%{args.concurrent}',str(CASE/'submit_micro.sh'),args.name]
    # Record intent before submitting; a timeout is not proof of rejection.
    record=folder/f'submission_{args.start:03d}.json'
    write_json(record,dict(state='submitting',command=command))
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                          universal_newlines=True,check=False)
    data=dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
    record.write_text(json.dumps(data,indent=2)+'\n')
    if result.returncode:raise ValueError('Slurm submission rejected: '+result.stderr.strip())
    if not re.fullmatch(r'\d+(?:;\S+)?\n?',result.stdout):raise ValueError('Ambiguous Slurm response; inspect submission record before retrying')
    print('SUBMITTED',result.stdout.strip())


def collect(args):
    folder,data=prepared(args.name);check=checker();raw=[];windows=[];summaries=[]
    for entry in data['cases']:
        out=folder/f"case_{entry['index']:03d}"
        if (out/'FAILED.json').exists():raise ValueError('Case has a failure record')
        identity=json.loads((out/'CASE.json').read_text())
        if any(identity[k]!=v for k,v in entry.items()) or identity['prepared_sha256']!=digest(folder/'PREPARED.json'):raise ValueError('Wrong case identity')
        launches=json.loads((out/'LAUNCHES.json').read_text())
        wanted=[(entry['configuration'],entry['tau'],p,f'p{p}') for p in entry['piters']]
        if data['stage']=='pilot':wanted += [(100,4,p,f'n100_p{p}') for p in (4,40)]
        if [(r['configuration'],r['tau'],r['piter'],r['label']) for r in launches]!=wanted or any(r['exit_code'] for r in launches):raise ValueError('MPI launch inventory failed')
        summaries.append(read_case(out,entry,data))
        raw.extend(csv_rows(out/'profiles_raw.csv'));windows.extend(csv_rows(out/'su3_windows.csv'))
    target=folder/'collection';target.mkdir(exist_ok=False)
    check.write_csv(target/'profiles_raw.csv',check.RAW_COLUMNS,raw)
    check.write_csv(target/'su3_windows.csv',check.SU3_COLUMNS,windows)
    write_json(target/'SUMMARY.json',dict(success=True,stage=data['stage'],cases=summaries,raw_records=len(raw),window_records=len(windows)))
    archive=ROOT/f'qx_projection_spectrum_{args.name}_review.zip'
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as bundle:
        files=[p for p in folder.rglob('*') if p.is_file()]
        files += [p for p in CASE.iterdir() if p.is_file() and p.suffix not in ('.out','.err')]
        for p in sorted(set(files)):
            if p.suffix in ('.h5','.a','.o','.zip') or '__pycache__' in p.parts or p.stat().st_size>8*1024**2:raise ValueError('Unexpected review file')
            bundle.write(p,p.relative_to(ROOT).as_posix())
    print('REVIEW_ARCHIVE',archive)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action')
    p=sub.add_parser('prepare');p.add_argument('--name',required=True);p.add_argument('--stage',choices=['pilot','intermediate13','extended'],required=True)
    p.add_argument('--configs',default=','.join(map(str,range(110,231,10))))
    p=sub.add_parser('run');p.add_argument('--name',required=True)
    p=sub.add_parser('submit');p.add_argument('--name',required=True);p.add_argument('--pilot',default='pilot');p.add_argument('--start',type=int,default=0)
    p.add_argument('--count',type=int,default=39);p.add_argument('--concurrent',type=int,default=4)
    p=sub.add_parser('collect');p.add_argument('--name',required=True)
    args=parser.parse_args()
    if args.action is None:parser.error('Select prepare, run, submit or collect')
    globals()[args.action](args)


if __name__=='__main__':main()
