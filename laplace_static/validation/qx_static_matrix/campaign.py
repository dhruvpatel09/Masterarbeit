#!/usr/bin/env python3
"""Prepare, run and replay the three-configuration static-matrix pilot."""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import time
import zipfile
from pathlib import Path
import check_outputs as check
from matrix_plan import PILOT,CUTOFFS,RADII,TAUS,AUDIT_TIMES

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BUILD=HERE/'build'
GAUGES=Path('/hppfs/scratch/00/di98qes3/qX')
RUNS=ROOT/'laplace_static/validation/qx_static_extension/inputs/RUNS.json'
REFERENCE=ROOT/'laplace_static/analysis/qx_static_extension_71configs/inputs/job5958012/campaigns/late_modes_full71/correlators.csv'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def write(path,value,exclusive=False):
    path=Path(path)
    if exclusive and path.exists():raise ValueError('Receipt already exists: '+str(path))
    tmp=path.with_name(path.name+'.pending')
    with tmp.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def csv_write(path,values):
    path=Path(path);tmp=path.with_name(path.name+'.pending')
    with tmp.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(values[0]),lineterminator='\n');writer.writeheader();writer.writerows(values);f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def verify_sources():
    source=json.loads((HERE/'INPUTS.json').read_text())
    for name,expected in source['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Pinned source/reference changed: '+name)


def record_build():
    verify_sources()
    paths=[BUILD/'scan',BUILD/'audit_selected',ROOT/'qcd/lib/libqcd.a',
        ROOT/'solver/build_supermuc_qx_p1x4x4x4/libsolver.a',ROOT/'solver/build_supermuc_qx_p1x4x4x4/include/global.h',
        ROOT/'dist/disteigvecshdf5/c/src/DistEigvecsHdf5Reader.c',BUILD/'audit_selected.c']
    for folder in ('qcd/include','solver/include','dist/disteigvecshdf5/c/include'):
        paths+=sorted((ROOT/folder).glob('*.h'))
    paths+=sorted(BUILD.glob('*.h'))+sorted(BUILD.glob('*.txt'))
    write(BUILD/'MANIFEST.json',dict(inputs_sha256=digest(HERE/'INPUTS.json'),files={str(p.relative_to(ROOT)):digest(p) for p in paths}),exclusive=True)
    print('PASS isolated matrix build manifest')


def build_manifest():
    manifest=json.loads((BUILD/'MANIFEST.json').read_text())
    if manifest['inputs_sha256']!=digest(HERE/'INPUTS.json'):raise ValueError('Build source pins changed')
    for name,expected in manifest['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Compiled input changed: '+name)
    return digest(BUILD/'MANIFEST.json')


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def fingerprint(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size<=0:raise ValueError('Missing or empty physical input: '+str(path))
    stat=path.stat()
    return dict(path=str(path),resolved=str(path.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns)


def input_record(cfg):
    entries=json.loads(RUNS.read_text());entries=entries['runs'] if isinstance(entries,dict) else entries
    entry=next(e for e in entries if e['configuration']==cfg)
    if cfg==640:
        ev=Path(entry['eigenvectors']);values=Path(entry['eigenvalues'])
        if digest(values)!=entry['eigenvalue_sha256']:raise ValueError('Recovered n640 receipt changed')
    else:
        run=Path(entry['run'])
        if digest(run/'COMPLETED.json')!=entry['completion_marker_sha256']:raise ValueError('Saved-mode admission marker changed')
        ev=run/'eigenvectors';candidates=sorted((run/'eigenvalues').glob('*.h5'))
        if len(candidates)!=1:raise ValueError('Expected one complete eigenvalue table')
        values=candidates[0]
    return dict(configuration=cfg,eigenvectors=str(ev),eigenvalues=str(values),eigenvalue_sha256=digest(values),
        gauge=fingerprint(GAUGES/('qX_p2n%d'%cfg)),vectors=[fingerprint(ev/('qX_p2n%d_evec_t%d.h5'%(cfg,t))) for t in AUDIT_TIMES])


def validate_inputs(entry):
    cfg=entry['configuration']
    expected=['qX_p2n%d_evec_t%d.h5'%(cfg,t) for t in AUDIT_TIMES]
    if [Path(r['path']).name for r in entry['vectors']]!=expected:raise ValueError('Wrong saved-mode time inventory')
    if Path(entry['gauge']['path']).name!='qX_p2n%d'%cfg:raise ValueError('Wrong gauge configuration')
    if digest(entry['eigenvalues'])!=entry['eigenvalue_sha256']:raise ValueError('Eigenvalue file changed')
    for record in [entry['gauge']]+entry['vectors']:
        if fingerprint(record['path'])!=record:raise ValueError('Prepared physical input changed: '+record['path'])


def prepared(name,physical=False):
    verify_sources();target=folder(name);data=json.loads((target/'PREPARED.json').read_text())
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared plan checksum mismatch')
    if data['build_manifest_sha256']!=build_manifest():raise ValueError('Build differs from prepared plan')
    if [e['configuration'] for e in data['cases']]!=list(PILOT) or [e['index'] for e in data['cases']]!=list(range(3)):raise ValueError('Pilot inventory changed')
    if any(e['piters']!=([10,40] if e['configuration']==650 else [10]) for e in data['cases']):raise ValueError('Pilot projection controls changed')
    if data['geometry']!=geometry() or data['stage']!='technical_pilot' or data['target_production_configurations']!=list(range(100,801,10)):raise ValueError('Pilot geometry changed')
    if any(e['input']['configuration']!=e['configuration'] for e in data['cases']):raise ValueError('Physical input belongs to another configuration')
    if physical:
        for entry in data['cases']:validate_inputs(entry['input'])
    return target,data


def geometry():
    return dict(t0=90,taus=list(TAUS),R=list(RADII),Nv_src=list(CUTOFFS),Nv_sink=list(CUTOFFS),saved_Nv=100)


def prepare(args):
    verify_sources();manifest=build_manifest();check.load_references(REFERENCE,PILOT)
    cases=[dict(index=i,configuration=cfg,piters=[10,40] if cfg==650 else [10],input=input_record(cfg)) for i,cfg in enumerate(PILOT)]
    target=folder(args.name);target.mkdir(parents=True,exist_ok=False)
    data=dict(stage='technical_pilot',target_production_configurations=list(range(100,801,10)),cases=cases,
        build_manifest_sha256=manifest,geometry=geometry(),
        interpretation='All 16 source/sink combinations; no per-configuration Hermitian symmetrization at nonzero time; statistical unit is one configuration.')
    write(target/'PREPARED.json',data,exclusive=True);(target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared matrix pilot: n170, n640, n650; three audits, three p10 scans, one n650 p40 control; three array tasks')


def complete_case(target,entry):
    marker=target/('case_%03d'%entry['index'])/'LATEST.json'
    if not marker.is_file():return None
    record=json.loads(marker.read_text())
    if not re.fullmatch(r'attempt_\d{4}',record['attempt']):raise ValueError('Invalid attempt path')
    out=marker.parent/record['attempt']
    if out.is_symlink():raise ValueError('Attempt directory is a symlink')
    summary=json.loads((out/'SUMMARY.json').read_text())
    if digest(out/'SUMMARY.json')!=record['summary_sha256'] or summary.get('success') is not True:raise ValueError('Changed validation marker')
    if digest(out/'LAUNCHES.json')!=summary['launches_sha256'] or digest(out/'CASE.json')!=summary['case_sha256']:raise ValueError('Changed case provenance')
    for name in ('matrix_correlators.csv','rank_blocks.csv'):
        if digest(out/name)!=summary['outputs_sha256'][name]:raise ValueError('Changed validated numerical export')
    if json.loads((out/'CASE.json').read_text())!=dict(entry,prepared_sha256=digest(target/'PREPARED.json')):raise ValueError('Case belongs to another plan')
    return out


def replay_case(out,entry,references):
    cfg=entry['configuration'];launches=json.loads((out/'LAUNCHES.json').read_text())
    if [(r['kind'],r['piter'],r['exit_code'],r['stem']) for r in launches]!=[('audit',0,0,'audit')]+[('scan',p,0,'p%d'%p) for p in entry['piters']]:raise ValueError('Incomplete MPI launch inventory')
    for r in launches:
        if r['configuration']!=cfg or not math.isfinite(r['seconds']) or r['seconds']<=0:raise ValueError('Invalid launch metadata or timing')
        for suffix in ('out','err'):
            if digest(out/(r['stem']+'.'+suffix))!=r[suffix+'_sha256']:raise ValueError('Archived output changed')
    quality=check.audit(out/'audit.out',cfg);values=[];blocks=[];scans={};summaries=[]
    for p in entry['piters']:
        rows,parts,summary=check.scan(out/('p%d.out'%p),cfg,p,references)
        values+=rows;blocks+=parts;scans[p]=rows;summaries.append(summary)
    difference=check.projection_pair(scans[10],scans[40]) if 40 in scans else None
    return values,blocks,dict(configuration=cfg,audit=quality,scans=summaries,p10_p40_maximum_component_difference=difference)


def live_previous(target):
    receipts=sorted(target.glob('SUBMISSION_*.json'))
    if not receipts:return
    previous=[]
    for path in receipts:
        record=json.loads(path.read_text());match=re.fullmatch(r'(\d+)(?:;\S+)?\n?',record.get('stdout',''))
        if record.get('returncode')!=0 or not match:raise ValueError('Ambiguous prior submission; inspect its receipt')
        previous.append(match.group(1))
    result=subprocess.run(['squeue','-h','-u',os.environ['USER'],'-o','%A %F'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if result.returncode or set(previous)&set(re.findall(r'\d+',result.stdout)):raise ValueError('Prior submission may still be active; no duplicate submission made')


def submit(args):
    target,data=prepared(args.name,physical=True)
    if not 1<=args.concurrent<=3:raise ValueError('Pilot concurrency must be 1..3')
    live_previous(target);indices=[e['index'] for e in data['cases'] if complete_case(target,e) is None]
    if not indices:print('Pilot tasks already complete; collect the results');return
    command=['sbatch','--parsable','--array='+','.join(map(str,indices))+'%'+str(args.concurrent),str(HERE/'submit_micro.sh'),args.name]
    receipt=target/('SUBMISSION_%03d.json'%(len(list(target.glob('SUBMISSION_*.json')))+1))
    write(receipt,dict(state='submitting',command=command),exclusive=True)
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    write(receipt,dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?',result.stdout):raise ValueError('Submission failed or returned an ambiguous job ID')
    print('SUBMITTED',result.stdout.strip(),len(indices),'matrix pilot tasks')


def launch(out,kind,piter,entry,launches):
    cfg=entry['configuration'];source=entry['input'];stem='audit' if kind=='audit' else 'p%d'%piter
    command=['mpiexec','-n','64',str(BUILD/('audit_selected' if kind=='audit' else 'scan'))]
    if kind=='audit':command += [str(cfg),source['gauge']['path'],source['eigenvectors'],source['eigenvalues']]
    env=os.environ.copy();env.update(QX_CFG=str(cfg),HYP_PITER=str(piter),EVROOT=source['eigenvectors'],QX_ROOT=str(Path(source['gauge']['path']).parent))
    start=time.monotonic()
    with (out/(stem+'.out')).open('x') as stdout,(out/(stem+'.err')).open('x') as stderr:
        result=subprocess.run(command,env=env,stdout=stdout,stderr=stderr)
    record=dict(kind=kind,piter=piter,configuration=cfg,stem=stem,seconds=time.monotonic()-start,exit_code=result.returncode,
        command=command,job_id=os.environ['SLURM_JOB_ID'],array_task_id=os.environ['SLURM_ARRAY_TASK_ID'],out_sha256=digest(out/(stem+'.out')),err_sha256=digest(out/(stem+'.err')))
    launches.append(record);write(out/'LAUNCHES.json',launches)
    if result.returncode:raise ValueError('MPI launch failed: '+stem)


def run(args):
    target,data=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:raise ValueError('Eight-node, 64-rank Slurm allocation required')
    index=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0<=index<len(data['cases']):raise ValueError('Unknown array task')
    entry=data['cases'][index];validate_inputs(entry['input']);references=check.load_references(REFERENCE,[entry['configuration']])
    old=complete_case(target,entry)
    if old:replay_case(old,entry,references);print('REUSED matrix n%d'%entry['configuration']);return
    case=target/('case_%03d'%index);case.mkdir(exist_ok=True)
    numbers=[int(p.name.split('_')[1]) for p in case.glob('attempt_*') if re.fullmatch(r'attempt_\d{4}',p.name)]
    out=case/('attempt_%04d'%(max(numbers,default=0)+1));out.mkdir();launches=[]
    write(out/'CASE.json',dict(entry,prepared_sha256=digest(target/'PREPARED.json')),exclusive=True)
    try:
        launch(out,'audit',0,entry,launches);check.audit(out/'audit.out',entry['configuration'])
        for piter in entry['piters']:launch(out,'scan',piter,entry,launches)
        values,blocks,summary=replay_case(out,entry,references)
        csv_write(out/'matrix_correlators.csv',values);csv_write(out/'rank_blocks.csv',blocks)
        validate_inputs(entry['input'])
        write(out/'SUMMARY.json',dict(summary,success=True,launches_sha256=digest(out/'LAUNCHES.json'),case_sha256=digest(out/'CASE.json'),
            outputs_sha256={name:digest(out/name) for name in ('matrix_correlators.csv','rank_blocks.csv')}),exclusive=True)
        write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=digest(out/'SUMMARY.json')))
    except Exception as exc:
        write(out/'FAILED.json',dict(success=False,error=str(exc),launches=launches));raise
    print('PASS matrix pilot n%d'%entry['configuration'])


def collect(args):
    target,data=prepared(args.name);references=check.load_references(REFERENCE,PILOT);values=[];blocks=[];cases=[];missing=[];timing=[]
    for entry in data['cases']:
        out=complete_case(target,entry)
        if out is None:missing.append(entry['configuration']);continue
        rows,parts,summary=replay_case(out,entry,references);values+=rows;blocks+=parts;cases.append(summary)
        timing.append(dict(configuration=entry['configuration'],launches=json.loads((out/'LAUNCHES.json').read_text())))
    success=not missing
    if values:csv_write(target/'matrix_correlators.csv',values);csv_write(target/'rank_blocks.csv',blocks)
    summary=dict(success=success,technical_pilot_configurations=list(PILOT),validated_configurations=[r['configuration'] for r in cases],missing_configurations=missing,
        production_target_configurations=list(range(100,801,10)),matrix_rows=len(values),rank_block_rows=len(blocks),cases=cases,timing=timing,
        interpretation='Technical validation on three configurations. No GEVP spectrum or ensemble-level physics claim; full production target remains all 71.')
    write(target/'SUMMARY.json',summary)
    archive=ROOT/('qx_static_matrix_'+args.name+'_review.zip');tmp=archive.with_name(archive.name+'.pending')
    files=[p for p in target.rglob('*') if p.is_file()]
    files += [p for p in HERE.rglob('*') if p.is_file() and 'build' not in p.relative_to(HERE).parts and 'campaigns' not in p.relative_to(HERE).parts and '__pycache__' not in p.parts]
    files += [p for p in BUILD.iterdir() if p.is_file() and (p.suffix in ('.h','.txt') or p.name in ('MANIFEST.json','audit_selected.c'))]
    if tmp.exists():raise ValueError('Unfinished prior review archive; inspect '+str(tmp))
    with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(set(files)):
            if p.is_symlink() or p.name.endswith('.pending'):raise ValueError('Unfinished or linked campaign output '+str(p))
            z.write(p,str(p.relative_to(ROOT)))
    with tmp.open('rb') as f:os.fsync(f.fileno())
    tmp.replace(archive)
    print(('PASS' if success else 'INCOMPLETE')+' matrix pilot; validated '+str(len(cases))+'/3; missing '+str(missing))
    print('REVIEW_ARCHIVE',archive)
    if not success:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action');sub.required=True
    for name in ('check-inputs','record-build'):sub.add_parser(name)
    for name in ('prepare','submit','run','collect'):
        p=sub.add_parser(name);p.add_argument('--name',required=True)
        if name=='submit':p.add_argument('--concurrent',type=int,default=3)
    args=parser.parse_args()
    if args.action=='check-inputs':verify_sources();print('PASS matrix source and reference pins')
    elif args.action=='record-build':record_build()
    else:globals()[args.action](args)


if __name__=='__main__':main()
