#!/usr/bin/env python3
"""Fixed-source qX correlators through R/a=32, targeting all 71 configurations."""
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

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=list(range(100,801,10))
PILOT=[170,640,650]
BUILD=HERE/'build'
GAUGES=Path('/hppfs/scratch/00/di98qes3/qX')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write(path,value,exclusive=True):
    with Path(path).open('x' if exclusive else 'w') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def csv_write(path,values):
    with Path(path).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(values)


def verify_sources():
    for name,expected in json.loads((HERE/'INPUTS.json').read_text())['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Pinned source changed: '+name)


def record_build():
    files=[BUILD/'scan',BUILD/'audit_selected',ROOT/'qcd/lib/libqcd.a',
        ROOT/'solver/build_supermuc_qx_p1x4x4x4/libsolver.a',
        ROOT/'solver/build_supermuc_qx_p1x4x4x4/include/global.h']
    files+=sorted((ROOT/'qcd/include').glob('*.h'))+sorted((ROOT/'solver/include').glob('*.h'))
    files+=sorted((ROOT/'dist/disteigvecshdf5/c/include').glob('*.h'))
    files+=[ROOT/'dist/disteigvecshdf5/c/src/DistEigvecsHdf5Reader.c']
    files+=sorted(BUILD.glob('*.h'))+[BUILD/'audit_selected.c']
    write(BUILD/'MANIFEST.json',dict(inputs_sha256=digest(HERE/'INPUTS.json'),
        files={str(path.relative_to(ROOT)):digest(path) for path in files}))
    print('PASS isolated static-correlator build manifest')


def build_manifest():
    manifest=json.loads((BUILD/'MANIFEST.json').read_text())
    if manifest['inputs_sha256']!=digest(HERE/'INPUTS.json'):raise ValueError('Build uses different source pins')
    for name,expected in manifest['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Compiled input changed: '+name)
    return digest(BUILD/'MANIFEST.json')


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def fingerprint(path):
    path=Path(path)
    if not path.is_file():raise ValueError('Missing input: '+str(path))
    stat=path.stat()
    if stat.st_size<=0:raise ValueError('Empty input: '+str(path))
    return dict(path=str(path),resolved=str(path.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns)


def input_record(cfg):
    entries=json.loads((HERE/'inputs/RUNS.json').read_text())
    # RUNS is a list of admitted runs, including the recovered n640 path.
    if isinstance(entries,dict):entries=entries['runs']
    entry=next(r for r in entries if r['configuration']==cfg)
    if cfg==640:
        ev=Path(entry['eigenvectors']);values=Path(entry['eigenvalues'])
        if digest(values)!=entry['eigenvalue_sha256']:raise ValueError('n640 eigenvalue receipt changed')
    else:
        run=Path(entry['run'])
        if digest(run/'COMPLETED.json')!=entry['completion_marker_sha256']:raise ValueError('Run admission marker changed')
        ev=run/'eigenvectors';candidates=sorted((run/'eigenvalues').glob('*.h5'))
        if len(candidates)!=1:raise ValueError('Expected one complete eigenvalue table for n%d'%cfg)
        values=candidates[0]
    return dict(configuration=cfg,eigenvectors=str(ev),eigenvalues=str(values),
        eigenvalue_sha256=digest(values),gauge=fingerprint(GAUGES/('qX_p2n%d'%cfg)),
        vectors=[fingerprint(ev/('qX_p2n%d_evec_t%d.h5'%(cfg,t))) for t in range(90,107)])


def validate_inputs(entry):
    if digest(entry['eigenvalues'])!=entry['eigenvalue_sha256']:raise ValueError('Eigenvalue file changed')
    for old in [entry['gauge']]+entry['vectors']:
        if fingerprint(old['path'])!=old:raise ValueError('Prepared physical input changed: '+old['path'])


def prepared(name,physical=False):
    verify_sources();target=folder(name)
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared plan checksum mismatch')
    data=json.loads((target/'PREPARED.json').read_text())
    if data['build_manifest_sha256']!=build_manifest():raise ValueError('Build differs from prepared plan')
    if physical:
        for entry in data['cases']:validate_inputs(entry['input'])
    return target,data


def references():
    with (HERE/'inputs/reference_correlators.csv').open(newline='') as stream:
        return {tuple(int(r[k]) for k in ('configuration','piter','tau','R','axis')):r for r in csv.DictReader(stream)}


def complete_case(target,entry):
    marker=target/('case_%03d'%entry['index'])/'LATEST.json'
    if not marker.is_file():return None
    record=json.loads(marker.read_text())
    if not re.fullmatch(r'attempt_\d{4}',record['attempt']):raise ValueError('Invalid attempt path')
    out=marker.parent/record['attempt']
    summary=json.loads((out/'SUMMARY.json').read_text())
    if digest(out/'SUMMARY.json')!=record['summary_sha256'] or not summary['success']:
        raise ValueError('Changed case validation')
    if digest(out/'LAUNCHES.json')!=summary['launches_sha256'] or digest(out/'CASE.json')!=summary['case_sha256']:
        raise ValueError('Case provenance changed')
    case=json.loads((out/'CASE.json').read_text())
    if case!=dict(entry,prepared_sha256=digest(target/'PREPARED.json')):raise ValueError('Case belongs to another prepared plan')
    return out


def replay_case(out,entry,refs):
    cfg=entry['configuration'];launches=json.loads((out/'LAUNCHES.json').read_text())
    if [(r['kind'],r['piter'],r['exit_code']) for r in launches]!=[('audit',0,0)]+[('scan',p,0) for p in entry['piters']]:
        raise ValueError('Incomplete MPI launch inventory')
    for launch in launches:
        if not math.isfinite(launch['seconds']) or launch['seconds']<=0:raise ValueError('Invalid elapsed time')
        for suffix in ('out','err'):
            if digest(out/(launch['stem']+'.'+suffix))!=launch[suffix+'_sha256']:raise ValueError('Archived log changed')
    quality=check.audit(out/'audit.out',cfg);values=[];blocks=[];scans={};summaries=[]
    for p in entry['piters']:
        rows,parts,summary=check.scan(out/('p%d.out'%p),cfg,p,refs)
        values+=rows;blocks+=parts;scans[p]=rows;summaries.append(summary)
    difference=check.projection_pair(scans[10],scans[40]) if 40 in scans else None
    return values,blocks,dict(configuration=cfg,quality=quality,scans=summaries,p10_p40_maximum_absolute_difference=difference)


def pilot_cases(name):
    target,data=prepared(name)
    if data['stage']!='pilot' or [e['configuration'] for e in data['cases']]!=PILOT:raise ValueError('Expected the matched extended pilot')
    refs=references();times=[]
    for entry in data['cases']:
        out=complete_case(target,entry)
        if out is None:raise ValueError('Pilot is incomplete; collect and inspect it before the full run')
        replay_case(out,entry,refs)
        launches=json.loads((out/'LAUNCHES.json').read_text())
        times.append(sum(r['seconds'] for r in launches if r['piter'] in (0,10)))
    return target,data,times


def case_plan(stage,reuse=(),per_task=1):
    configs=PILOT if stage=='pilot' else [c for c in TARGET if c not in reuse]
    cases=[dict(index=i,configuration=c,piters=[10,40] if stage=='pilot' else [10]) for i,c in enumerate(configs)]
    return cases,[list(range(i,min(i+per_task,len(cases)))) for i in range(0,len(cases),per_task)]


def prepare(args):
    verify_sources();manifest=build_manifest();reuse=[];pilot_record=None;size=1
    if args.stage=='full':
        pilot_target,pilot_data,times=pilot_cases(args.pilot)
        # At least 50% timing margin inside 50 minutes of the one-hour allocation.
        size=2 if 1.5*max(times)*2<=3000 else 1
        if 1.5*max(times)>3000:raise ValueError('Pilot timing leaves insufficient margin; review batching before production')
        reuse=PILOT[:]
        pilot_record=dict(name=args.pilot,prepared_sha256=digest(pilot_target/'PREPARED.json'),
            cases=[dict(configuration=e['configuration'],attempt=str(complete_case(pilot_target,e))) for e in pilot_data['cases']],
            maximum_audit_plus_p10_seconds=max(times))
    cases,slots=case_plan(args.stage,reuse,size)
    for entry in cases:entry['input']=input_record(entry['configuration'])
    target=folder(args.name);target.mkdir(parents=True,exist_ok=False)
    data=dict(stage=args.stage,target_configurations=TARGET,new_configurations=[e['configuration'] for e in cases],
        reused_configurations=reuse,pilot=pilot_record,cases=cases,slots=slots,cases_per_task=size,
        build_manifest_sha256=manifest,geometry=dict(t0=90,taus=list(range(17)),R=list(range(1,33)),Nv=100),
        statistical_unit='One configuration after source and spatial-axis averaging; block outputs are correlated regional diagnostics.')
    write(target/'PREPARED.json',data);(target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared %s: %d new configurations; %d reused; %d array tasks; target all 71.'%(args.stage,len(cases),len(reuse),len(slots)))


def live_previous(target):
    records=sorted(target.glob('SUBMISSION_*.json'))
    if not records:return
    previous=[]
    for path in records:
        record=json.loads(path.read_text());match=re.fullmatch(r'(\d+)(?:;\S+)?\n?',record.get('stdout',''))
        if record.get('returncode')!=0 or not match:raise ValueError('Ambiguous prior submission; inspect its receipt')
        previous.append(match.group(1))
    result=subprocess.run(['squeue','-h','-u',os.environ['USER'],'-o','%A %F'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if result.returncode or set(previous)&set(re.findall(r'\d+',result.stdout)):
        raise ValueError('A prior submission may still be active; no duplicate submission made')


def submit(args):
    target,data=prepared(args.name,physical=True)
    if not 1<=args.concurrent<=4 or not 1<=args.count<=39:raise ValueError('Use concurrency 1..4 and at most 39 array tasks per submission')
    if not 0<=args.start<len(data['slots']):raise ValueError('Array start outside plan')
    live_previous(target)
    slots=[s for s in range(args.start,min(args.start+args.count,len(data['slots'])))
        if any(complete_case(target,data['cases'][i]) is None for i in data['slots'][s])]
    if not slots:print('Selected tasks are complete; collect their results');return
    command=['sbatch','--parsable','--array='+','.join(map(str,slots))+'%'+str(args.concurrent),str(HERE/'submit_micro.sh'),args.name]
    number=len(list(target.glob('SUBMISSION_*.json')))+1;receipt=target/('SUBMISSION_%03d.json'%number)
    write(receipt,dict(command=command,state='submitting'))
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    write(receipt,dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr),exclusive=False)
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?',result.stdout):raise ValueError('Submission failed or returned an ambiguous job ID')
    print('SUBMITTED',result.stdout.strip(),len(slots),'array tasks')


def launch(out,kind,piter,entry,launches):
    cfg=entry['configuration'];source=entry['input'];stem='audit' if kind=='audit' else 'p%d'%piter
    if kind=='audit':command=['mpiexec','-n','64',str(BUILD/'audit_selected'),str(cfg),source['gauge']['path'],source['eigenvectors'],source['eigenvalues']]
    else:command=['mpiexec','-n','64',str(BUILD/'scan')]
    env=os.environ.copy();env.update(QX_CFG=str(cfg),HYP_PITER=str(piter),EVROOT=source['eigenvectors'],QX_ROOT=str(Path(source['gauge']['path']).parent))
    start=time.monotonic()
    with (out/(stem+'.out')).open('x') as stdout,(out/(stem+'.err')).open('x') as stderr:
        result=subprocess.run(command,env=env,stdout=stdout,stderr=stderr)
    record=dict(kind=kind,piter=piter,configuration=cfg,stem=stem,seconds=time.monotonic()-start,exit_code=result.returncode,
        command=command,job_id=os.environ['SLURM_JOB_ID'],out_sha256=digest(out/(stem+'.out')),err_sha256=digest(out/(stem+'.err')))
    launches.append(record);write(out/'LAUNCHES.json',launches,exclusive=False)
    if result.returncode:raise ValueError('MPI launch failed: '+stem)


def run(args):
    target,data=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:
        raise ValueError('An eight-node, 64-rank array allocation is required')
    slot=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0<=slot<len(data['slots']):raise ValueError('Unknown array task')
    refs=references()
    for index in data['slots'][slot]:
        entry=data['cases'][index];validate_inputs(entry['input'])
        old=complete_case(target,entry)
        if old:replay_case(old,entry,refs);print('REUSED n%d'%entry['configuration']);continue
        case=target/('case_%03d'%index);case.mkdir(exist_ok=True)
        numbers=[int(p.name.split('_')[1]) for p in case.glob('attempt_*') if re.fullmatch(r'attempt_\d{4}',p.name)]
        out=case/('attempt_%04d'%(max(numbers,default=0)+1));out.mkdir();launches=[]
        write(out/'CASE.json',dict(entry,prepared_sha256=digest(target/'PREPARED.json')))
        try:
            launch(out,'audit',0,entry,launches);check.audit(out/'audit.out',entry['configuration'])
            for p in entry['piters']:launch(out,'scan',p,entry,launches)
            values,blocks,summary=replay_case(out,entry,refs)
            csv_write(out/'correlators.csv',values);csv_write(out/'spatial_blocks.csv',blocks)
            validate_inputs(entry['input'])
            write(out/'SUMMARY.json',dict(summary,success=True,launches_sha256=digest(out/'LAUNCHES.json'),case_sha256=digest(out/'CASE.json')))
            tmp=case/'LATEST.tmp';write(tmp,dict(attempt=out.name,summary_sha256=digest(out/'SUMMARY.json')),exclusive=False);tmp.replace(case/'LATEST.json')
        except Exception as exc:
            write(out/'FAILED.json',dict(success=False,error=str(exc),launches=launches));raise
        print('PASS fixed-source scan n%d'%entry['configuration'])


def collect(args):
    target,data=prepared(args.name);refs=references();values=[];blocks=[];summaries=[];missing=[];attempts=[]
    for entry in data['cases']:
        out=complete_case(target,entry)
        if out is None:missing.append(entry['configuration']);continue
        rows,parts,summary=replay_case(out,entry,refs);values+=rows;blocks+=parts;summaries.append(summary);attempts.append(out)
    if data['stage']=='full':
        pilot_target,pilot_data,_=pilot_cases(data['pilot']['name'])
        if digest(pilot_target/'PREPARED.json')!=data['pilot']['prepared_sha256']:raise ValueError('Reused pilot plan changed')
        for entry in pilot_data['cases']:
            out=complete_case(pilot_target,entry);rows,parts,summary=replay_case(out,entry,refs)
            values += [r for r in rows if int(r['piter'])==10]
            blocks += [r for r in parts if int(r['piter'])==10]
            summaries.append(summary);attempts.append(out)
    expected=PILOT if data['stage']=='pilot' else TARGET
    actual=sorted({int(r['cfg']) for r in values})
    success=not missing and actual==expected
    if values:
        csv_write(target/'correlators.csv',values);csv_write(target/'spatial_blocks.csv',blocks)
    summary=dict(success=success,stage=data['stage'],target_configurations=expected,validated_configurations=actual,
        missing_configurations=missing,correlator_rows=len(values),spatial_block_rows=len(blocks),cases=summaries,
        timing=[dict(configuration=json.loads((p/'CASE.json').read_text())['configuration'],launches=json.loads((p/'LAUNCHES.json').read_text())) for p in attempts],
        interpretation='Correlators, not fitted potentials. Ensemble ratio/log estimates, positivity, signal quality and plateau tests follow collection.')
    write(target/'SUMMARY.json',summary,exclusive=False)
    archive=ROOT/('qx_static_potential_'+args.name+'_review.zip');temporary=archive.with_suffix('.zip.tmp')
    files=[p for p in target.rglob('*') if p.is_file()]
    for out in attempts:
        files += [p for p in out.rglob('*') if p.is_file()]
    if data['stage']=='full':
        files += [p for p in pilot_target.glob('*') if p.is_file()]
    files += [p for p in HERE.rglob('*') if p.is_file() and 'campaigns' not in p.relative_to(HERE).parts and '__pycache__' not in p.parts and 'build' not in p.relative_to(HERE).parts]
    files += [p for p in BUILD.glob('*') if p.is_file() and p.suffix in ('.json','.txt','.h','.c')]
    with zipfile.ZipFile(str(temporary),'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted(set(files)):z.write(str(path),str(path.relative_to(ROOT)))
    temporary.replace(archive)
    print('PASS' if success else 'INCOMPLETE','static scan;',len(actual),'configurations; missing:',','.join(map(str,missing)) or 'none')
    print('REVIEW_ARCHIVE',archive)
    if not success:raise SystemExit(2)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action');sub.required=True
    sub.add_parser('check-inputs');sub.add_parser('record-build')
    p=sub.add_parser('prepare');p.add_argument('--name',required=True);p.add_argument('--stage',choices=('pilot','full'),default='pilot');p.add_argument('--pilot',default='pilot')
    p=sub.add_parser('submit');p.add_argument('--name',required=True);p.add_argument('--start',type=int,default=0);p.add_argument('--count',type=int,default=39);p.add_argument('--concurrent',type=int,default=3)
    for name in ('run','collect'):
        p=sub.add_parser(name);p.add_argument('--name',required=True)
    args=parser.parse_args()
    if args.action=='check-inputs':verify_sources();print('PASS static-correlator source pins')
    elif args.action=='record-build':verify_sources();record_build()
    else:globals()[args.action](args)


if __name__=='__main__':main()
