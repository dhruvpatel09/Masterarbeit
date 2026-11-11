#!/usr/bin/env python3
"""Build receipts, physical admission, three-case native pilot and strict replay."""
import argparse,json,math,os,re,subprocess,time,zipfile
from pathlib import Path
from parents import HERE,ROOT,DESIGN,inventory,design
from pilot_plan import PILOT,execution_scans,geometry
from receipts import digest,read,write,csv_write,csv_equal
import check_outputs as check
from fixed_projections import derive

BUILD=HERE/'build'
RUNS=ROOT/'laplace_static/validation/qx_static_extension/inputs/RUNS.json'
REFERENCE=DESIGN/'inputs/reference_t90_axis_average.csv'
SOURCES=DESIGN/'inputs/FIXED_SOURCES.json'
GAUGES=Path('/hppfs/scratch/00/di98qes3/qX')


def verify_sources():
    from verify_files import verify
    verify()


def record_build():
    verify_sources()
    paths=[HERE/'scan.c',ROOT/'qcd/lib/libqcd.a',ROOT/'solver/build_supermuc_qx_p1x4x4x4/libsolver.a',
        ROOT/'solver/build_supermuc_qx_p1x4x4x4/include/global.h',ROOT/'dist/disteigvecshdf5/c/src/DistEigvecsHdf5Reader.c']
    paths += [BUILD/name for name in ('scan','audit66','audit97')]
    for directory in ('qcd/include','solver/include','dist/disteigvecshdf5/c/include'):
        paths+=sorted((ROOT/directory).glob('*.h'))
    paths+=[p for p in BUILD.iterdir() if p.is_file() and p.suffix in ('.c','.h','.txt')]
    write(BUILD/'MANIFEST.json',dict(package_manifest_sha256=digest(HERE/'FILE_MANIFEST.json'),
        files={str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths))}),exclusive=True)
    print('PASS isolated source-time native build manifest')


def build_manifest():
    manifest=read(BUILD/'MANIFEST.json')
    if manifest['package_manifest_sha256']!=digest(HERE/'FILE_MANIFEST.json'):raise ValueError('Native build differs from source package')
    for name,expected in manifest['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Compiled input changed: '+name)
    return digest(BUILD/'MANIFEST.json')


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def physical_cases():
    report=inventory.inspect_inputs(read(RUNS),GAUGES)
    if not report['success']:raise ValueError('Physical preflight blocked: '+json.dumps([{k:c[k] for k in ('configuration','errors')} for c in report['cases'] if not c['available']]))
    return {c['configuration']:c for c in report['cases']}


def validate_inputs(entry):
    cfg=entry['configuration'];data=entry['input']
    if data['configuration']!=cfg or data['required_times']!=design.audit_times(cfg) or [v['time'] for v in data['vectors']]!=design.audit_times(cfg):raise ValueError('Prepared physical time inventory changed')
    if any(Path(v['path']).name!='qX_p2n%d_evec_t%d.h5'%(cfg,v['time']) for v in data['vectors']):raise ValueError('Prepared eigenvector identity changed')
    if Path(data['gauge']['path']).name!='qX_p2n%d'%cfg:raise ValueError('Prepared gauge identity changed')
    directories={str(Path(v['path']).parent) for v in data['vectors']}
    if len(directories)!=1:raise ValueError('More than one physical eigenvector root')
    for record in [data['gauge'],data['eigenvalues']]+data['vectors']:
        current=inventory.fingerprint(record['path'])
        if current!={k:record[k] for k in current}:raise ValueError('Prepared physical input changed: '+record['path'])
    if digest(data['eigenvalues']['path'])!=data['eigenvalues']['sha256']:raise ValueError('Prepared eigenvalue table changed')
    admitted=data['admission']
    if admitted['kind']=='saved_mode_completion_marker':
        if digest(admitted['marker'])!=admitted['sha256']:raise ValueError('Prepared saved-mode marker changed')
    elif cfg!=640 or admitted['kind']!='recovered_n640_eigenvalue_checksum' or admitted['expected_sha256']!=data['eigenvalues']['sha256']:
        raise ValueError('Wrong recovered-mode admission receipt')


def prepared(name,physical=False):
    verify_sources();target=folder(name);data=read(target/'PREPARED.json')
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared plan checksum mismatch')
    if data['build_manifest_sha256']!=build_manifest() or data['design_plan_sha256']!=digest(DESIGN/'PLAN.json'):raise ValueError('Prepared build/design changed')
    if data['geometry']!=geometry() or data['stage']!='source_time_technical_pilot' or [e['configuration'] for e in data['cases']]!=list(PILOT):raise ValueError('Prepared native pilot inventory changed')
    for i,entry in enumerate(data['cases']):
        if entry['index']!=i or entry['scans']!=execution_scans(entry['configuration']):raise ValueError('Prepared scan order changed')
        if physical:validate_inputs(entry)
    return target,data


def prepare(args):
    verify_sources();manifest=build_manifest();check.load_references(REFERENCE);physical=physical_cases()
    cases=[dict(index=i,configuration=cfg,scans=execution_scans(cfg),input=physical[cfg]) for i,cfg in enumerate(PILOT)]
    for entry in cases:validate_inputs(entry)
    target=folder(args.name);target.mkdir(parents=True,exist_ok=False)
    write(target/'PREPARED.json',dict(stage='source_time_technical_pilot',cases=cases,geometry=geometry(),
        build_manifest_sha256=manifest,design_plan_sha256=digest(DESIGN/'PLAN.json'),
        source_vectors_sha256=digest(SOURCES),reference_sha256=digest(REFERENCE),
        interpretation='Three technical configurations; frozen 68-configuration training sample; no ensemble potential fits.'),exclusive=True)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared source-time pilot: three tasks; 66/66/97-slice audits; eight scans; eleven MPI launches.')
    print('n650 execution order: fresh t90 regression, shifted t58, shifted t122, t122 p40 control.')


def complete_case(target,entry):
    case=target/('case_%03d'%entry['index']);marker=case/'LATEST.json'
    if not marker.exists():return None
    record=read(marker)
    if marker.is_symlink() or not re.fullmatch(r'attempt_\d{4}',record['attempt']):raise ValueError('Invalid native attempt marker')
    out=case/record['attempt']
    if out.is_symlink():raise ValueError('Linked native attempt directory')
    summary=read(out/'SUMMARY.json')
    if digest(out/'SUMMARY.json')!=record['summary_sha256'] or summary.get('success') is not True:raise ValueError('Changed validated native marker')
    if digest(out/'CASE.json')!=summary['case_sha256'] or digest(out/'LAUNCHES.json')!=summary['launches_sha256']:raise ValueError('Changed native case provenance')
    if read(out/'CASE.json')!=dict(entry,prepared_sha256=digest(target/'PREPARED.json')):raise ValueError('Native case belongs to another prepared plan')
    for name in ('matrix_correlators.csv','rank_blocks.csv'):
        if digest(out/name)!=summary['outputs_sha256'][name]:raise ValueError('Changed validated native export')
    return out


def replay_case(out,entry,reference):
    cfg=entry['configuration'];launches=read(out/'LAUNCHES.json');size=len(design.audit_times(cfg))
    expected=[('audit',0,None,'audit%d'%size)]+[('scan',s['piter'],s['t0'],s['label']) for s in entry['scans']]
    if [(r['kind'],r['piter'],r['t0'],r['stem']) for r in launches]!=expected:raise ValueError('Incomplete native launch order/inventory')
    for record in launches:
        if record['configuration']!=cfg or record['exit_code']!=0 or not math.isfinite(record['seconds']) or record['seconds']<=0:raise ValueError('Invalid native MPI launch result')
        for suffix in ('out','err'):
            if digest(out/(record['stem']+'.'+suffix))!=record[suffix+'_sha256']:raise ValueError('Archived native stdout/stderr changed')
    quality=check.audit(out/('audit%d.out'%size),cfg);values=[];blocks=[];scans={};summaries=[]
    for spec in entry['scans']:
        rows,parts,summary=check.scan(out/(spec['label']+'.out'),cfg,spec,reference)
        values+=rows;blocks+=parts;scans[spec['label']]=rows;summaries.append(summary)
    control=check.projection_pair(scans['shift122_p10'],scans['shift122_p40']) if cfg==650 else None
    return values,blocks,dict(configuration=cfg,audit=quality,scans=summaries,p10_p40_maximum_component_difference=control)


def live_previous(target):
    receipts=sorted(target.glob('SUBMISSION_*.json'))
    if not receipts:return
    jobs=[]
    for path in receipts:
        record=read(path);match=re.fullmatch(r'(\d+)(?:;\S+)?\n?',record.get('stdout',''))
        if record.get('returncode')!=0 or not match:raise ValueError('Ambiguous prior submission; inspect its preserved receipt')
        jobs.append(match.group(1))
    result=subprocess.run(['squeue','-h','-u',os.environ['USER'],'-o','%A %F'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if result.returncode or set(jobs)&set(re.findall(r'\d+',result.stdout)):raise ValueError('Prior native submission may still be active; no duplicate submitted')


def submit(args):
    target,data=prepared(args.name,physical=True)
    if not 1<=args.concurrent<=3:raise ValueError('Native pilot concurrency must be 1..3')
    live_previous(target);indices=[]
    for entry in data['cases']:
        if (target/('case_%03d'%entry['index'])/'RUNNING.json').exists():raise ValueError('Preserved running lock; inspect job state before retry')
        if complete_case(target,entry) is None:indices.append(entry['index'])
    if not indices:print('All native pilot cases already validated; collect them');return
    command=['sbatch','--parsable','--array='+','.join(map(str,indices))+'%'+str(args.concurrent),str(HERE/'submit_micro.sh'),args.name]
    receipt=target/('SUBMISSION_%03d.json'%(len(list(target.glob('SUBMISSION_*.json')))+1))
    write(receipt,dict(state='submitting',command=command),exclusive=True)
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    write(receipt,dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?',result.stdout):raise ValueError('Native submission failed or returned an ambiguous job ID')
    print('SUBMITTED',result.stdout.strip(),len(indices),'source-time tasks; concurrency',args.concurrent)


def recover_lock(args):
    target,data=prepared(args.name)
    if not 0<=args.index<3:raise ValueError('Unknown native case index')
    case=target/('case_%03d'%args.index);lock=case/'RUNNING.json';record=read(lock)
    job=str(record['job_id'])
    if not re.fullmatch(r'\d+',job) or record['array_task_id']!=args.index:raise ValueError('Invalid running-lock identity')
    result=subprocess.run(['squeue','-h','-j',job,'-o','%i'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if result.returncode or result.stdout.strip():raise ValueError('Cannot establish that the recorded job has stopped; lock retained')
    if read(lock)!=record:raise ValueError('Running lock changed during recovery')
    number=len(list(case.glob('RECOVERED_RUNNING_*.json')))+1
    write(case/('RECOVERED_RUNNING_%04d.json'%number),dict(record,squeue_stdout=result.stdout,squeue_stderr=result.stderr),exclusive=True)
    lock.unlink();print('Recovered stopped native running lock; attempt files retained.')


def launch(out,entry,launches,spec=None):
    cfg=entry['configuration'];source=entry['input'];size=len(design.audit_times(cfg))
    kind='audit' if spec is None else 'scan';stem='audit%d'%size if spec is None else spec['label']
    piter=0 if spec is None else spec['piter'];t0=None if spec is None else spec['t0']
    command=['mpiexec','-n','64',str(BUILD/('audit%d'%size if spec is None else 'scan'))]
    evroot=str(Path(source['vectors'][0]['path']).parent)
    if spec is None:command += [str(cfg),source['gauge']['path'],evroot,source['eigenvalues']['path']]
    env=os.environ.copy();env.update(QX_CFG=str(cfg),HYP_PITER=str(piter),QX_T0=str(t0 or 0),EVROOT=evroot,QX_ROOT=str(Path(source['gauge']['path']).parent))
    start=time.monotonic()
    with (out/(stem+'.out')).open('x') as stdout,(out/(stem+'.err')).open('x') as stderr:
        result=subprocess.run(command,env=env,stdout=stdout,stderr=stderr)
    record=dict(kind=kind,piter=piter,t0=t0,configuration=cfg,stem=stem,seconds=time.monotonic()-start,exit_code=result.returncode,
        command=command,job_id=os.environ['SLURM_JOB_ID'],array_task_id=os.environ['SLURM_ARRAY_TASK_ID'],
        out_sha256=digest(out/(stem+'.out')),err_sha256=digest(out/(stem+'.err')))
    launches.append(record);write(out/'LAUNCHES.json',launches)
    if result.returncode:raise ValueError('Native MPI launch failed: '+stem)


def run(args):
    target,data=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:raise ValueError('Eight-node, 64-rank Slurm allocation required')
    index=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0<=index<3:raise ValueError('Unknown source-time array task')
    entry=data['cases'][index];validate_inputs(entry);reference=check.load_references(REFERENCE)
    old=complete_case(target,entry)
    if old:
        values,blocks,summary=replay_case(old,entry,reference);csv_equal(old/'matrix_correlators.csv',values);csv_equal(old/'rank_blocks.csv',blocks)
        print('REUSED validated native case',entry['configuration']);return
    case=target/('case_%03d'%index);case.mkdir(exist_ok=True);lock=case/'RUNNING.json'
    with lock.open('x') as stream:json.dump(dict(job_id=os.environ['SLURM_JOB_ID'],array_task_id=index),stream)
    try:
        numbers=[int(p.name.split('_')[1]) for p in case.glob('attempt_*') if re.fullmatch(r'attempt_\d{4}',p.name)]
        out=case/('attempt_%04d'%(max(numbers,default=0)+1));out.mkdir();launches=[]
        write(out/'CASE.json',dict(entry,prepared_sha256=digest(target/'PREPARED.json')),exclusive=True)
        try:
            launch(out,entry,launches);check.audit(out/('audit%d.out'%len(design.audit_times(entry['configuration']))),entry['configuration'])
            for spec in entry['scans']:
                launch(out,entry,launches,spec);check.scan(out/(spec['label']+'.out'),entry['configuration'],spec,reference)
            values,blocks,summary=replay_case(out,entry,reference)
            csv_write(out/'matrix_correlators.csv',values);csv_write(out/'rank_blocks.csv',blocks);validate_inputs(entry)
            write(out/'SUMMARY.json',dict(summary,success=True,case_sha256=digest(out/'CASE.json'),launches_sha256=digest(out/'LAUNCHES.json'),
                outputs_sha256={name:digest(out/name) for name in ('matrix_correlators.csv','rank_blocks.csv')}),exclusive=True)
            write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=digest(out/'SUMMARY.json')))
        except Exception as exc:
            write(out/'FAILED.json',dict(success=False,error=str(exc),launches=launches));raise
    finally:lock.unlink()
    print('PASS source-time native case',entry['configuration'])


def status(args):
    target,data=prepared(args.name);done=[e['configuration'] for e in data['cases'] if complete_case(target,e) is not None]
    print('Validated source-time cases: {}/3; missing {}'.format(len(done),[c for c in PILOT if c not in done]))


def collect(args):
    target,data=prepared(args.name);reference=check.load_references(REFERENCE);values=[];blocks=[];cases=[];missing=[];timings=[]
    if any(target.glob('case_*/RUNNING.json')):raise ValueError('Native tasks have running locks; use status and resolve stopped-job locks before collection')
    for entry in data['cases']:
        out=complete_case(target,entry)
        if out is None:missing.append(entry['configuration']);continue
        rows,parts,summary=replay_case(out,entry,reference);csv_equal(out/'matrix_correlators.csv',rows);csv_equal(out/'rank_blocks.csv',parts)
        values+=rows;blocks+=parts;cases.append(summary);timings.append(dict(configuration=entry['configuration'],launches=read(out/'LAUNCHES.json')))
    projections=[];contrasts=[]
    if values:
        csv_write(target/'matrix_correlators.csv',values);csv_write(target/'rank_blocks.csv',blocks)
        projections,contrasts=derive(values,reference,read(SOURCES))
        csv_write(target/'fixed_source_correlators.csv',projections);csv_write(target/'paired_fixed_source_contrasts.csv',contrasts)
    success=not missing
    write(target/'SUMMARY.json',dict(success=success,validated_configurations=[r['configuration'] for r in cases],missing_configurations=missing,
        matrix_rows=len(values),rank_block_rows=len(blocks),fixed_source_rows=len(projections),paired_fixed_source_rows=len(contrasts),cases=cases,timing=timings,
        statistical_configurations=len(cases),source_training_configurations=68,source_vectors_sha256=digest(SOURCES),
        interpretation='Technical source-time pilot only. No ensemble correlated fit, final potential, boundary-effect significance claim or automatic full71 submission.'))
    archive=ROOT/('qx_static_source_time_'+args.name+'_review.zip');temporary=archive.with_name(archive.name+'.pending')
    paths=[p for p in target.rglob('*') if p.is_file()]
    paths += [HERE/name for name in read(HERE/'FILE_MANIFEST.json')['files']]+[HERE/'FILE_MANIFEST.json']
    paths += [DESIGN/name for name in read(DESIGN/'FILE_MANIFEST.json')['files']]+[DESIGN/'FILE_MANIFEST.json']
    paths += [p for p in BUILD.iterdir() if p.is_file() and (p.suffix in ('.c','.h','.txt') or p.name=='MANIFEST.json')]
    for p in paths:
        if p.is_symlink():raise ValueError('Linked native output: '+str(p))
    with zipfile.ZipFile(temporary,'x',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(set(paths)):
            z.write(p,str(p.relative_to(ROOT)))
    with temporary.open('rb') as stream:os.fsync(stream.fileno())
    temporary.replace(archive)
    print(('PASS' if success else 'INCOMPLETE')+' source-time pilot; validated '+str(len(cases))+'/3; missing '+str(missing))
    print('REVIEW_ARCHIVE',archive)
    if not success:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action');sub.required=True
    for name in ('check-inputs','record-build'):sub.add_parser(name)
    for name in ('prepare','submit','run','status','collect','recover-lock'):
        p=sub.add_parser(name);p.add_argument('--name',required=True)
        if name=='submit':p.add_argument('--concurrent',type=int,default=3)
        if name=='recover-lock':p.add_argument('--index',type=int,required=True)
    args=parser.parse_args()
    if args.action=='check-inputs':verify_sources()
    elif args.action=='record-build':record_build()
    elif args.action=='recover-lock':recover_lock(args)
    else:globals()[args.action](args)


if __name__=='__main__':main()
