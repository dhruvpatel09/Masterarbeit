#!/usr/bin/env python3
"""Full71 orchestration using the unchanged, validated static-matrix binary."""
import argparse
import contextlib
import csv
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = list(range(100,801,10))
PILOT = [170,640,650]
NEW = [c for c in TARGET if c not in PILOT]
MATRIX_ROWS = 16*33*4*4
BLOCK_ROWS = 16*2*2*3*64
NATIVE = ROOT/'laplace_static/validation/qx_static_matrix'
BASE = None


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path,value,exclusive=False):
    path=Path(path)
    if exclusive and path.exists():raise ValueError('Receipt already exists: '+str(path))
    pending=path.with_name(path.name+'.pending')
    with pending.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
        stream.flush();os.fsync(stream.fileno())
    pending.replace(path)


def geometry():
    return dict(t0=90,taus=list(range(33)),R=[4,12,24,32],
                Nv_src=[40,60,80,100],Nv_sink=[40,60,80,100],saved_Nv=100)


def verify_sources():
    manifest=read(HERE/'FILE_MANIFEST.json')
    for group,base in (('files',HERE),('dependencies',ROOT)):
        for name,receipt in manifest[group].items():
            p=base/name
            if p.stat().st_size!=receipt['bytes'] or digest(p)!=receipt['sha256']:
                raise ValueError('Pinned workflow/source changed: '+name)
    source=read(HERE/'INPUTS.json')
    if source['target_configurations']!=TARGET or source['pilot_configurations']!=PILOT or source['geometry']!=geometry():
        raise ValueError('Full71 input specification changed')


def load_base():
    """Keep the native pilot module and its generic import names isolated."""
    global BASE
    verify_sources()
    if BASE is None:
        names=('check_outputs','matrix_plan','measurement_plan')
        old={name:sys.modules.get(name) for name in names}
        previous=sys.path[:]
        for name in names:sys.modules.pop(name,None)
        sys.path.insert(0,str(NATIVE))
        try:
            spec=importlib.util.spec_from_file_location('qx_matrix_full71_native',str(NATIVE/'campaign.py'))
            module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
            BASE=module
        finally:
            sys.path[:]=previous
            for name in names:
                sys.modules.pop(name,None)
                if old[name] is not None:sys.modules[name]=old[name]
    return BASE


def folder(name):
    if not re.fullmatch(r'[a-z][a-z0-9_]{0,40}',name):
        raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def slots():
    return [list(range(i,min(i+2,len(NEW)))) for i in range(0,len(NEW),2)]


def csv_rows(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def replay(out,entry,references):
    """Reparse unchanged native checks and require exact per-case CSV agreement."""
    base=load_base()
    rows,blocks,summary=base.replay_case(out,entry,references)
    if csv_rows(out/'matrix_correlators.csv')!=rows or csv_rows(out/'rank_blocks.csv')!=blocks:
        raise ValueError('Per-case CSV differs from validated stdout replay')
    return rows,blocks,summary


def admit_pilot(name,physical=False,reparse=False):
    base=load_base()
    target,plan=base.prepared(name,physical=physical)
    published=read(HERE/'INPUTS.json')['pilot_receipts']
    if name!=published['name'] or digest(target/'PREPARED.json')!=published['prepared_sha256']:
        raise ValueError('Pilot preparation differs from the published job5961383 receipt')
    if (plan.get('stage')!='technical_pilot' or plan.get('geometry')!=geometry()
            or [e['configuration'] for e in plan['cases']]!=PILOT):
        raise ValueError('Pilot geometry or configuration inventory differs')
    references=base.check.load_references(base.REFERENCE,PILOT) if reparse else None
    reused=[]
    for entry in plan['cases']:
        expected=[10,40] if entry['configuration']==650 else [10]
        if entry['piters']!=expected:raise ValueError('Pilot controls differ')
        out=base.complete_case(target,entry)
        if out is None:raise ValueError('Pilot configuration is not validated: '+str(entry['configuration']))
        launches=read(out/'LAUNCHES.json')
        if [(r['kind'],r['piter'],r['exit_code']) for r in launches]!=[('audit',0,0)]+[('scan',p,0) for p in expected]:
            raise ValueError('Pilot launch inventory changed')
        for launch in launches:
            for suffix in ('out','err'):
                if digest(out/(launch['stem']+'.'+suffix))!=launch[suffix+'_sha256']:
                    raise ValueError('Validated pilot stdout/stderr changed')
        if reparse:replay(out,entry,references)
        reused.append(dict(configuration=entry['configuration'],
            attempt=out.relative_to(ROOT).as_posix(),summary_sha256=digest(out/'SUMMARY.json'),
            case_sha256=digest(out/'CASE.json'),launches_sha256=digest(out/'LAUNCHES.json')))
    if reused!=published['cases']:
        raise ValueError('Pilot case receipts differ from the published job5961383 evidence')
    return target,plan,reused


def prepare(args):
    base=load_base();build=base.build_manifest()
    pilot,pp,reused=admit_pilot(args.pilot,physical=True,reparse=True)
    if pp['build_manifest_sha256']!=build:raise ValueError('Pilot uses a different native build')
    base.check.load_references(base.REFERENCE,TARGET)
    cases=[dict(index=i,configuration=cfg,piters=[10],input=base.input_record(cfg))
           for i,cfg in enumerate(NEW)]
    target=folder(args.name);target.mkdir(parents=True,exist_ok=False)
    plan=dict(stage='full71_matrix',target_configurations=TARGET,reused_configurations=PILOT,
        cases=cases,slots=slots(),geometry=geometry(),piter=10,cases_per_task=2,
        native_build_manifest_sha256=build,workflow_manifest_sha256=digest(HERE/'FILE_MANIFEST.json'),
        pilot=dict(name=args.pilot,prepared_sha256=digest(pilot/'PREPARED.json'),cases=reused),
        reference_sha256=digest(base.REFERENCE))
    write(target/'PREPARED.json',plan,exclusive=True)
    (target/'prepared.sha256').write_text(digest(target/'PREPARED.json')+'\n')
    print('Prepared full71 matrix: 68 new configurations; 3 validated pilot configurations reused; 34 array tasks, two cases each.')
    print('Unchanged p10 binary: R=4,12,24,32; all 16 ordered source/sink pairs; correlator times 0..32.')


def prepared(name,physical=False):
    base=load_base();target=folder(name);plan=read(target/'PREPARED.json')
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():
        raise ValueError('Prepared full71 plan checksum mismatch')
    if plan.get('stage')!='full71_matrix' or plan['target_configurations']!=TARGET or plan['reused_configurations']!=PILOT:
        raise ValueError('Full71 cohort changed')
    if plan['geometry']!=geometry() or plan['piter']!=10 or plan['cases_per_task']!=2 or plan['slots']!=slots():
        raise ValueError('Full71 geometry or slot inventory changed')
    if plan['workflow_manifest_sha256']!=digest(HERE/'FILE_MANIFEST.json'):
        raise ValueError('Orchestration differs from prepared plan')
    if plan['native_build_manifest_sha256']!=base.build_manifest():
        raise ValueError('Native build differs from prepared plan')
    if plan['reference_sha256']!=digest(base.REFERENCE):
        raise ValueError('Published diagonal reference changed')
    if len(plan['cases'])!=len(NEW):raise ValueError('Wrong number of new configurations')
    for i,(entry,cfg) in enumerate(zip(plan['cases'],NEW)):
        if (set(entry)!=set(('index','configuration','piters','input')) or
                entry['index']!=i or entry['configuration']!=cfg or entry['piters']!=[10]
                or entry['input']['configuration']!=cfg):
            raise ValueError('Invalid full71 case inventory')
        if physical:base.validate_inputs(entry['input'])
    pilot,_,reused=admit_pilot(plan['pilot']['name'],physical=physical,reparse=False)
    if plan['pilot']['prepared_sha256']!=digest(pilot/'PREPARED.json') or plan['pilot']['cases']!=reused:
        raise ValueError('Validated pilot reuse receipts changed')
    return target,plan


def selected_slots(args,plan):
    start=args.start
    count=len(plan['slots'])-start if args.count is None else args.count
    if not isinstance(start,int) or not isinstance(count,int) or start<0 or count<1 or start+count>len(plan['slots']):
        raise ValueError('Invalid array task range')
    return list(range(start,start+count))


def status(args):
    base=load_base();target,plan=prepared(args.name)
    done={e['index']:base.complete_case(target,e) is not None for e in plan['cases']}
    completed=sum(done.values());finished=sum(all(done[i] for i in slot) for slot in plan['slots'])
    print('Validated markers: {} new + 3 reused = {}/71 configurations; {}/34 array tasks complete.'.format(completed,completed+3,finished))
    pending=[e['configuration'] for e in plan['cases'] if not done[e['index']]]
    print('Missing configurations:',pending)
    for lock in sorted(target.glob('case_*/RUNNING.json')):
        print('RUNNING receipt:',lock.relative_to(ROOT).as_posix(),json.dumps(read(lock),sort_keys=True))


def submit(args):
    base=load_base();target,plan=prepared(args.name,physical=True)
    if not 1<=args.concurrent<=4:raise ValueError('Concurrency must be 1..4')
    base.live_previous(target)
    indices=[s for s in selected_slots(args,plan)
             if any(base.complete_case(target,plan['cases'][i]) is None for i in plan['slots'][s])]
    if not indices:
        print('Selected tasks are validated already; collect the campaign.');return
    for s in indices:
        for i in plan['slots'][s]:
            lock=target/('case_%03d'%i)/'RUNNING.json'
            if lock.exists() and base.complete_case(target,plan['cases'][i]) is None:
                raise ValueError('Existing RUNNING receipt; inspect its recorded job before resubmitting: '+str(lock))
    command=['sbatch','--parsable','--array='+','.join(map(str,indices))+'%'+str(args.concurrent),
             str(HERE/'submit_micro.sh'),args.name]
    receipt=target/('SUBMISSION_%03d.json'%(len(list(target.glob('SUBMISSION_*.json')))+1))
    write(receipt,dict(state='submitting',command=command),exclusive=True)
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    write(receipt,dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr))
    if result.returncode or not re.fullmatch(r'\d+(?:;\S+)?\n?',result.stdout):
        raise ValueError('Submission failed or returned an ambiguous job ID')
    print('SUBMITTED',result.stdout.strip(),len(indices),'full71 matrix tasks; at most',args.concurrent,'concurrent')


@contextlib.contextmanager
def case_lock(case):
    path=case/'RUNNING.json'
    try:
        with path.open('x') as stream:
            json.dump(dict(job_id=os.environ.get('SLURM_JOB_ID'),
                array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'),pid=os.getpid()),stream)
            stream.write('\n');stream.flush();os.fsync(stream.fileno())
    except FileExistsError:
        raise ValueError('Case has a RUNNING receipt: '+str(path))
    try:yield
    finally:path.unlink()


def run_case(target,entry,references):
    base=load_base();old=base.complete_case(target,entry)
    if old is not None:
        replay(old,entry,references);print('REUSED validated full71 n%d'%entry['configuration']);return
    base.validate_inputs(entry['input'])
    case=target/('case_%03d'%entry['index']);case.mkdir(exist_ok=True)
    with case_lock(case):
        old=base.complete_case(target,entry)
        if old is not None:
            replay(old,entry,references);return
        numbers=[int(p.name.split('_')[1]) for p in case.glob('attempt_*') if re.fullmatch(r'attempt_\d{4}',p.name)]
        out=case/('attempt_%04d'%(max(numbers,default=0)+1));out.mkdir();launches=[]
        write(out/'CASE.json',dict(entry,prepared_sha256=digest(target/'PREPARED.json')),exclusive=True)
        try:
            base.launch(out,'audit',0,entry,launches)
            base.check.audit(out/'audit.out',entry['configuration'])
            base.launch(out,'scan',10,entry,launches)
            rows,blocks,summary=base.replay_case(out,entry,references)
            if len(rows)!=MATRIX_ROWS or len(blocks)!=BLOCK_ROWS:
                raise ValueError('Native output inventory differs from the validated matrix geometry')
            base.csv_write(out/'matrix_correlators.csv',rows)
            base.csv_write(out/'rank_blocks.csv',blocks)
            base.validate_inputs(entry['input'])
            write(out/'SUMMARY.json',dict(summary,success=True,
                launches_sha256=digest(out/'LAUNCHES.json'),case_sha256=digest(out/'CASE.json'),
                outputs_sha256={n:digest(out/n) for n in ('matrix_correlators.csv','rank_blocks.csv')}),exclusive=True)
            write(case/'LATEST.json',dict(attempt=out.name,summary_sha256=digest(out/'SUMMARY.json')))
        except Exception as exc:
            write(out/'FAILED.json',dict(success=False,error=str(exc),launches=launches));raise
    print('PASS full71 matrix n%d'%entry['configuration'])


def run(args):
    base=load_base();target,plan=prepared(args.name)
    if int(os.environ.get('SLURM_NTASKS','0'))!=64 or int(os.environ.get('SLURM_JOB_NUM_NODES','0'))!=8:
        raise ValueError('Eight nodes and 64 MPI ranks required')
    slot=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0<=slot<len(plan['slots']):raise ValueError('Unknown full71 array task')
    entries=[plan['cases'][i] for i in plan['slots'][slot]]
    refs=base.check.load_references(base.REFERENCE,[e['configuration'] for e in entries])
    for entry in entries:run_case(target,entry,refs)


class CSVStream:
    def __init__(self,path):
        self.path=Path(path);self.pending=self.path.with_name(self.path.name+'.pending')
        self.stream=self.pending.open('x',newline='')
        self.writer=None;self.fields=None;self.count=0
    def add(self,rows):
        if not rows:return
        if self.writer is None:
            self.fields=list(rows[0]);self.writer=csv.DictWriter(self.stream,fieldnames=self.fields,lineterminator='\n')
            self.writer.writeheader()
        if any(list(row)!=self.fields for row in rows):raise ValueError('Merged CSV schema changed')
        self.writer.writerows(rows);self.count+=len(rows)
    def commit(self):
        self.stream.flush();os.fsync(self.stream.fileno());self.stream.close();self.pending.replace(self.path)
    def abort(self):
        self.stream.close()
        if self.pending.exists():self.pending.unlink()


def resolved_cases(target,plan):
    base=load_base();result={}
    for entry in plan['cases']:
        out=base.complete_case(target,entry)
        if out is None:raise ValueError('Incomplete full71 campaign; missing n%d; use status'%entry['configuration'])
        result[entry['configuration']]=(out,entry,False)
    pilot,pp,reused=admit_pilot(plan['pilot']['name'],reparse=False)
    for entry in pp['cases']:
        out=base.complete_case(pilot,entry)
        result[entry['configuration']]=(out,entry,True)
    if sorted(result)!=TARGET:raise ValueError('Assembled cohort is not exactly full71')
    return result,pilot


def review_archive(target,pilot):
    base=load_base();files=set();omitted=[];native_omitted=[]
    source=read(HERE/'INPUTS.json')
    candidates=list(target.rglob('*'))+list(pilot.rglob('*'))
    candidates += [HERE/name for name in read(HERE/'FILE_MANIFEST.json')['files']]+[HERE/'FILE_MANIFEST.json']
    published_pilot=ROOT/source['published_pilot_archive']
    candidates += [ROOT/name for name in source['files'] if ROOT/name not in (base.REFERENCE,published_pilot)]
    candidates += [p for p in base.BUILD.iterdir()
                   if p.is_file() and (p.suffix in ('.h','.txt') or p.name in ('MANIFEST.json','audit_selected.c'))]
    for name,expected in read(base.BUILD/'MANIFEST.json')['files'].items():
        p=ROOT/name
        if p.suffix in ('.h','.c','.txt','.json'):
            candidates.append(p)
        else:
            native_omitted.append(dict(path=name,sha256=expected,reason='native executable or library; original build receipt retained'))
    for p in candidates:
        if p.is_symlink():raise ValueError('Linked provenance file: '+str(p))
        if not p.is_file() or p.name=='ARCHIVE_FILES.json':continue
        if p.name.endswith('.pending') or p.name=='RUNNING.json':
            omitted.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p),reason='incomplete or running receipt'))
        else:files.add(p)
    inventory=dict(files=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(files)],
        omitted=omitted,published_diagonal_reference=dict(path=base.REFERENCE.relative_to(ROOT).as_posix(),sha256=digest(base.REFERENCE)),
        published_pilot_archive=dict(path=published_pilot.relative_to(ROOT).as_posix(),sha256=digest(published_pilot)),
        omitted_native_inputs=native_omitted,
        native_binary_policy='Native binary/library hashes are in the unchanged build manifest; executables and libraries are not copied.',
        inventory_self_excluded=True)
    write(target/'ARCHIVE_FILES.json',inventory);files.add(target/'ARCHIVE_FILES.json')
    archive=ROOT/('qx_static_matrix_'+target.name+'_review.zip')
    pending=archive.with_name(archive.name+'.pending')
    if pending.exists():raise ValueError('Unfinished prior archive: '+str(pending))
    with zipfile.ZipFile(pending,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(files):z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(pending) as z:
        if z.testzip() is not None:raise ValueError('Review archive CRC check failed')
    with pending.open('rb') as stream:os.fsync(stream.fileno())
    pending.replace(archive)
    return archive


def collect(args):
    base=load_base();target,plan=prepared(args.name)
    if list(target.glob('case_*/RUNNING.json')):
        raise ValueError('RUNNING receipts remain; inspect status before collecting')
    resolved,pilot=resolved_cases(target,plan)
    refs=base.check.load_references(base.REFERENCE,TARGET)
    matrices=CSVStream(target/'matrix_correlators.csv')
    try:blocks=CSVStream(target/'rank_blocks.csv')
    except Exception:matrices.abort();raise
    summaries=[];timing=[]
    try:
        for cfg in TARGET:
            out,entry,reused=resolved[cfg]
            rows,parts,summary=replay(out,entry,refs)
            primary=[r for r in rows if r['piter']=='10']
            rank_primary=[r for r in parts if r['piter']=='10']
            if len(primary)!=MATRIX_ROWS or len(rank_primary)!=BLOCK_ROWS:
                raise ValueError('Incomplete p10 production geometry at n%d'%cfg)
            matrices.add(primary);blocks.add(rank_primary)
            summaries.append(dict(summary,reused_from_pilot=reused,attempt=out.relative_to(ROOT).as_posix()))
            timing.append(dict(configuration=cfg,reused_from_pilot=reused,launches=read(out/'LAUNCHES.json')))
        if matrices.count!=71*MATRIX_ROWS or blocks.count!=71*BLOCK_ROWS:
            raise ValueError('Full71 merged row count changed')
        matrices.commit();blocks.commit()
    except Exception:
        matrices.abort();blocks.abort();raise
    summary=dict(success=True,configuration_labels=TARGET,reused_configurations=PILOT,new_configurations=NEW,
        validated_configurations=TARGET,missing_configurations=[],geometry=geometry(),piter=10,
        matrix_rows=matrices.count,rank_block_rows=blocks.count,
        prepared_sha256=digest(target/'PREPARED.json'),workflow_manifest_sha256=digest(HERE/'FILE_MANIFEST.json'),
        native_build_manifest_sha256=plan['native_build_manifest_sha256'],
        outputs_sha256={name:digest(target/name) for name in ('matrix_correlators.csv','rank_blocks.csv')},
        cases=summaries,timing=timing,
        interpretation='71 paired configuration matrices, retaining all ordered complex entries; no automatic Hermitian projection or physical energy extraction.')
    write(target/'SUMMARY.json',summary)
    archive=review_archive(target,pilot)
    print('PASS full71 matrix: 71 configurations; 68 new + 3 reused; {} p10 matrix rows; {} rank-block rows.'.format(matrices.count,blocks.count))
    print('REVIEW_ARCHIVE',archive)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action');sub.required=True
    sub.add_parser('check-inputs');sub.add_parser('plan')
    for action in ('prepare','submit','run','collect','status'):
        p=sub.add_parser(action);p.add_argument('--name',required=True)
        if action=='prepare':p.add_argument('--pilot',default='pilot')
        if action=='submit':
            p.add_argument('--start',type=int,default=0);p.add_argument('--count',type=int)
            p.add_argument('--concurrent',type=int,default=4)
    args=parser.parse_args()
    if args.action=='check-inputs':
        verify_sources();load_base().verify_sources()
        print('PASS full71 matrix workflow and unchanged native source/reference pins')
    elif args.action=='plan':
        print(json.dumps(dict(target_configurations=TARGET,reused_configurations=PILOT,new_configurations=NEW,
            array_tasks=len(slots()),cases_per_task=2,nodes_per_task=8,ranks_per_task=64,
            maximum_concurrent_tasks=4,geometry=geometry()),indent=2))
    else:globals()[args.action](args)


if __name__=='__main__':main()
