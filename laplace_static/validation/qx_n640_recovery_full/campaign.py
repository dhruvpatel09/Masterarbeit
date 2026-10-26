#!/usr/bin/env python3
"""Complete n640 using disjoint saved-mode blocks and a full read-back audit."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BUILD=HERE/'build'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def package_pins():
    result={}
    for p in sorted(HERE.rglob('*')):
        rel=p.relative_to(HERE)
        if not p.is_file() or '__pycache__' in rel.parts:continue
        if rel.parts[0]!='inputs' and (any(x in ('build','campaigns') for x in rel.parts) or p.suffix in ('.out','.err','.log','.o','.zip')):continue
        result[str(rel)]=sha(p)
    return result


def pilot_module():
    path=HERE/'inputs/pilot/workflow/campaign.py'
    spec=importlib.util.spec_from_file_location('qx_recovery_saved_pilot',str(path))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def pilot_replay():
    folder=HERE/'inputs/pilot/campaign';m=pilot_module()
    plan=read(folder/'PREPARED.json');launches=read(folder/'LAUNCHES.json')
    if [r['label'] for r in launches]!=['generate','audit']:raise ValueError('Pilot launch inventory changed')
    for row in launches:
        if row['exit_code']!=0:raise ValueError('Pilot MPI launch failed')
        for ext,key in (('out','stdout_sha256'),('err','stderr_sha256')):
            if sha(folder/'outputs'/(row['label']+'.'+ext))!=row[key]:raise ValueError('Pilot output changed')
    summary=m.replay_logs((folder/'outputs/generate.out').read_text(),(folder/'outputs/audit.out').read_text(),plan['selected_times'])
    if summary!=read(folder/'SUMMARY.json'):raise ValueError('Pilot summary is not reproducible')
    return plan,summary


def check_spec(spec):
    reused=[0,1,32,64,95,96,97,98,99,100,101,128,160,191]
    expected=[[t for t in range(lo,lo+12) if t not in reused] for lo in range(0,192,12)]
    if (spec['configuration']!=640 or spec['times']!=list(range(192)) or spec['reused_times']!=reused or
        spec['slots']!=expected or spec['parameters']['res']!=1e-9):raise ValueError('Full recovery scope changed')


def check_inputs():
    spec=read(HERE/'INPUTS.json');check_spec(spec)
    for name,digest in spec['preserved_inputs'].items():
        if sha(HERE/name)!=digest:raise ValueError('Preserved input changed: '+name)
    for name,digest in spec['repository_references'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Production dependency changed: '+name)
    plan,summary=pilot_replay()
    if plan!=spec['pilot_plan']:raise ValueError('Pilot plan changed')
    return spec


def record_build():
    check_inputs()
    names=['mental_recovery','check_saved','assemble_values','mental_recovery.c','check_saved.c','selection.h',
        'mental_recovery.o','check_saved.o','assemble_values.o','DistRedelemsHdf5Writer.o','DistEigvecsHdf5Writer.o',
        'DistEigvecsHdf5Reader.o','DistEigvalsHdf5Writer.o','linked_generator.txt','linked_checker.txt',
        'linked_assembler.txt','modules.txt','compiler.txt']
    write(BUILD/'MANIFEST.json',{'files':{n:sha(BUILD/n) for n in names},'package_files':package_pins()})


def verify_build():
    check_inputs();manifest=read(BUILD/'MANIFEST.json')
    if package_pins()!=manifest['package_files']:raise ValueError('Package changed after build')
    for name,digest in manifest['files'].items():
        if sha(BUILD/name)!=digest:raise ValueError('Build changed: '+name)


def stat_record(path):
    p=Path(path);st=p.stat()
    return {'path':str(p),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}


def check_saved_inventory(inventory):
    for row in inventory:
        if stat_record(row['path'])!=row:raise ValueError('Saved data changed or are missing: '+row['path'])


def directory(name):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def prepare(args):
    spec=check_inputs();verify_build();folder=directory(args.name)
    if folder.exists():raise ValueError('Campaign already exists; preserve it')
    pilot=spec['pilot_plan'];m=pilot_module()
    if m.gauge_record(pilot['gauge']['path'])!=pilot['gauge']:raise ValueError('Gauge changed since pilot')
    inventory=read(HERE/'inputs/pilot/campaign/SAVED_FILES.json');check_saved_inventory(inventory)
    storage=Path(spec['storage_root']);storage.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(str(storage)).free<230*1024**3:raise ValueError('Need at least 230 GiB free for the 178 remaining times')
    run=storage/(args.name+'_'+time.strftime('%Y%m%d_%H%M%S',time.gmtime()));run.mkdir()
    folder.mkdir(parents=True);slots=[]
    for index,times in enumerate(spec['slots']):
        slot=folder/('slot_%02d'%index);slot.mkdir()
        saved=run/('slot_%02d'%index);saved.mkdir()
        for name in ('log','dat','eigenvectors','eigenvalues'):(saved/name).mkdir()
        (slot/'mental.in').write_text(m.input_text(saved))
        slots.append({'index':index,'times':times,'run_directory':str(saved),'input_sha256':sha(slot/'mental.in'),
            'eigenvalue_file':str(saved/'eigenvalues/qX_p2_n640_selected_recovery_n640_eval.h5')})
    plan={'name':args.name,'configuration':640,'times':spec['times'],'reused_times':spec['reused_times'],
        'slots':slots,'run_directory':str(run),'gauge':pilot['gauge'],'pilot_inventory':inventory,
        'pilot_eigenvalue_sha256':sha(pilot['eigenvalue_file']),'pilot_plan':pilot,
        'build_manifest_sha256':sha(BUILD/'MANIFEST.json'),'scope':spec['scope']}
    write(folder/'PREPARED.json',plan);(folder/'prepared.sha256').write_text(sha(folder/'PREPARED.json')+'\n')
    write(run/'RECOVERY_SCOPE.json',{'complete_configuration':False,'reused_times':spec['reused_times'],'new_times':178})
    print('Prepared full n640: reuse 14 times; generate 178 in 16 disjoint tasks; final audit over 192 times')


def prepared(name):
    folder=directory(name);plan=read(folder/'PREPARED.json')
    if sha(folder/'PREPARED.json')!=(folder/'prepared.sha256').read_text().strip():raise ValueError('Plan changed')
    if sha(BUILD/'MANIFEST.json')!=plan['build_manifest_sha256']:raise ValueError('Prepared build changed')
    verify_build()
    if pilot_module().gauge_record(plan['gauge']['path'])!=plan['gauge']:raise ValueError('Gauge changed')
    check_saved_inventory(plan['pilot_inventory'])
    if sha(plan['pilot_plan']['eigenvalue_file'])!=plan['pilot_eigenvalue_sha256']:raise ValueError('Pilot eigenvalue container changed')
    for row in plan['slots']:
        if sha(folder/('slot_%02d'%row['index'])/'mental.in')!=row['input_sha256']:raise ValueError('Slot input changed')
    return folder,plan


def sbatch(command):
    r=subprocess.run(command,cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if r.returncode:raise ValueError(r.stderr)
    job=r.stdout.strip().split(';')[0]
    if not job.isdigit():raise ValueError('Unexpected sbatch response: '+r.stdout)
    return job


def submit_final(args):
    folder,plan=prepared(args.name);record=read(folder/'SUBMISSION.json')
    if record.get('final_job'):raise ValueError('Final job already submitted')
    job=sbatch(['sbatch','--parsable','--dependency=afterok:'+record['array_job'],str(HERE/'submit_final.sh'),args.name])
    record['final_job']=job;write(folder/'SUBMISSION.json',record)
    print('SUBMITTED '+job+' final assembly and full read-back audit; waits for array '+record['array_job'])


def submit(args):
    folder,plan=prepared(args.name)
    if not 1<=args.concurrent<=4:raise ValueError('Concurrency must be 1..4')
    if (folder/'SUBMISSION.json').exists():raise ValueError('Array already submitted; do not duplicate it')
    job=sbatch(['sbatch','--parsable','--array=0-15%'+str(args.concurrent),str(HERE/'submit_micro.sh'),args.name])
    write(folder/'SUBMISSION.json',{'array_job':job,'prepared_sha256':sha(folder/'PREPARED.json')})
    print('SUBMITTED '+job+' 16 tasks, 178 new times, concurrency '+str(args.concurrent))
    submit_final(args)


def records(text,kind):
    return [dict(x.split('=',1) for x in line.split()[1:] if '=' in x) for line in text.splitlines() if line.startswith(kind+' ')]


def validate_generation(text,times):
    begin=records(text,'RECOVERY_START');solves=records(text,'RECOVERY_SOLVE');written=records(text,'RECOVERY_WRITTEN')
    expected={'cfg':'640','ranks':'64','P':'1,4,4,4','selected':str(len(times)),'res':'1.00000000000000006e-09'}
    if begin!=[expected]:raise ValueError('Changed generation setup')
    for rows in (solves,written):
        if [int(r['t']) for r in rows]!=times or any(r['cfg']!='640' for r in rows):raise ValueError('Incomplete generated time coverage')
    for r in solves:
        if (r['code']!='0' or not 1<=int(r['nRed'])<=30 or int(r['restarts'])!=int(r['nRed'])-1 or
            not math.isfinite(float(r['maxRes'])) or not 0<=float(r['maxRes'])<=1e-9 or
            not 0<float(r['lambda_Nv'])<0.15):raise ValueError('Generated solve failed original target')
    if any(r['Nv']!='100' for r in written):raise ValueError('Wrong mode count')
    return solves


def validate_audit(text,times,expected_values=None):
    begin=records(text,'QUALITY_START');end=records(text,'QUALITY_END')
    if (len(begin)!=1 or begin[0]['cfg']!='640' or begin[0]['grid']!='1x4x4x4' or begin[0]['Nv']!='100' or
        begin[0]['APE']!='20,0.5' or int(begin[0]['FUNNELED'])<1):raise ValueError('Changed audit setup')
    if end!=[{'cfg':'640','slices':str(len(times)),'status':'PASS'}]:raise ValueError('Saved-mode audit did not pass')
    modes=records(text,'MODE');slices=records(text,'SLICE');values={}
    if len(modes)!=100*len(times) or [int(r['t']) for r in slices]!=times:raise ValueError('Incomplete audit inventory')
    for t,row in zip(times,slices):
        group=[m for m in modes if int(m['t'])==t]
        if [int(m['i']) for m in group]!=list(range(100)):raise ValueError('Duplicate or missing mode')
        for i,m in enumerate(group):
            if m['cfg']!='640' or any(not math.isfinite(float(m[k])) for k in ('lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im')):
                raise ValueError('Nonfinite or wrong saved mode')
            x=float(m['lambda']);norm=float(m['norm']);res=float(m['residual_relative'])
            if not 0<x<0.15 or abs(norm-1)>1e-9 or not 0<=res<=1e-9 or float(m['residual_l2'])<0:
                raise ValueError('Saved pair fails original target')
            if not math.isclose(float(m['residual_l2'])/(x*norm),res,rel_tol=1e-12,abs_tol=1e-16):raise ValueError('Residual formula mismatch')
            if i and x+1e-12<float(group[i-1]['lambda']):raise ValueError('Unordered spectrum')
            if expected_values is not None and x!=expected_values[(t,i)]:raise ValueError('Assembled eigenvalue mismatch')
            values[(t,i)]=x
        if row['cfg']!='640' or row['Nv']!='100' or row['status']!='PASS':raise ValueError('Slice failed')
        limits={'max_relative_residual':1e-9,'max_norm_defect':1e-9,'gram_diag_defect':1e-9,'gram_offdiag':1e-9,'max_rayleigh_im':1e-10}
        for key,limit in limits.items():
            x=float(row[key])
            if not math.isfinite(x) or not 0<=x<=limit:raise ValueError('Audit limit failed: '+key)
        for key,x in [('max_relative_residual',max(float(m['residual_relative']) for m in group)),
                      ('max_norm_defect',max(abs(float(m['norm'])-1) for m in group))]:
            if not math.isclose(float(row[key]),x,rel_tol=1e-12,abs_tol=1e-15):raise ValueError('Audit aggregate mismatch')
        if float(row['lambda_min'])!=values[(t,0)] or float(row['lambda_max'])!=values[(t,99)]:raise ValueError('Audit spectrum range mismatch')
    return values,slices


def launch(folder,label,command,slot,launches):
    output=folder/'outputs';output.mkdir(exist_ok=True)
    if (output/(label+'.out')).exists():raise ValueError('Launch output already exists; preserve the attempt')
    env=dict(os.environ);env['QX_RECOVERY_SLOT']=str(slot);env.pop('MENTAL_EVAL_ONLY',None)
    began=time.time()
    with (output/(label+'.out')).open('w') as out,(output/(label+'.err')).open('w') as err:
        r=subprocess.run(command,cwd=str(ROOT),env=env,stdout=out,stderr=err)
    launches.append({'label':label,'command':command,'slot':slot,'exit_code':r.returncode,'seconds':time.time()-began,
        'stdout_sha256':sha(output/(label+'.out')),'stderr_sha256':sha(output/(label+'.err'))})
    write(folder/'LAUNCHES.json',launches)
    if r.returncode:raise ValueError('MPI launch failed: '+label)


def verify_launches(folder,labels):
    launches=read(folder/'LAUNCHES.json')
    if [r['label'] for r in launches]!=labels:raise ValueError('Incomplete launch inventory')
    for row in launches:
        if row['exit_code']!=0:raise ValueError('Failed MPI launch')
        for ext,key in (('out','stdout_sha256'),('err','stderr_sha256')):
            if sha(folder/'outputs'/(row['label']+'.'+ext))!=row[key]:raise ValueError('Changed MPI output')
    return launches


def slot_inventory(row):
    saved=Path(row['run_directory'])/'eigenvectors'
    expected=['qX_p2n640_evec_t%d.h5'%t for t in row['times']]
    if sorted(p.name for p in saved.glob('*.h5'))!=sorted(expected):raise ValueError('Unexpected slot vector inventory')
    inventory=[stat_record(saved/name) for name in expected]+[stat_record(row['eigenvalue_file'])]
    if any(x['bytes']<1024 for x in inventory):raise ValueError('Empty or truncated HDF5 output')
    return inventory


def replay_slot(folder,row):
    verify_launches(folder,['generate','audit'])
    solves=validate_generation((folder/'outputs/generate.out').read_text(),row['times'])
    values,slices=validate_audit((folder/'outputs/audit.out').read_text(),row['times'])
    for solve in solves:
        if abs(float(solve['lambda_Nv'])-values[(int(solve['t']),99)])>1e-10:raise ValueError('Saved spectrum differs from generation')
    inventory=slot_inventory(row)
    if (folder/'SAVED_FILES.json').exists() and inventory!=read(folder/'SAVED_FILES.json'):raise ValueError('Slot files changed since audit')
    return values,slices,inventory


def run_slot(args):
    folder,plan=prepared(args.name)
    if not 0<=args.slot<len(plan['slots']):raise ValueError('Invalid slot')
    row=plan['slots'][args.slot];slot=folder/('slot_%02d'%args.slot);launches=[]
    launch(slot,'generate',['mpiexec','-n','64','-ppn','8',str(BUILD/'mental_recovery'),'-i',str(slot/'mental.in')],args.slot,launches)
    launch(slot,'audit',['mpiexec','-n','64','-ppn','8',str(BUILD/'check_saved'),'640',plan['gauge']['path'],
        str(Path(row['run_directory'])/'eigenvectors'),row['eigenvalue_file']],args.slot,launches)
    values,slices,inventory=replay_slot(slot,row)
    write(slot/'SAVED_FILES.json',inventory)
    write(slot/'SUMMARY.json',{'success':True,'times':row['times'],'audited_modes':len(values),'slices':slices,'complete_configuration':False})
    print('PASS slot %d: %d saved times; %d modes'%(args.slot,len(row['times']),len(values)))


def joined_data(folder,plan):
    pilot_replay();pilot=HERE/'inputs/pilot/campaign'
    values,slices=validate_audit((pilot/'outputs/audit.out').read_text(),plan['reused_times'])
    paths={t:str(Path(plan['pilot_plan']['run_directory'])/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)) for t in plan['reused_times']}
    streams=[{'source':'pilot','audit_sha256':sha(pilot/'outputs/audit.out'),'times':plan['reused_times']}]
    for row in plan['slots']:
        slot=folder/('slot_%02d'%row['index']);current,current_slices,inventory=replay_slot(slot,row)
        if set(current)&set(values):raise ValueError('Overlapping time sources')
        values.update(current);slices.extend(current_slices)
        paths.update({t:str(Path(row['run_directory'])/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)) for t in row['times']})
        streams.append({'source':'slot_%02d'%row['index'],'audit_sha256':sha(slot/'outputs/audit.out'),'times':row['times']})
    if set(values)!={(t,i) for t in range(192) for i in range(100)} or set(paths)!=set(range(192)):
        raise ValueError('Incomplete full-configuration coverage')
    return values,paths,streams


def assembled_inventory(final):
    return [stat_record(final/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)) for t in range(192)]+[
        stat_record(final/'eigenvalues/qX_p2_n640_full192_recovery_eval.h5')]


def run_final(args):
    folder,plan=prepared(args.name);values,paths,streams=joined_data(folder,plan)
    control=folder/'final';control.mkdir()
    final=Path(plan['run_directory'])/'assembled';final.mkdir()
    (final/'eigenvectors').mkdir();(final/'eigenvalues').mkdir()
    for t,path in paths.items():(final/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t)).symlink_to(path)
    csv_path=control/'eigenvalues_192times.csv'
    with csv_path.open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(['t','i','lambda'])
        for t in range(192):
            for i in range(100):writer.writerow([t,i,'%.17e'%values[(t,i)]])
    source_inventory=[stat_record(paths[t]) for t in range(192)]
    provenance={'sources':streams,'vector_paths':{str(t):p for t,p in paths.items()},'values_sha256':sha(csv_path),
        'template_eigenvalue_file':plan['pilot_plan']['eigenvalue_file'],'template_metadata':'Preserved from pilot; complete numerical table is assembled from per-time saved-mode audits.',
        'complete_configuration':False,'production_configuration_admitted':False}
    write(control/'ASSEMBLY.json',provenance)
    eigenvalues=final/'eigenvalues/qX_p2_n640_full192_recovery_eval.h5'
    shutil.copyfile(plan['pilot_plan']['eigenvalue_file'],eigenvalues)
    launches=[]
    launch(control,'assemble',['mpiexec','-n','1','-ppn','1',str(BUILD/'assemble_values'),str(eigenvalues),str(csv_path),sha(csv_path)],'all',launches)
    launch(control,'audit_all',['mpiexec','-n','64','-ppn','8',str(BUILD/'check_saved'),'640',plan['gauge']['path'],str(final/'eigenvectors'),str(eigenvalues)],'all',launches)
    audit_values,slices=validate_audit((control/'outputs/audit_all.out').read_text(),list(range(192)),values)
    check_saved_inventory(source_inventory)
    write(control/'SAVED_FILES.json',assembled_inventory(final))
    write(control/'SUMMARY.json',{'complete_eigenpair_set_valid':True,'configuration':640,'times':192,'modes_per_time':100,
        'audited_modes':len(audit_values),'slices':slices,'production_configuration_admitted':False,'profile_measurements_added':False,
        'eigenvector_directory':str(final/'eigenvectors'),'eigenvalue_file':str(eigenvalues),
        'eigenvalue_file_sha256':sha(eigenvalues),'values_sha256':sha(csv_path),'assembly_manifest_sha256':sha(control/'ASSEMBLY.json'),
        'max_relative_residual':max(float(r['max_relative_residual']) for r in slices),
        'scope':'Complete saved-mode set, independently read back at all 192 times. Profile admission is a separate step.'})
    write(final/'VALIDATED_RECOVERY.json',read(control/'SUMMARY.json'))
    print('PASS n640 complete saved-mode set: 192 times; 19200 independently audited modes')


def final_replay(folder,plan):
    values,paths,streams=joined_data(folder,plan);control=folder/'final'
    verify_launches(control,['assemble','audit_all'])
    if records((control/'outputs/assemble.out').read_text(),'ASSEMBLED_VALUES')!=[{'times':'192','Nv':'100','rows':'19200'}]:
        raise ValueError('Eigenvalue assembly marker missing')
    audit_values,slices=validate_audit((control/'outputs/audit_all.out').read_text(),list(range(192)),values)
    summary=read(control/'SUMMARY.json')
    if (not summary['complete_eigenpair_set_valid'] or summary['configuration']!=640 or summary['times']!=192 or
        summary['modes_per_time']!=100 or summary['audited_modes']!=19200 or summary['slices']!=slices or
        summary['production_configuration_admitted'] or summary['profile_measurements_added']):raise ValueError('Final summary mismatch')
    if (sha(control/'ASSEMBLY.json')!=summary['assembly_manifest_sha256'] or
        sha(control/'eigenvalues_192times.csv')!=summary['values_sha256'] or
        sha(summary['eigenvalue_file'])!=summary['eigenvalue_file_sha256']):raise ValueError('Assembly inputs or eigenvalue container changed')
    assembly=read(control/'ASSEMBLY.json')
    if (assembly['sources']!=streams or assembly['vector_paths']!={str(t):p for t,p in paths.items()} or
        assembly['values_sha256']!=summary['values_sha256']):raise ValueError('Assembly provenance mismatch')
    with (control/'eigenvalues_192times.csv').open(newline='') as f:
        rows=list(csv.reader(f))
    if rows[0]!=['t','i','lambda'] or len(rows)!=19201:raise ValueError('Incomplete assembled eigenvalue table')
    for index,row in enumerate(rows[1:]):
        t,i=divmod(index,100)
        if len(row)!=3 or int(row[0])!=t or int(row[1])!=i or float(row[2])!=values[(t,i)]:
            raise ValueError('Assembled eigenvalue table differs from source audits')
    check_saved_inventory(read(control/'SAVED_FILES.json'))
    for t,path in paths.items():
        if (Path(summary['eigenvector_directory'])/('qX_p2n640_evec_t%d.h5'%t)).resolve()!=Path(path).resolve():
            raise ValueError('Assembled vector link changed')
    return summary


def collect(args):
    folder=directory(args.name)
    if not folder.is_dir():raise ValueError('Campaign does not exist')
    try:
        folder,plan=prepared(args.name);summary=final_replay(folder,plan)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as error:
        summary={'complete_eigenpair_set_valid':False,'error':str(error),'production_configuration_admitted':False}
    write(folder/'SUMMARY.json',summary)
    archive=ROOT/('qx_n640_recovery_'+args.name+'_review.zip')
    if archive.exists():raise ValueError('Review archive exists; preserve or rename it first')
    with tempfile.TemporaryDirectory(dir=str(ROOT),prefix='qx_n640_full_export_') as tmp:
        path=Path(tmp)/'review.zip'
        with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(folder.rglob('*')):
                if p.is_file():z.write(p,'campaign/'+str(p.relative_to(folder)))
            for name in package_pins():z.write(HERE/name,'workflow/'+name)
            for p in sorted(BUILD.iterdir()):
                if p.is_file() and p.suffix in ('.c','.h','.txt','.json'):z.write(p,'build_provenance/'+p.name)
        with zipfile.ZipFile(path) as z:
            if z.testzip():raise ValueError('Review ZIP integrity check failed')
        with archive.open('xb') as target,path.open('rb') as source:shutil.copyfileobj(source,target)
    print('PASS full n640 recovery' if summary['complete_eigenpair_set_valid'] else 'INCOMPLETE recovery: '+summary['error'])
    print('REVIEW_ARCHIVE '+str(archive))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check-inputs','verify-pilot','record-build','prepare','submit','submit-final','run-slot','run-final','collect'])
    parser.add_argument('--name',default='full192');parser.add_argument('--concurrent',type=int,default=4)
    parser.add_argument('--slot',type=int,default=-1);args=parser.parse_args()
    if args.action=='check-inputs':check_inputs();print('PASS source pins and selected14 replay')
    elif args.action=='verify-pilot':pilot_replay();print('PASS 14 saved times and 1400 audited modes')
    elif args.action=='record-build':record_build();print('PASS isolated full recovery build')
    else:globals()[args.action.replace('-','_')](args)


if __name__=='__main__':main()
