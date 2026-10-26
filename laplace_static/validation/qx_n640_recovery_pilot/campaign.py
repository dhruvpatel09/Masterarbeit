#!/usr/bin/env python3
"""Prepare, run and export an isolated n640 selected-time saved-mode pilot."""
import argparse
import hashlib
import importlib.util
import json
import math
import re
import shutil
import struct
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BUILD = HERE/'build'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())
def write(path, obj):Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def package_pins():
    return {str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.rglob('*'))
            if p.is_file() and '__pycache__' not in p.relative_to(HERE).parts
            and (p.relative_to(HERE).parts[0]=='inputs' or
                 (not any(x in ('build','campaigns') for x in p.relative_to(HERE).parts)
                  and p.suffix not in ('.out','.err','.log','.o','.zip')))}


def baseline():
    base = HERE/'inputs/baseline'
    spec = read(base/'workflow/INPUTS.json')
    module_spec = importlib.util.spec_from_file_location('qx_baseline_replay',str(base/'workflow/campaign.py'))
    module = importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(module)
    reports = []
    launches = read(base/'campaigns/baseline/LAUNCHES.json')
    if len(launches)!=3:raise ValueError('Incomplete baseline')
    for case, launch in zip(spec['cases'], launches):
        label = 'n%d_t%d_%s'%(case['configuration'],case['time'],case['variant'])
        for ext, key in (('out','stdout_sha256'),('err','stderr_sha256')):
            if sha(base/'campaigns/baseline/outputs'/(label+'.'+ext))!=launch[key]:
                raise ValueError('Baseline output checksum changed')
        if launch['exit_code']!=0:raise ValueError('Baseline MPI launch failed')
        report = module.parse_log((base/'campaigns/baseline/outputs'/(label+'.out')).read_text(),case,spec['parameters'])
        if not report['selected_slice_valid']:raise ValueError('Baseline selected slice is invalid')
        if report!=read(base/'campaigns/baseline/outputs'/(label+'.json')):
            raise ValueError('Baseline report mismatch')
        reports.append(report)
    comparison = module.compare_controls(reports[:2],spec['control_eigenvalue_absolute_tolerance'])
    summary = read(base/'campaigns/baseline/SUMMARY.json')
    if (summary['cases']!=reports or summary['control_comparison']!=comparison or
        not summary['diagnostic_completed'] or summary['production_configuration_admitted']):
        raise ValueError('Baseline summary mismatch')
    return reports[-1]


def check_inputs():
    spec = read(HERE/'INPUTS.json')
    if (spec['configuration']!=640 or spec['parameters']['res']!=1e-9 or
        spec['selected_times']!=[0,1,32,64,95,96,97,98,99,100,101,128,160,191]):
        raise ValueError('Recovery scope changed')
    for name, digest in spec['preserved_inputs'].items():
        if sha(HERE/name)!=digest:raise ValueError('Changed preserved input: '+name)
    for name, digest in spec['repository_references'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Changed production dependency: '+name)
    baseline()
    return spec


def record_build():
    spec = check_inputs()
    names = ['mental_recovery','check_saved','mental_recovery.c','check_saved.c','selection.h',
             'mental_recovery.o','check_saved.o','DistRedelemsHdf5Writer.o','DistEigvecsHdf5Writer.o',
             'DistEigvecsHdf5Reader.o','DistEigvalsHdf5Writer.o','linked_generator.txt',
             'linked_checker.txt','modules.txt','compiler.txt']
    write(BUILD/'MANIFEST.json',{'files':{n:sha(BUILD/n) for n in names},
         'package_files':package_pins(),'repository_references':spec['repository_references']})


def verify_build():
    check_inputs()
    manifest = read(BUILD/'MANIFEST.json')
    if manifest['package_files']!=package_pins():raise ValueError('Package changed after build')
    for name, digest in manifest['files'].items():
        if sha(BUILD/name)!=digest:raise ValueError('Build changed: '+name)
    return manifest


def directory(name):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def input_text(run):
    return '''[Run name]
name qX_p2
output qX_p2_n640_selected_recovery
[Directories]
log_dir {run}/log
cnfg_dir /hppfs/scratch/00/di98qes3/qX
dat_dir {run}/dat
evec_dir {run}/eigenvectors
eval_dir {run}/eigenvalues
[Lattice]
glob_size 192 64 64 64
proc_grid 1 4 4 4
[Configurations]
first 640
last 640
step 10
[Measurements]
nset 1
Nv 100
compute_ev 1
write_ev 1
laplace_nAPE 20
laplace_alphaAPE 0.5
[Lanczos]
res 1e-9
redTol 1e-12
deg 10
upper 0.15
tolReortho 1e-10
maxRestarts 30
nLanczosFactor 7
redFreqFactor 10
[set 0]
p2max 0
nder 0
'''.format(run=run)


def gauge_record(path):
    path = Path(path);st = path.stat()
    with path.open('rb') as f:head = f.read(16)
    if len(head)!=16 or struct.unpack('<4i',head)!=(192,64,64,64):raise ValueError('Gauge header mismatch')
    if st.st_size<192*64**3*8:raise ValueError('Gauge file is unexpectedly small')
    return {'path':str(path),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns,'header_hex':head.hex()}


def prepare(args):
    spec = check_inputs();verify_build()
    folder = directory(args.name)
    if folder.exists():raise ValueError('Campaign already exists; preserve it')
    gauge = gauge_record(spec['gauge_path'])
    storage = Path(spec['storage_root']);storage.mkdir(parents=True,exist_ok=True)
    if shutil.disk_usage(str(storage)).free<25*1024**3:raise ValueError('Need at least 25 GiB free for this selected-time pilot')
    stamp = time.strftime('%Y%m%d_%H%M%S',time.gmtime())
    run = storage/(args.name+'_'+stamp)
    run.mkdir()
    for name in ('log','dat','eigenvectors','eigenvalues'):(run/name).mkdir()
    folder.mkdir(parents=True)
    (folder/'mental.in').write_text(input_text(run))
    (run/'SELECTED_TIME_PILOT.json').write_text(json.dumps({'selected_times':spec['selected_times'],
        'complete_configuration':False,'production_configuration_admitted':False,'scope':spec['scope']},indent=2)+'\n')
    plan = {'name':args.name,'selected_times':spec['selected_times'],'configuration':640,'run_directory':str(run),
        'gauge':gauge,'build_manifest_sha256':sha(BUILD/'MANIFEST.json'),'input_sha256':sha(folder/'mental.in'),
        'eigenvalue_file':str(run/'eigenvalues/qX_p2_n640_selected_recovery_n640_eval.h5'),
        'production_configuration_admitted':False,'scope':spec['scope']}
    write(folder/'PREPARED.json',plan)
    (folder/'prepared.sha256').write_text(sha(folder/'PREPARED.json')+'\n')
    print('Prepared 14 selected times, 1400 saved modes; 64 ranks, 8 nodes; generator then independent read-back audit')


def prepared(name):
    folder = directory(name);plan = read(folder/'PREPARED.json')
    if sha(folder/'PREPARED.json')!=(folder/'prepared.sha256').read_text().strip():raise ValueError('Plan changed')
    if sha(folder/'mental.in')!=plan['input_sha256']:raise ValueError('Input changed')
    if sha(BUILD/'MANIFEST.json')!=plan['build_manifest_sha256']:raise ValueError('Prepared build changed')
    verify_build()
    if gauge_record(plan['gauge']['path'])!=plan['gauge']:raise ValueError('Gauge changed')
    return folder,plan


def submit(args):
    folder,plan = prepared(args.name)
    if (folder/'SUBMISSION.json').exists():raise ValueError('Campaign already submitted')
    r = subprocess.run(['sbatch','--parsable',str(HERE/'submit_micro.sh'),args.name],cwd=str(ROOT),
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if r.returncode:raise ValueError(r.stderr)
    job = r.stdout.strip().split(';')[0]
    if not job.isdigit():raise ValueError('Unexpected sbatch output: '+r.stdout)
    write(folder/'SUBMISSION.json',{'job_id':job,'prepared_sha256':sha(folder/'PREPARED.json')})
    print('SUBMITTED '+job+' n640 selected saved-mode pilot')


def records(text,kind):
    return [dict(x.split('=',1) for x in line.split()[1:] if '=' in x)
            for line in text.splitlines() if line.startswith(kind+' ')]


def replay_logs(generator, checker, times):
    starts = records(generator,'RECOVERY_START');solves = records(generator,'RECOVERY_SOLVE')
    written = records(generator,'RECOVERY_WRITTEN')
    if len(starts)!=1 or starts[0]!={'cfg':'640','ranks':'64','P':'1,4,4,4','selected':str(len(times)),'res':'1.00000000000000006e-09'}:
        raise ValueError('Missing or changed generation setup')
    for rows in (solves,written):
        if [int(r['t']) for r in rows]!=times or any(r['cfg']!='640' for r in rows):
            raise ValueError('Incomplete or wrong generated time inventory')
    for r in solves:
        if (int(r['code'])!=0 or not 1<=int(r['nRed'])<=30 or int(r['restarts'])!=int(r['nRed'])-1 or
            not math.isfinite(float(r['maxRes'])) or not 0<=float(r['maxRes'])<=1e-9 or
            not 0<float(r['lambda_Nv'])<0.15):raise ValueError('Solver did not meet original target')
    if any(r['Nv']!='100' for r in written):raise ValueError('Wrong saved mode count')
    begin = records(checker,'QUALITY_START');end = records(checker,'QUALITY_END')
    if len(begin)!=1 or begin[0]['cfg']!='640' or begin[0]['grid']!='1x4x4x4' or begin[0]['Nv']!='100' or begin[0]['APE']!='20,0.5' or int(begin[0]['FUNNELED'])<1:
        raise ValueError('Invalid saved-mode audit setup')
    if len(end)!=1 or end[0]!={'cfg':'640','slices':str(len(times)),'status':'PASS'}:
        raise ValueError('Saved-mode audit did not pass')
    modes = records(checker,'MODE');slices = records(checker,'SLICE')
    if len(modes)!=100*len(times) or [int(r['t']) for r in slices]!=times:raise ValueError('Saved-mode audit inventory is incomplete')
    baseline_modes = baseline()['audited_modes']
    t1difference = 0.
    for t,row in zip(times,slices):
        group = [m for m in modes if int(m['t'])==t]
        if [int(m['i']) for m in group]!=list(range(100)):raise ValueError('Duplicate or missing saved mode')
        for m in group:
            if m['cfg']!='640' or any(not math.isfinite(float(m[k])) for k in ('lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im')):
                raise ValueError('Nonfinite or wrong saved mode')
            if (not 0<float(m['lambda'])<0.15 or abs(float(m['norm'])-1)>1e-9 or
                not 0<=float(m['residual_relative'])<=1e-9 or float(m['residual_l2'])<0):raise ValueError('Saved pair fails residual/norm/cutoff target')
            derived = float(m['residual_l2'])/(abs(float(m['lambda']))*float(m['norm']))
            if not math.isclose(derived,float(m['residual_relative']),rel_tol=1e-12,abs_tol=1e-16):raise ValueError('Residual formula mismatch')
        if any(float(group[i]['lambda'])+1e-12<float(group[i-1]['lambda']) for i in range(1,100)):
            raise ValueError('Unordered saved eigenvalues')
        generated = solves[times.index(t)]
        if abs(float(generated['lambda_Nv'])-float(group[-1]['lambda']))>1e-10:
            raise ValueError('Saved largest eigenvalue disagrees with generator')
        if row['cfg']!='640' or row['Nv']!='100' or row['status']!='PASS':raise ValueError('Slice audit failed')
        limits = {'max_relative_residual':1e-9,'max_norm_defect':1e-9,'gram_diag_defect':1e-9,'gram_offdiag':1e-9,'max_rayleigh_im':1e-10}
        for key,limit in limits.items():
            x = float(row[key])
            if not math.isfinite(x) or not 0<=x<=limit:raise ValueError('Slice audit bound failed: '+key)
        for key,actual in [('max_relative_residual',max(float(m['residual_relative']) for m in group)),
                           ('max_norm_defect',max(abs(float(m['norm'])-1) for m in group))]:
            if not math.isclose(float(row[key]),actual,rel_tol=1e-12,abs_tol=1e-15):raise ValueError('Slice audit aggregate mismatch')
        if (float(row['lambda_min'])!=float(group[0]['lambda']) or
            float(row['lambda_max'])!=float(group[-1]['lambda'])):raise ValueError('Eigenvalue range mismatch')
        if t==1:
            t1difference = max(abs(float(m['lambda'])-float(b['lambda'])) for m,b in zip(group,baseline_modes))
    if t1difference>1e-10:raise ValueError('Saved t1 eigenvalues disagree with diagnostic baseline')
    return {'selected_saved_modes_valid':True,'selected_times':times,'modes_per_time':100,
        'audited_modes':len(modes),'slices':slices,'t1_baseline_max_absolute_eigenvalue_difference':t1difference,
        'complete_configuration':False,'production_configuration_admitted':False,
        'scope':'Selected-time pilot; complete 192-time coverage and full-ensemble admission are separate steps.'}


def run(args):
    folder,plan = prepared(args.name);output = folder/'outputs';output.mkdir()
    commands = [
        ('generate',['mpiexec','-n','64','-ppn','8',str(BUILD/'mental_recovery'),'-i',str(folder/'mental.in')]),
        ('audit',['mpiexec','-n','64','-ppn','8',str(BUILD/'check_saved'),'640',plan['gauge']['path'],
                  str(Path(plan['run_directory'])/'eigenvectors'),plan['eigenvalue_file']])]
    launches = []
    for label,command in commands:
        began = time.time()
        with (output/(label+'.out')).open('w') as out,(output/(label+'.err')).open('w') as err:
            r = subprocess.run(command,cwd=plan['run_directory'],stdout=out,stderr=err)
        launches.append({'label':label,'command':command,'exit_code':r.returncode,'seconds':time.time()-began,
            'stdout_sha256':sha(output/(label+'.out')),'stderr_sha256':sha(output/(label+'.err'))})
        write(folder/'LAUNCHES.json',launches)
        if r.returncode:raise ValueError('MPI launch failed: '+label)
    print('Generation and read-back completed; collect to replay and export')


def file_inventory(plan):
    run = Path(plan['run_directory'])
    paths = [run/'eigenvectors'/('qX_p2n640_evec_t%d.h5'%t) for t in plan['selected_times']]+[Path(plan['eigenvalue_file'])]
    if sorted(p.name for p in (run/'eigenvectors').glob('*.h5'))!=sorted(p.name for p in paths[:-1]):
        raise ValueError('Unexpected eigenvector file inventory')
    inventory = []
    for p in paths:
        st = p.stat()
        if st.st_size<1024:raise ValueError('Empty or truncated HDF5 file')
        inventory.append({'path':str(p),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns})
    return inventory


def collect(args):
    folder = directory(args.name)
    if not folder.is_dir():raise ValueError('Campaign does not exist')
    try:
        folder,plan = prepared(args.name);launches = read(folder/'LAUNCHES.json')
        if [x['label'] for x in launches]!=['generate','audit']:raise ValueError('Incomplete launch inventory')
        for x in launches:
            if x['exit_code']!=0:raise ValueError('MPI launch failed: '+x['label'])
            for ext,key in (('out','stdout_sha256'),('err','stderr_sha256')):
                if sha(folder/'outputs'/(x['label']+'.'+ext))!=x[key]:raise ValueError('Changed MPI output')
        summary = replay_logs((folder/'outputs/generate.out').read_text(),(folder/'outputs/audit.out').read_text(),plan['selected_times'])
        inventory = file_inventory(plan)
        old = folder/'SAVED_FILES.json'
        if old.exists() and read(old)!=inventory:raise ValueError('Saved data changed after first collection')
        write(old,inventory)
    except (OSError,ValueError,KeyError,TypeError,IndexError) as error:
        summary = {'selected_saved_modes_valid':False,'error':str(error),'complete_configuration':False,
                   'production_configuration_admitted':False}
    write(folder/'SUMMARY.json',summary)
    archive = ROOT/('qx_n640_recovery_pilot_'+args.name+'_review.zip')
    if archive.exists():raise ValueError('Review ZIP exists; preserve or rename it first')
    with tempfile.TemporaryDirectory(dir=str(ROOT),prefix='qx_n640_recovery_export_') as tmp:
        zip_path = Path(tmp)/'review.zip'
        with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(folder.rglob('*')):
                if p.is_file():z.write(p,'campaign/'+str(p.relative_to(folder)))
            for name in package_pins():z.write(HERE/name,'workflow/'+name)
            for p in sorted(BUILD.iterdir()):
                if p.is_file() and (p.suffix in ('.c','.h','.json','.txt')):z.write(p,'build_provenance/'+p.name)
            plan_file = folder/'PREPARED.json'
            if plan_file.exists():
                run = Path(read(plan_file)['run_directory'])
                for p in sorted((run/'log').glob('*')):
                    if p.is_file() and p.stat().st_size<10*1024**2:z.write(p,'run_log/'+p.name)
        with zipfile.ZipFile(zip_path) as z:
            if z.testzip():raise ValueError('Review ZIP integrity check failed')
        with archive.open('xb') as target,zip_path.open('rb') as src:shutil.copyfileobj(src,target)
    print('PASS selected saved-mode pilot' if summary['selected_saved_modes_valid'] else 'INCOMPLETE pilot: '+summary['error'])
    print('REVIEW_ARCHIVE '+str(archive))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check-inputs','verify-baseline','record-build','prepare','submit','run','collect'])
    parser.add_argument('--name',default='selected14')
    args = parser.parse_args()
    if args.action=='check-inputs':check_inputs();print('PASS recovery dependencies and baseline replay')
    elif args.action=='verify-baseline':baseline();print('PASS baseline: n650 controls exact; n640 t1 valid; no full-configuration admission')
    elif args.action=='record-build':record_build();print('PASS isolated writer/checker build manifest')
    else:globals()[args.action](args)


if __name__=='__main__':main()
