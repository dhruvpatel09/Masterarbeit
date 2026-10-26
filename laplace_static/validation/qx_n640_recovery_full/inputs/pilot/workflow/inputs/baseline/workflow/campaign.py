#!/usr/bin/env python3
"""Build provenance, run and collect the selected-slice n640 diagnostic."""
import argparse
import hashlib
import json
import math
import re
import struct
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BUILD=HERE/'build'
CODES={0:'converged',1:'allocation_failure',2:'residual_nonconvergence',
       3:'DSYEVR_failure',4:'DSTEMR_failure',5:'upper_cutoff_failure'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read_json(path):return json.loads(Path(path).read_text())
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def package_pins():
    return {str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.rglob('*'))
            if p.is_file() and p.suffix not in ('.out','.err','.log','.o','.zip')
            and not any(x in ('build','campaigns','__pycache__') for x in p.relative_to(HERE).parts)}


def check_inputs():
    spec=read_json(HERE/'INPUTS.json')
    if spec['parameters']['res']!=1e-9:
        raise ValueError('The diagnostic retains the production 1e-9 residual requirement')
    for name,digest in spec['preserved_inputs'].items():
        if sha(HERE/name)!=digest:raise ValueError('Preserved input changed: '+name)
    for name,digest in spec['repository_references'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Production source/library changed: '+name)
    return spec


def config_header(spec):
    p=spec['parameters']
    names={'NV':'Nv','RES':'res','RED_TOL':'redTol','DEG':'deg','UPPER':'upper',
        'REORTHO':'tolReortho','MAX_RESTARTS':'maxRestarts','LANCZOS_FACTOR':'nLanczosFactor',
        'RED_FREQ_FACTOR':'redFreqFactor','NAPE':'laplace_nAPE','ALPHA_APE':'laplace_alphaAPE'}
    return ''.join('#define QX_%s %s\n'%(macro,repr(p[key])) for macro,key in names.items())


def record_build():
    check_inputs()
    names=['test_n640','test_n640_reference','test_n640_reference.o','test_n640.o','qcd_lanczos_trace.o','qcd_lanczos_trace.c',
           'diagnostic_config.h','compiler.txt','modules.txt','linked_libraries.txt','link.map']
    report={'files':{str(BUILD/name):sha(BUILD/name) for name in names},
            'package_files':package_pins(),'repository_references':read_json(HERE/'INPUTS.json')['repository_references']}
    write_json(BUILD/'MANIFEST.json',report)


def verify_build():
    report=read_json(BUILD/'MANIFEST.json')
    for path,digest in report['files'].items():
        if sha(path)!=digest:raise ValueError('Diagnostic build changed: '+path)
    if report['package_files']!=package_pins():raise ValueError('Diagnostic package changed after build')
    check_inputs()
    return report


def directory(name):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',name):raise ValueError('Invalid campaign name')
    return HERE/'campaigns'/name


def prepared(name):
    folder=directory(name)
    plan=read_json(folder/'PREPARED.json')
    if sha(folder/'PREPARED.json')!=(folder/'prepared.sha256').read_text().strip():
        raise ValueError('Prepared plan changed')
    if plan['build_manifest_sha256']!=sha(BUILD/'MANIFEST.json'):raise ValueError('Prepared build changed')
    verify_build()
    return folder,plan


def prepare(args):
    spec=check_inputs();verify_build()
    folder=directory(args.name)
    if folder.exists():raise ValueError('Campaign already exists; use a new name')
    cases=[]
    for case in spec['cases']:
        path=Path(spec['gauge_directory'])/(spec['gauge_prefix']+str(case['configuration']))
        with path.open('rb') as stream:header=stream.read(16)
        if len(header)!=16 or list(struct.unpack('<4i',header))!=spec['geometry']['L']:
            raise ValueError('Unexpected gauge header: '+str(path))
        metadata=path.stat()
        cases.append(dict(case,gauge_path=str(path),gauge_bytes=metadata.st_size,gauge_mtime_ns=metadata.st_mtime_ns,
                          gauge_header_hex=header.hex()))
    folder.mkdir(parents=True)
    plan={'name':args.name,'cases':cases,'geometry':spec['geometry'],'parameters':spec['parameters'],
        'build_manifest_sha256':sha(BUILD/'MANIFEST.json'),'package_files':package_pins(),
        'scope':spec['scope'],'control_gate':'Both n650 controls must pass selected-slice audits and agree in eigenvalues before the n640 target is launched.',
        'control_eigenvalue_absolute_tolerance':spec['control_eigenvalue_absolute_tolerance']}
    write_json(folder/'PREPARED.json',plan)
    (folder/'prepared.sha256').write_text(sha(folder/'PREPARED.json')+'\n')
    print('Prepared original-library and traced n650 controls, then n640 target at global t=1; 3 launches, 64 ranks, 8 nodes')


def submit(args):
    folder,plan=prepared(args.name)
    if (folder/'SUBMISSION.json').exists():raise ValueError('Already submitted; collect or prepare a new campaign')
    result=subprocess.run(['sbatch','--parsable',str(HERE/'submit_micro.sh'),args.name],cwd=str(ROOT),
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    if result.returncode:raise RuntimeError(result.stderr.strip())
    job=result.stdout.strip().split(';')[0]
    if not job.isdigit():raise ValueError('Unexpected sbatch response: '+result.stdout)
    write_json(folder/'SUBMISSION.json',{'job_id':job,'prepared_sha256':sha(folder/'PREPARED.json'),
                                     'submitted_epoch':time.time()})
    print('SUBMITTED '+job+' n650 control and n640 t1 diagnostic')


def parse_log(text,case,parameters):
    if parameters['res']!=1e-9:raise ValueError('Production residual target must remain 1e-9')
    records={}
    for line in text.splitlines():
        if line.startswith('QX_'):
            parts=line.split();kind=parts[0]
            values=dict(part.split('=',1) for part in parts[1:] if '=' in part)
            records.setdefault(kind,[]).append(values)
    def one(kind):
        values=records.get(kind,[])
        if len(values)!=1:raise ValueError('Expected exactly one '+kind)
        return values[0]
    start,end,returned=one('QX_CASE'),one('QX_COMPLETE'),one('QX_RETURN')
    for row in (start,end):
        if (int(row['configuration']),int(row['time']))!=(case['configuration'],case['time']):
            raise ValueError('Wrong configuration/time')
    if (int(start['Nv']),int(start['ranks']),int(start['spatial_ranks']),start['P'])!=(100,64,64,'1,4,4,4'):
        raise ValueError('Wrong diagnostic geometry')
    for key in ('res','redTol','upper','tolReortho'):
        if float(start[key])!=parameters[key]:raise ValueError('Solver parameter changed: '+key)
    for key in ('deg','maxRestarts','nLanczosFactor','redFreqFactor'):
        if int(start[key])!=parameters[key]:raise ValueError('Solver parameter changed: '+key)
    if int(start['nAPE'])!=20 or float(start['alphaAPE'])!=.5 or int(start['mpi_thread_min'])<1:
        raise ValueError('APE/thread setup changed')
    code=int(returned['code'])
    if code not in CODES or any(int(returned[k])!=code for k in ('code_min','code_max')) or int(end['solver_code'])!=code:
        raise ValueError('Inconsistent solver return codes')
    gauge=one('QX_GAUGE');ape=one('QX_APE')
    if not math.isfinite(float(gauge['plaquette'])) or int(ape['steps'])!=20 or float(ape['alpha'])!=.5:
        raise ValueError('Invalid gauge/APE record')
    links=records.get('QX_LINKS',[])
    if len(links)!=2 or {r['stage'] for r in links}!={'raw','ape'}:
        raise ValueError('Missing spatial link checks')
    for r in links:
        if int(r['count'])!=3*64**3 or int(r['nonfinite_links'])!=0:
            raise ValueError('Invalid spatial link inventory')
    if start.get('variant')!=case['variant']:raise ValueError('Wrong executable variant')
    if case['variant']=='traced':
        begins=records.get('QX_CYCLE_BEGIN',[]);ends=records.get('QX_CYCLE_END',[])
        expected_m=parameters['Nv']*parameters['nLanczosFactor']
        if (not begins and code!=1) or [int(r['cycle']) for r in begins]!=list(range(1,len(begins)+1)) or len(begins)>parameters['maxRestarts']:
            raise ValueError('Invalid solver-cycle inventory')
        if any(int(r['m'])!=expected_m for r in begins):raise ValueError('Basis size changed')
        if code in (0,2,5) and len(ends)!=len(begins):raise ValueError('Incomplete completed-cycle inventory')
        if code in (3,4) and len(ends)!=len(begins)-1:raise ValueError('Invalid reduced-solver failure history')
        if code==1 and (begins or ends):raise ValueError('Unexpected history before allocation failure')
        for index,row in enumerate(ends,1):
            if int(row['cycle'])!=index or int(row['restarts'])!=index-1 or int(row['nRed'])!=index or int(row['maxRed'])!=expected_m:
                raise ValueError('Inconsistent cycle counters')
            ritz=[r for r in records.get('QX_RITZ',[]) if int(r['cycle'])==index]
            count=int(row['checked'])
            if not 1<=count<=100 or [int(r['index']) for r in ritz]!=list(range(99,99-count,-1)):
                raise ValueError('Residual trace does not match checked-mode count')
            residuals=[float(r['relative_residual']) for r in ritz]
            traced_max=max([0.]+[x for x in residuals if not math.isnan(x)])
            if traced_max!=float(row['max_checked_res']):raise ValueError('Checked residual maximum mismatch')
        if len(records.get('QX_RITZ',[]))!=sum(int(r['checked']) for r in ends):
            raise ValueError('Unassigned or duplicate residual trace')
        reduced=records.get('QX_REDUCED',[])
        if len(reduced)!=len(begins):raise ValueError('Missing reduced solver status')
        for i,row in enumerate(reduced,1):
            expected_method='DSTEMR' if i==1 else 'DSYEVR'
            if int(row['cycle'])!=i or row['method']!=expected_method or int(row['m'])!=expected_m:
                raise ValueError('Wrong reduced-solver trace')
            if i<=len(ends) and int(row['info'])!=0:raise ValueError('Completed cycle has failed reduced solver')
        if code in (3,4) and int(reduced[-1]['info'])==0:raise ValueError('Missing reduced-solver error code')
        if code in (0,2,5):
            if (int(returned['nRed'])!=len(ends) or int(returned['restarts'])!=len(ends)-1 or
                int(returned['maxRed'])!=expected_m or float(returned['maxRes'])!=float(ends[-1]['max_checked_res'])):
                raise ValueError('Final status disagrees with cycle history')
            if int(ends[-1]['checked'])!=100:raise ValueError('Final cycle did not check all requested modes')
            if code==2 and (len(ends)!=parameters['maxRestarts'] or float(returned['maxRes'])<=parameters['res']):
                raise ValueError('Nonconvergence return without exhausted restart budget or failing residual')
            if code in (0,5) and float(returned['maxRes'])>parameters['res']:
                raise ValueError('Convergence/cutoff return with failing residual')
            if code==5 and float(ends[-1]['lambda_Nv'])<parameters['upper']:
                raise ValueError('Cutoff return with eigenvalue below cutoff')
            if code==0 and float(ends[-1]['lambda_Nv'])>=parameters['upper']:
                raise ValueError('Successful return with eigenvalue above cutoff')
    elif any(records.get(k) for k in ('QX_CYCLE_BEGIN','QX_CYCLE_END','QX_RITZ','QX_REDUCED')):
        raise ValueError('Original library unexpectedly contains diagnostic instrumentation')
    else:ends=[]
    audit=one('QX_AUDIT');valid=False
    if code in (0,2,5):
        modes=records.get('QX_AUDIT_MODE',[])
        if len(modes)!=100 or [int(r['index']) for r in modes]!=list(range(100)) or int(audit['modes'])!=100:
            raise ValueError('Incomplete independent residual audit')
        finite=all(int(r['finite_positive'])==1 and all(math.isfinite(float(r[k])) for k in ('lambda','norm','relative_residual'))
                   and float(r['lambda'])>0 and float(r['norm'])>0 for r in modes)
        if finite:
            residual=max(float(r['relative_residual']) for r in modes)
            norm=max(abs(float(r['norm'])-1) for r in modes)
            if residual!=float(audit['max_relative_residual']) or abs(norm-float(audit['max_norm_defect']))>1e-15:
                raise ValueError('Independent audit aggregate mismatch')
        valid=(code==0 and finite and int(audit['invalid'])==0 and
               float(audit['max_relative_residual'])<=parameters['res'] and
               float(audit['max_norm_defect'])<=1e-10 and float(audit['max_gram_offdiag'])<=1e-8)
        if int(audit['valid_selected_slice'])!=int(valid):raise ValueError('Selected-slice validity flag mismatch')
    elif audit.get('skipped')!='1':raise ValueError('Unsafe audit after early solver failure')
    return {'diagnostic_completed':True,'configuration':case['configuration'],'time':case['time'],
        'solver_code':code,'solver_result':CODES[code],'selected_slice_valid':valid,
        'production_configuration_admitted':False,'cycles':len(ends) if case['variant']=='traced' else int(returned['nRed']),'cycle_history':ends,
        'variant':case['variant'],'return_status':returned,'spatial_link_checks':links,'independent_audit':audit,
        'audited_modes':records.get('QX_AUDIT_MODE',[])}


def compare_controls(reports,tolerance):
    if len(reports)!=2 or [r['variant'] for r in reports]!=['library','traced'] or any(
        r['configuration']!=650 or r['time']!=1 or not r['selected_slice_valid'] for r in reports):
        raise ValueError('Both original-library and traced n650 controls must validate')
    differences=[abs(float(a['lambda'])-float(b['lambda']))
                 for a,b in zip(reports[0]['audited_modes'],reports[1]['audited_modes'])]
    if len(differences)!=100 or max(differences)>tolerance:
        raise ValueError('Original-library and traced control eigenvalues disagree; target not admitted')
    return {'success':True,'configuration':650,'time':1,'maximum_absolute_eigenvalue_difference':max(differences),
            'absolute_tolerance':tolerance,'scope':'One selected-slice control comparison; not a proof of source identity.'}


def run(args):
    folder,plan=prepared(args.name)
    output=folder/'outputs';output.mkdir()
    launches=[];reports=[]
    for case in plan['cases']:
        gauge=Path(case['gauge_path']);stat=gauge.stat()
        if stat.st_size!=case['gauge_bytes'] or stat.st_mtime_ns!=case['gauge_mtime_ns']:
            raise ValueError('Gauge file metadata changed after preparation')
        label='n%d_t%d_%s'%(case['configuration'],case['time'],case['variant'])
        command=['mpiexec','-n','64','-ppn','8',str(BUILD/('test_n640_reference' if case['variant']=='library' else 'test_n640')),str(gauge),str(case['configuration']),str(case['time'])]
        start=time.time()
        with (output/(label+'.out')).open('w') as out,(output/(label+'.err')).open('w') as err:
            process=subprocess.run(command,cwd=str(ROOT),stdout=out,stderr=err)
        launch=dict(case,command=command,exit_code=process.returncode,seconds=time.time()-start,
                    stdout_sha256=sha(output/(label+'.out')),stderr_sha256=sha(output/(label+'.err')))
        launches.append(launch);write_json(folder/'LAUNCHES.json',launches)
        if process.returncode:raise ValueError('MPI diagnostic did not complete: '+label)
        report=parse_log((output/(label+'.out')).read_text(),case,plan['parameters'])
        write_json(output/(label+'.json'),report)
        reports.append(report)
        if case['role']!='target' and not report['selected_slice_valid']:
            raise ValueError('Control failed; target was not launched')
        if case['role']=='trace_control':
            write_json(folder/'CONTROL_COMPARISON.json',compare_controls(reports[:2],plan['control_eigenvalue_absolute_tolerance']))
    print('Diagnostic launches completed; collect to replay and export')


def collect(args):
    folder=directory(args.name)
    if not folder.is_dir():raise ValueError('Campaign directory does not exist')
    try:
        folder,plan=prepared(args.name)
        launches=read_json(folder/'LAUNCHES.json')
        if len(launches)!=len(plan['cases']):raise ValueError('Incomplete launch inventory')
        reports=[]
        for case,launch in zip(plan['cases'],launches):
            label='n%d_t%d_%s'%(case['configuration'],case['time'],case['variant'])
            if launch['configuration']!=case['configuration'] or launch['time']!=case['time'] or launch['exit_code']!=0:
                raise ValueError('Failed or wrong MPI launch')
            for suffix,key in (('out','stdout_sha256'),('err','stderr_sha256')):
                if sha(folder/'outputs'/(label+'.'+suffix))!=launch[key]:raise ValueError('Changed diagnostic output')
            reports.append(parse_log((folder/'outputs'/(label+'.out')).read_text(),case,plan['parameters']))
        control_comparison=compare_controls(reports[:2],plan['control_eigenvalue_absolute_tolerance'])
        summary={'diagnostic_completed':True,'cases':reports,'control_comparison':control_comparison,
            'production_configuration_admitted':False,
            'scope':'Selected global time slice only; full n640 generation and independent saved-mode audits remain required.'}
    except (OSError,ValueError,KeyError,TypeError,IndexError) as error:
        summary={'diagnostic_completed':False,'error':str(error),'production_configuration_admitted':False,
                 'scope':'Incomplete or invalid diagnostic evidence; logs are exported for inspection.'}
    write_json(folder/'SUMMARY.json',summary)
    archive=ROOT/('qx_n640_diagnostic_'+args.name+'_review.zip')
    if archive.exists():raise ValueError('Review archive already exists; preserve or rename it before exporting again')
    with tempfile.TemporaryDirectory(dir=str(ROOT),prefix='qx_n640_export_') as temporary:
        path=Path(temporary)/'review.zip'
        with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in sorted(folder.rglob('*')):
                if p.is_file():z.write(p,str(p.relative_to(HERE)))
            for name in ('MANIFEST.json','modules.txt','compiler.txt','linked_libraries.txt','link.map','qcd_lanczos_trace.c','diagnostic_config.h'):
                if (BUILD/name).is_file():z.write(BUILD/name,'build_provenance/'+name)
            for name in package_pins():z.write(HERE/name,'workflow/'+name)
        with zipfile.ZipFile(path) as z:
            if z.testzip() is not None:raise ValueError('Review ZIP failed integrity check')
        with archive.open('xb') as target,path.open('rb') as source:
            for chunk in iter(lambda:source.read(1024*1024),b''):target.write(chunk)
    if summary['diagnostic_completed']:
        print('PASS selected-slice diagnostic replay; n650 controls valid; n640 solver result: '+reports[-1]['solver_result'])
    else:print('INCOMPLETE diagnostic: '+summary['error'])
    print('REVIEW_ARCHIVE '+str(archive))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check-inputs','config-header','record-build','prepare','submit','run','collect'])
    parser.add_argument('--name',default='baseline')
    args=parser.parse_args()
    if args.action=='check-inputs':check_inputs();print('PASS diagnostic source/library pins')
    elif args.action=='config-header':print(config_header(read_json(HERE/'INPUTS.json')),end='')
    elif args.action=='record-build':record_build();print('PASS isolated diagnostic build manifest')
    else:globals()[args.action](args)


if __name__=='__main__':main()
