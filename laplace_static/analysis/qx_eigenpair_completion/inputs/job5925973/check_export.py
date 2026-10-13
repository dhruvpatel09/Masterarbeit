#!/usr/bin/env python3
import csv,json,math,sys
from pathlib import Path
CONFIGS=[800]
TIMES=[0,32,64,95,96,97,98,99,100,101,128,160,191]
MODE_FIELDS=['lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im']
SLICE_FIELDS=['max_relative_residual','max_norm_defect','gram_diag_defect','gram_offdiag','max_rayleigh_im','lambda_min','lambda_max']

def parse(path,cfg):
    modes={};slices={};end=None;started=0
    for line in path.read_text().splitlines():
        if not line.startswith(('MODE ','SLICE ','QUALITY_START ','QUALITY_END ')):continue
        tag,*bits=line.split();d={}
        for bit in bits:
            k,v=bit.split('=',1)
            if k in d:raise ValueError('Duplicate field: '+k)
            d[k]=v
        if int(d['cfg'])!=cfg:raise ValueError('Wrong configuration')
        if tag=='QUALITY_START':started+=1;continue
        if tag=='QUALITY_END':
            if end is not None:raise ValueError('Duplicate end marker')
            end=d;continue
        t=int(d['t'])
        if t not in TIMES:raise ValueError('Unexpected time')
        if tag=='MODE':
            i=int(d['i']);key=(t,i)
            if not 0<=i<100 or key in modes:raise ValueError('Invalid/duplicate mode')
            row={'configuration':cfg,'t':t,'mode':i}
            row.update({k:float(d[k]) for k in MODE_FIELDS})
            if not all(math.isfinite(row[k]) for k in MODE_FIELDS):raise ValueError('Nonfinite mode')
            if row['lambda']<=0 or row['norm']<=0 or row['residual_l2']<0:raise ValueError('Invalid mode numbers')
            q=row['residual_l2']/(abs(row['lambda'])*row['norm'])
            if not math.isclose(q,row['residual_relative'],rel_tol=5e-13,abs_tol=1e-30):raise ValueError('Residual definition mismatch')
            modes[key]=row
        else:
            if t in slices or int(d['Nv'])!=100:raise ValueError('Duplicate slice/wrong Nv')
            row={'configuration':cfg,'t':t,'status':d['status']}
            row.update({k:float(d[k]) for k in SLICE_FIELDS})
            if not all(math.isfinite(row[k]) for k in SLICE_FIELDS):raise ValueError('Nonfinite slice')
            if any(row[k]<0 for k in SLICE_FIELDS):raise ValueError('Negative slice diagnostic')
            slices[t]=row
    if started!=1 or end is None or int(end['slices'])!=len(TIMES):raise ValueError('Missing start/end coverage')
    if set(modes)!={(t,i) for t in TIMES for i in range(100)} or set(slices)!=set(TIMES):raise ValueError('Incomplete mode/slice inventory')
    success=True
    for t,row in slices.items():
        entries=[modes[t,i] for i in range(100)]
        maxima={'max_relative_residual':max(r['residual_relative'] for r in entries),
            'max_norm_defect':max(abs(r['norm']-1) for r in entries),
            'max_rayleigh_im':max(abs(r['rayleigh_im']) for r in entries),
            'lambda_min':entries[0]['lambda'],'lambda_max':entries[-1]['lambda']}
        for k,v in maxima.items():
            if not math.isclose(v,row[k],rel_tol=5e-13,abs_tol=1e-24):raise ValueError('Slice reduction mismatch: '+k)
        if any(entries[i]['lambda']+1e-12<entries[i-1]['lambda'] for i in range(1,100)):raise ValueError('Unordered selected spectrum')
        passed=row['max_relative_residual']<=1.1e-9 and row['max_norm_defect']<=1e-9 and row['gram_diag_defect']<=1e-9 and row['gram_offdiag']<=1e-9 and row['max_rayleigh_im']<=1e-10
        if row['status']!=('PASS' if passed else 'FAIL'):raise ValueError('Wrong slice status')
        success=success and passed
    if end['status']!=('PASS' if success else 'FAIL'):raise ValueError('Wrong end status')
    return success,[modes[k] for k in sorted(modes)],[slices[k] for k in sorted(slices)]

def export(result):
    launches=json.loads((result/'LAUNCHES.json').read_text());all_modes=[];all_slices=[];runs=[]
    if [r['configuration'] for r in launches]!=list(CONFIGS):raise ValueError('Unexpected launch coverage')
    for launch in launches:
        cfg=launch['configuration'];entry={'configuration':cfg,'mpi_exit_code':launch['mpi_exit_code'],'success':False}
        try:
            passed,modes,slices=parse(result/('n%d.out'%cfg),cfg)
            all_modes.extend(modes);all_slices.extend(slices)
            entry['success']=passed and launch['mpi_exit_code']==0
            entry['maximum_relative_residual']=max(r['residual_relative'] for r in modes)
            entry['modes_above_requested_1e_minus9']=sum(r['residual_relative']>1e-9 for r in modes)
        except (ValueError,KeyError,OSError) as exc:entry['error']=str(exc)
        runs.append(entry)
    for filename,rows,fields in [('modes.csv',all_modes,['configuration','t','mode']+MODE_FIELDS),('slices.csv',all_slices,['configuration','t','status']+SLICE_FIELDS)]:
        with (result/filename).open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    success=len(runs)==len(CONFIGS) and all(r['success'] for r in runs)
    summary={'success':success,'runs':runs,'mode_records':len(all_modes),'slice_records':len(all_slices),
        'scope':'Saved-pair residual/normalization/all-pair Gram check on 13 selected slices in each audited configuration, using existing qcd Laplacian. Eigenvalue metadata, coverage and spectra checked for all times. Not full eigenvector coverage or independent operator implementation.',
        'residual_tolerance':1.1e-9,'requested_solver_res':1e-9,'normalization_and_gram_tolerance':1e-9,'rayleigh_im_tolerance':1e-10}
    with (result/'SUMMARY.json').open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
    print(('PASS' if success else 'FAIL'),'saved-eigenmode audit',json.dumps(summary));return success
if __name__=='__main__':sys.exit(0 if export(Path(sys.argv[1])) else 1)
