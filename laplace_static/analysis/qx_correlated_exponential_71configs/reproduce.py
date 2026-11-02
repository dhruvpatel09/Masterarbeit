#!/usr/bin/env python3
"""Replay numerical exports from archived solutions without rerunning optimizers."""
import argparse
import csv
import json
import math
import numpy as np
import fit_models as fm
from run_analysis import HERE,RESULTS,CONFIGS,CUTOFFS,RADII,TWO_R,load_data,parameters,make_row,key
from derived_results import load_rows,window_for,derive


def same(actual,expected,label,rtol=3e-5,atol=1e-8):
    if isinstance(expected,(bool,str)) or expected is None:
        if actual!=expected:raise ValueError('Mismatch '+label+': '+str((actual,expected)))
    elif not math.isfinite(float(actual)) or not np.isclose(float(actual),float(expected),rtol=rtol,atol=atol):
        raise ValueError('Numerical mismatch '+label+': '+str((actual,expected)))


def check_csv(path,expected):
    with path.open(newline='') as f:actual=list(csv.DictReader(f))
    if len(actual)!=len(expected):raise ValueError('Row count '+path.name)
    for i,(a,e) in enumerate(zip(actual,expected)):
        if set(a)!=set(e):raise ValueError('Columns '+path.name)
        for k,v in e.items():
            value=None if v is None and a[k]=='' else a[k]=='True' if isinstance(v,bool) else a[k] if isinstance(v,str) else float(a[k])
            same(value,v,path.name+':'+str(i)+':'+k)


def verify():
    a=load_data();rows=load_rows();solutions=json.loads((RESULTS/'solutions.json').read_text())
    expected=set(key(n,R,lo,hi,s) for n in CUTOFFS for R in RADII for lo in (4,8,12,16) for hi in (24,32)
                 for s in ((1,2) if R in TWO_R else (1,)))
    if set(solutions)!=expected or {r['fit_id'] for r in rows}!=expected or len(rows)!=544:raise ValueError('Fit inventory')
    bootstrap_indices=json.loads((RESULTS/'bootstrap_indices.json').read_text())
    if bootstrap_indices['configuration_labels']!=CONFIGS or bootstrap_indices['seed']!=20261102:raise ValueError('Bootstrap labels')
    chosen=np.array(bootstrap_indices['indices'],dtype=int)
    np.testing.assert_array_equal(chosen,np.random.RandomState(20261102).randint(0,71,size=(399,71)))
    summary=json.loads((RESULTS/'SUMMARY.json').read_text());bootstrap=[];model_comparison=[]
    replayed_starts=0;replayed_jackknife=0;replayed_bootstrap=0
    lookup={(r['Nv'],r['R'],r['tau_min'],r['tau_max'],r['states']):r for r in rows}
    for r in rows:
        w=window_for(r,a);solution=solutions[r['fit_id']];candidates=solution['candidates']
        if len(candidates)!=(5 if r['states']==1 else 13):raise ValueError('Starting-point inventory')
        for c in candidates:
            if c['p'] is not None:
                residual=w.residual(np.array(c['p']),r['states'])
                same(float(residual@residual),c['q2'],r['fit_id']+':candidate objective',rtol=1e-7)
                replayed_starts+=1
        valid=[c for c in candidates if c['success'] and c['p'] is not None and math.isfinite(c['q2'])]
        best=min(valid,key=lambda c:c['q2']);best=dict(best,p=np.array(best['p']))
        same(best['start'],solution['best_start'],r['fit_id']+':best start')
        expected_row,d=make_row(r['Nv'],r['R'],r['tau_min'],r['tau_max'],r['states'],w,best,candidates)
        for field,value in expected_row.items():
            if field.startswith('jackknife_') or field=='E0_jackknife_se_formal':continue
            same(r[field],value,r['fit_id']+':'+field,rtol=3e-5 if r['identified'] else 5e-4)
        if d['local'] is None:
            same(solution['local_covariance'],None,r['fit_id']+':local covariance')
        else:np.testing.assert_allclose(solution['local_covariance'],d['local']['covariance'],rtol=5e-4,atol=1e-8)
        requested=r['Nv'] in (40,100) and r['R'] in (4,12,24,32) and r['tau_min'] in (8,12) and r['tau_max']==32
        same(r['jackknife_requested'],requested,r['fit_id']+':audit request')
        performed=requested and r['identified'];same(r['jackknife_performed'],performed,r['fit_id']+':audit performed')
        jk=solution['jackknife'];converged=0;resolved=0;energies=[]
        if len(jk)!=(71 if performed else 0):raise ValueError('Jackknife inventory '+r['fit_id'])
        for i,record in enumerate(jk):
            same(record['omitted_configuration'],CONFIGS[i],r['fit_id']+':omitted label')
            if not record['success']:continue
            mean=(71*w.mean-w.samples[i])/70;p=np.array(record['p']);residual=w.residual(p,r['states'],mean)
            same(float(residual@residual),record['q2'],r['fit_id']+':jackknife objective',rtol=1e-7)
            d2=fm.diagnostics(w,r['states'],p,record['q2'],mean)
            same(d2['identified'],record['identified'],r['fit_id']+':jackknife diagnostics')
            if d2['flags']!=record['flags']:raise ValueError('Jackknife flags '+r['fit_id'])
            converged+=1;resolved+=int(record['identified']);energies.append(p[1 if r['states']==1 else 2]);replayed_jackknife+=1
        same(r['jackknife_converged'],converged,r['fit_id']+':jackknife count')
        same(r['jackknife_resolved'],resolved,r['fit_id']+':resolved count')
        same(r['jackknife_reliable'],converged==71 and resolved==71,r['fit_id']+':reliable audit')
        se=float(np.sqrt(70/71*np.sum((np.array(energies)-np.mean(energies))**2))) if converged==71 else None
        same(r['E0_jackknife_se_formal'],se,r['fit_id']+':jackknife error')
        if r['states']==2:
            one=lookup[r['Nv'],r['R'],r['tau_min'],r['tau_max'],1]
            if r['q2']>one['q2']+1e-6:raise ValueError('Nested one-state control')
            model_comparison.append(dict(Nv=r['Nv'],R=r['R'],tau_min=r['tau_min'],tau_max=r['tau_max'],
                single_E0=one['E0'],two_E0=r['E0'],q2_single=one['q2'],q2_two=r['q2'],delta_q2=one['q2']-r['q2'],
                two_state_identified=r['identified'],two_state_flags=r['flags']))
        boot=solution.get('bootstrap',[])
        wanted=r['states']==1 and r['Nv'] in (40,100) and r['R'] in (4,12,24,32) and r['tau_min']==8 and r['tau_max']==32
        if len(boot)!=(399 if wanted else 0):raise ValueError('Bootstrap inventory '+r['fit_id'])
        if not wanted:continue
        null=w.samples-w.mean+fm.model(parameters(r),w.t,1);valid_q=[]
        np.testing.assert_allclose(null.mean(0),fm.model(parameters(r),w.t,1),rtol=1e-12,atol=1e-14)
        for draw,record in enumerate(boot):
            same(record['draw'],draw,r['fit_id']+':bootstrap draw')
            if not record['success']:continue
            sample=fm.Window(null[chosen[draw]],w.times,scale=1.)
            residual=sample.residual(np.array(record['p']),1)
            same(float(residual@residual),record['q2'],r['fit_id']+':bootstrap objective',rtol=1e-7)
            valid_q.append(record['q2']);replayed_bootstrap+=1
        value=(sum(q>=r['q2'] for q in valid_q)+1)/400 if len(valid_q)==399 else None
        bootstrap.append(dict(fit_id=r['fit_id'],Nv=r['Nv'],R=r['R'],tau_min=8,tau_max=32,draws=399,successful_draws=len(valid_q),
            p_recentered_bootstrap=value,monte_carlo_se=math.sqrt(value*(1-value)/400) if value is not None else None,
            p_gaussian_reference=r['p_gaussian_reference'],p_hotelling_linearized_reference=r['p_hotelling_linearized_reference']))
    check_csv(RESULTS/'bootstrap_goodness_71configs.csv',bootstrap)
    check_csv(RESULTS/'model_comparison_71configs.csv',model_comparison)
    exports=derive(a,rows,solutions)
    for name,records in exports.items():check_csv(RESULTS/name,records)
    expected_counts=dict(configuration_count=71,one_state_fits=384,two_state_fits=160,
        identified_two_state_fits=sum(r['states']==2 and r['identified'] for r in rows),
        nonlinear_jackknife_requested=sum(r['jackknife_requested'] for r in rows),
        nonlinear_jackknife_performed=sum(r['jackknife_performed'] for r in rows),
        nonlinear_jackknife_all_replicates_resolved=sum(r['jackknife_reliable'] for r in rows))
    for k,v in expected_counts.items():same(summary[k],v,'summary '+k)
    if summary['configuration_labels']!=CONFIGS or summary['bootstrap_null_tests']!=bootstrap:raise ValueError('Summary data')
    ratios=[r['jackknife_over_naive'] for r in exports['parameter_error_audit_71configs.csv'] if r['reliable']]
    np.testing.assert_allclose(summary['reliable_audit_jackknife_over_naive_range'],[min(ratios),max(ratios)])
    if summary['derived_exports']!={k:len(v) for k,v in exports.items()}:raise ValueError('Derived export counts')
    print('PASS full71 exponential replay: 544 fits; %d archived starts; %d jackknife fits; %d covariance-resampled null fits; paired contrasts reproduced.'%(replayed_starts,replayed_jackknife,replayed_bootstrap))
    print('Archived solver solutions replayed; optimizer search is not repeated by --verify.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',required=True)
    parser.parse_args();verify()
