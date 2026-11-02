#!/usr/bin/env python3
"""Replay archived joint-fit objectives, covariance propagation and resamples."""
import argparse
import csv
import json
import math
import numpy as np
import joint_models as jm
from analysis import RESULTS,CONFIGS,cases,load_data,window_for,selected_samples,fit_row,bootstrap_row,derived,audit_requested


def load_rows():
    integer={'R','tau_min','tau_max','time_stride','configurations','sources','time_points_per_source','covariance_dimension','dof',
        'nonpositive_means','successful_starts','total_starts','best_start','jackknife_converged','jackknife_resolved'}
    boolean={'local_checks_pass','jackknife_requested','jackknife_performed','jackknife_all_local_checks'}
    strings={'fit_id','case_id','design','model','flags'}
    with (RESULTS/'joint_fits_71configs.csv').open(newline='') as f:
        return [{k:(int(v) if k in integer else v=='True' if k in boolean else v if k in strings else float(v) if v else None)
            for k,v in r.items()} for r in csv.DictReader(f)]


def same(a,b,label,rtol=3e-5,atol=1e-8):
    if b is None or isinstance(b,(str,bool)):
        if a!=b:raise ValueError('Mismatch '+label+': '+str((a,b)))
    elif a is None or not math.isfinite(float(a)) or not np.isclose(float(a),float(b),rtol=rtol,atol=atol):
        raise ValueError('Numerical mismatch '+label+': '+str((a,b)))


def check_json(actual,expected,label='summary',rtol=1e-10,atol=1e-13):
    """Keep JSON structure/counts exact, allowing roundoff in computed floats."""
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or set(actual)!=set(expected):raise ValueError('Summary fields differ: '+label)
        for field,value in expected.items():check_json(actual[field],value,label+'.'+field,rtol,atol)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):raise ValueError('Summary list length differs: '+label)
        for i,(a,b) in enumerate(zip(actual,expected)):check_json(a,b,label+'['+str(i)+']',rtol,atol)
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(int,float)) or not math.isfinite(actual) or not math.isfinite(expected):
            raise ValueError('Invalid summary number: '+label)
        if not math.isclose(actual,expected,rel_tol=rtol,abs_tol=atol):
            raise ValueError('Summary numerical mismatch '+label+': '+str((actual,expected)))
    elif type(actual) is not type(expected) or actual!=expected:
        raise ValueError('Summary value differs '+label+': '+str((actual,expected)))


def check_csv(path,expected):
    with path.open(newline='') as f:actual=list(csv.DictReader(f))
    if len(actual)!=len(expected):raise ValueError('Row inventory '+path.name)
    for i,(a,b) in enumerate(zip(actual,expected)):
        if set(a)!=set(b):raise ValueError('Columns '+path.name)
        for k,v in b.items():
            x=None if v is None and a[k]=='' else a[k]=='True' if isinstance(v,bool) else a[k] if isinstance(v,str) else float(a[k])
            same(x,v,path.name+':'+str(i)+':'+k)


def verify():
    a,_=load_data();rows=load_rows();solutions=json.loads((RESULTS/'solutions.json').read_text())
    expected={case['case_id']+'_'+kind for case in cases() for kind in jm.KINDS}
    if len(rows)!=162 or {r['fit_id'] for r in rows}!=expected or set(solutions)!=expected:raise ValueError('Joint fit inventory')
    lookup={r['fit_id']:r for r in rows};case_models={};draw_count=0;refits=0;jk_count=0;starts=0;bootstrap=[]
    saved_indices=json.loads((RESULTS/'bootstrap_indices.json').read_text());indices=np.array(saved_indices['indices'])
    if saved_indices['seed']!=20261103 or saved_indices['configuration_labels']!=CONFIGS:raise ValueError('Resampling cohort')
    np.testing.assert_array_equal(indices,np.random.RandomState(20261103).randint(0,71,size=(399,71)))
    for case in cases():
        w=window_for(case,a);full={}
        # Explicit loops independently check the configuration/source/time indexing.
        raw=selected_samples(case,a)
        for i,n in enumerate(case['cutoffs']):
            from analysis import CUTOFFS,RADII
            np.testing.assert_array_equal(raw[:,i],np.array([[a[c,CUTOFFS.index(n),t,RADII.index(case['R'])] for t in case['times']] for c in range(71)]))
        if len(w.mean)>28:raise ValueError('Covariance dimension exceeded analysis plan')
        for kind in jm.KINDS:
            fit_id=case['case_id']+'_'+kind;r=lookup[fit_id];s=solutions[fit_id];candidates=s['candidates']
            expected_starts=5 if kind=='shared_one' else 6 if kind=='independent_one' else 14
            if len(candidates)!=expected_starts:raise ValueError('Candidate inventory '+fit_id)
            for c in candidates:
                if c['p'] is None:continue
                residual=w.residual(c['p'],kind);same(float(residual@residual),c['q2'],fit_id+':candidate objective',rtol=1e-7)
                starts+=1
            valid=[c for c in candidates if c['success'] and c['p'] is not None and np.isfinite(c['q2'])]
            best=min(valid,key=lambda c:c['q2']);np.testing.assert_array_equal(best['p'],s['p'])
            expected_row,d=fit_row(case,w,kind,best,candidates)
            for field,value in expected_row.items():
                if field.startswith('jackknife_') or field=='E0_jackknife_se_formal':continue
                same(r[field],value,fit_id+':'+field,rtol=3e-5 if r['local_checks_pass'] else 5e-3)
            if d['local'] is None:same(s['local_covariance'],None,fit_id+':curvature')
            else:np.testing.assert_allclose(s['local_covariance'],d['local']['covariance'],rtol=5e-3,atol=1e-8)
            requested=audit_requested(case,r)
            performed=requested and r['local_checks_pass']
            same(r['jackknife_requested'],requested,fit_id+':audit request');same(r['jackknife_performed'],performed,fit_id+':audit performed')
            records=s['jackknife'];energy=[];resolved=0;ei=2*w.m if kind=='shared_two' else w.m
            if len(records)!=(71 if performed else 0):raise ValueError('Audit inventory '+fit_id)
            for i,record in enumerate(records):
                same(record['omitted_configuration'],CONFIGS[i],fit_id+':omitted label')
                if not record['success']:continue
                mean=(71*w.mean-w.samples[i])/70;residual=w.residual(record['p'],kind,mean)
                same(float(residual@residual),record['q2'],fit_id+':jackknife objective',rtol=1e-7)
                d2=jm.diagnostics(w,kind,record['p'],record['q2'],mean)
                same(d2['local_checks_pass'],record['local_checks_pass'],fit_id+':jackknife check')
                if d2['flags']!=record['flags']:raise ValueError('Jackknife flags '+fit_id)
                energy.append(record['p'][ei]);resolved+=int(record['local_checks_pass']);jk_count+=1
            same(r['jackknife_converged'],len(energy),fit_id+':audit count');same(r['jackknife_resolved'],resolved,fit_id+':resolved count')
            same(r['jackknife_all_local_checks'],len(energy)==resolved==71,fit_id+':all audit checks')
            se=float(np.sqrt(70/71*np.sum((np.array(energy)-np.mean(energy))**2))) if len(energy)==71 else None
            same(r['E0_jackknife_se_formal'],se,fit_id+':audit error')
            full[kind]=r
        case_models[case['case_id']]=full
        for kind in ('independent_one','shared_two'):
            if full[kind]['q2']>full['shared_one']['q2']+1e-6:raise ValueError('Nested model objective')
        records=solutions[full['shared_one']['fit_id']].get('bootstrap',[])
        selected=case['tau_min']==12
        if len(records)!=(399 if selected else 0):raise ValueError('Null sample inventory '+case['case_id'])
        if not selected:continue
        p=solutions[full['shared_one']['fit_id']]['p'];prediction=jm.model(p,w.dt,w.m,'shared_one');null=w.samples-w.mean+prediction
        np.testing.assert_allclose(null.mean(0),prediction,rtol=1e-12,atol=1e-14)
        for draw,record in enumerate(records):
            same(record['draw'],draw,case['case_id']+':draw')
            if not record['success']:continue
            sample=jm.JointWindow(null[indices[draw]].reshape(71,w.m,w.nt),w.source_times,w.cutoffs,scales=np.ones(w.m))
            for label,kind in (('shared','shared_one'),('independent','independent_one')):
                residual=sample.residual(record[label]['p'],kind)
                same(float(residual@residual),record[label]['q2'],case['case_id']+':null objective',rtol=1e-7);refits+=1
            if record['independent']['q2']>record['shared']['q2']+1e-6:raise ValueError('Null nested control')
            draw_count+=1
        bootstrap.append(bootstrap_row(case,full,records))
    check_csv(RESULTS/'bootstrap_compatibility_71configs.csv',bootstrap)
    exports=derived(a,rows,solutions)
    for name,records in exports.items():check_csv(RESULTS/name,records)
    summary=json.loads((RESULTS/'SUMMARY.json').read_text());audits=[r for r in rows if r['jackknife_requested']]
    counts=dict(configuration_count=71,case_count=54,fit_count=162,
        source_time_covariance_max_dimension=max(r['covariance_dimension'] for r in rows),maximum_correlation_condition=max(r['correlation_condition'] for r in rows),
        locally_resolved_shared_two=sum(r['model']=='shared_two' and r['local_checks_pass'] for r in rows),shared_two_fit_count=54,
        jackknife_requested=len(audits),jackknife_performed=sum(r['jackknife_performed'] for r in audits),
        jackknife_all_local_checks=sum(r['jackknife_all_local_checks'] for r in audits),bootstrap_tests=len(bootstrap))
    for k,v in counts.items():
        if isinstance(v,int):check_json(summary[k],v,'summary.'+k)
        else:same(summary[k],v,'summary '+k)
    check_json(summary['configuration_labels'],CONFIGS,'summary.configuration_labels')
    check_json(summary['bootstrap_results'],bootstrap,'summary.bootstrap_results')
    check_json(summary['derived_exports'],{k:len(v) for k,v in exports.items()},'summary.derived_exports')
    ratios=[r['E0_jackknife_se_formal']/r['E0_naive_se'] for r in audits if r['jackknife_all_local_checks']]
    np.testing.assert_allclose(summary['jackknife_over_naive_range'],[min(ratios),max(ratios)])
    print('PASS shared-energy replay: 71 configurations; 54 joint cases; 162 fits; %d candidate solutions; %d jackknife solutions; %d null draws / %d refits.'%(starts,jk_count,draw_count,refits))
    print('All cross-source/time covariances and paired contrasts replayed; optimizer searches not repeated by --verify.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',required=True,action='store_true');parser.parse_args();verify()
