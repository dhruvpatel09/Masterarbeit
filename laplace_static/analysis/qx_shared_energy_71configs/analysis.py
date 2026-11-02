#!/usr/bin/env python3
"""Common-energy diagnostics on existing qX configuration measurements."""
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import beta,chi2
import joint_models as jm

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RESULTS=HERE/'results'
CONFIGS=list(range(100,801,10));CUTOFFS=[40,60,80,100]
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def write_json(path,value):
    path=Path(path);tmp=Path(str(path)+'.pending')
    with tmp.open('w') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def write_csv(path,rows):
    path=Path(path);tmp=Path(str(path)+'.pending')
    with tmp.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows);f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def load_data():
    source=json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name,expected in source['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Changed published dependency: '+name)
    raw=next(k for k in source['files'] if k.endswith('/correlators.csv'))
    a=np.full((71,4,33,12),np.nan);seen=set()
    with (ROOT/raw).open(newline='') as f:
        for r in csv.DictReader(f):
            if r['scheme']!='fixed' or r['axis']!='-1' or r['piter']!='10':continue
            key=tuple(int(r[k]) for k in ('cfg','Nv','tau','R'))
            if key in seen:raise ValueError('Duplicate combined primary')
            seen.add(key);c,n,t,R=key;a[CONFIGS.index(c),CUTOFFS.index(n),t,RADII.index(R)]=float(r['L_RE'])
    if len(seen)!=112464 or not np.all(np.isfinite(a)):raise ValueError('Incomplete 71-configuration source table')
    old=next(k for k in source['files'] if k.endswith('/fits_71configs.csv'))
    with (ROOT/old).open(newline='') as f:
        seeds={(int(r['Nv']),int(r['R']),int(r['tau_min'])):float(r['E0']) for r in csv.DictReader(f)
               if r['states']=='1' and r['tau_max']=='32'}
    return a,seeds


def cases():
    plan=json.loads((HERE/'ANALYSIS_PLAN.json').read_text());result=[]
    for design in plan['designs']:
        for R in plan['distances']:
            for lo in plan['starts']:
                result.append(dict(case_id='%s_R%d_t%d_32'%(design['name'],R,lo),design=design['name'],
                    cutoffs=design['cutoffs'],stride=design['time_stride'],R=R,tau_min=lo,tau_max=32,
                    times=list(range(lo,33,design['time_stride']))))
    return result


def selected_samples(case,a):
    # Separate advanced indexing operations preserve the configuration axis.
    return a[:,[CUTOFFS.index(n) for n in case['cutoffs']],:,:][:,:,case['times'],:][:,:,:,RADII.index(case['R'])]


def window_for(case,a):
    return jm.JointWindow(selected_samples(case,a),np.array(case['times']),case['cutoffs'])


def fit_row(case,w,kind,best,candidates):
    p=np.array(best['p']);d=jm.diagnostics(w,kind,p,best['q2'],candidates=candidates);local=d['local']
    ei=2*w.m if kind=='shared_two' else w.m;shared=kind!='independent_one'
    row=dict(fit_id=case['case_id']+'_'+kind,case_id=case['case_id'],design=case['design'],R=case['R'],
        tau_min=case['tau_min'],tau_max=32,time_stride=case['stride'],model=kind,configurations=71,
        sources=w.m,time_points_per_source=w.nt,covariance_dimension=len(w.mean),correlation_condition=w.condition,
        E0=float(p[ei]) if shared else None,E0_naive_se=float(local['se'][ei]) if local and shared else None,
        gap=float(p[-1]) if kind=='shared_two' else None,gap_naive_se=float(local['se'][-1]) if local and kind=='shared_two' else None,
        q2=best['q2'],dof=d['dof'],p_gaussian_reference=d['p_gaussian_reference'],local_checks_pass=d['local_checks_pass'],
        flags='|'.join(d['flags']) or 'none',near_optimal_E_span=d['near_E_span'],
        jacobian_condition=local['jacobian_condition'] if local else None,curvature_condition=local['curvature_condition'] if local else None,
        minimum_mean_over_se=float(np.min(w.mean/w.se)),nonpositive_means=int(np.sum(w.mean<=0)),
        successful_starts=sum(c['success'] for c in candidates),total_starts=len(candidates),best_start=best['start'],
        jackknife_requested=False,jackknife_performed=False,jackknife_converged=0,jackknife_resolved=0,
        E0_jackknife_se_formal=None,jackknife_all_local_checks=False)
    for n in CUTOFFS:
        if n in w.cutoffs:
            i=w.cutoffs.index(n);j=ei if shared else w.m+i
            row['energy_Nv%d'%n]=float(p[j]);row['energy_se_Nv%d'%n]=float(local['se'][j]) if local else None
        else:row['energy_Nv%d'%n]=row['energy_se_Nv%d'%n]=None
    return row,d


def tail_summary(values,observed):
    n=len(values);k=sum(v>=observed for v in values);value=(k+1)/(n+1)
    low=0. if k==0 else float(beta.ppf(.025,k,n-k+1))
    high=1. if k==n else float(beta.ppf(.975,k+1,n-k))
    return dict(exceedances=k,p_add_one=value,conditional_tail_low95=low,conditional_tail_high95=high,resolution=1/(n+1))


def ratio_estimate(numerator,denominator):
    numerator=np.asarray(numerator,float);denominator=np.asarray(denominator,float)
    mean=float(denominator.mean())
    if mean==0 or not np.isfinite(mean):raise ValueError('Undefined ratio reference')
    value=float(numerator.mean()/mean);fluctuation=(numerator-value*denominator)/mean
    return value,float(np.std(fluctuation,ddof=1)/np.sqrt(len(numerator)))


def bootstrap_row(case,full,records):
    row=dict(case_id=case['case_id'],design=case['design'],R=case['R'],tau_min=case['tau_min'],tau_max=32,
        time_stride=case['stride'],draws=399,successful_draws=sum(r['success'] for r in records))
    one=full['shared_one'];ind=full['independent_one'];delta=one['q2']-ind['q2']
    row.update(q2_shared=one['q2'],delta_q2=delta,shared_fit_local_checks=one['local_checks_pass'],independent_fit_local_checks=ind['local_checks_pass'],
        p_gaussian_shared=one['p_gaussian_reference'],p_gaussian_delta_reference=float(chi2.sf(max(delta,0.),len(case['cutoffs'])-1)) if one['local_checks_pass'] and ind['local_checks_pass'] else None)
    for label,values,observed in [('shared',[r['shared']['q2'] for r in records if r['success']],one['q2']),
                                  ('delta',[r['shared']['q2']-r['independent']['q2'] for r in records if r['success']],delta)]:
        summary=tail_summary(values,observed) if len(values)==399 else {k:None for k in ('exceedances','p_add_one','conditional_tail_low95','conditional_tail_high95','resolution')}
        for k,v in summary.items():row[label+'_'+k]=v
    return row


def derived(a,rows,solutions):
    bycase={};lookup={r['fit_id']:r for r in rows}
    for r in rows:bycase.setdefault(r['case_id'],{})[r['model']]=r
    contrasts=[];comparison=[];source_contrasts=[]
    for case in cases():
        full=bycase[case['case_id']];one,ind,two=(full[k] for k in jm.KINDS)
        comparison.append(dict(case_id=case['case_id'],design=case['design'],R=case['R'],tau_min=case['tau_min'],tau_max=32,
            time_stride=case['stride'],shared_one_E0=one['E0'],shared_one_se=one['E0_naive_se'],shared_two_E0=two['E0'],shared_two_se=two['E0_naive_se'],
            shared_two_local_checks=two['local_checks_pass'],shared_two_flags=two['flags'],q2_shared_one=one['q2'],q2_independent_one=ind['q2'],
            delta_common_constraint=one['q2']-ind['q2'],q2_shared_two=two['q2'],
            delta_extra_state=one['q2']-two['q2'],extra_state_p_value=None))
        w=window_for(case,a);fluct={}
        ip=np.array(solutions[ind['fit_id']]['p'])
        try:ilocal=w.local_errors(ip,'independent_one')
        except (ValueError,np.linalg.LinAlgError):ilocal=None
        for i,n in enumerate(case['cutoffs']):
            for j,n2 in enumerate(case['cutoffs']):
                if j<=i:continue
                difference=float(ip[w.m+i]-ip[w.m+j]);source_se=None
                if ilocal is not None:
                    f=ilocal['fluctuations'][:,w.m+i]-ilocal['fluctuations'][:,w.m+j]
                    source_se=float(np.std(f,ddof=1)/np.sqrt(71))
                source_contrasts.append(dict(case_id=case['case_id'],design=case['design'],R=case['R'],tau_min=case['tau_min'],tau_max=32,
                    Nv_first=n,Nv_second=n2,delta_energy=difference,paired_naive_se=source_se,
                    paired_z_diagnostic=difference/source_se if source_se else None,independent_model_local_checks=ind['local_checks_pass']))
        for kind in ('shared_one','shared_two'):
            p=np.array(solutions[full[kind]['fit_id']]['p']);ei=2*w.m if kind=='shared_two' else w.m
            try:fluct[kind]=w.local_errors(p,kind)['fluctuations'][:,ei]
            except (ValueError,np.linalg.LinAlgError):fluct[kind]=None
        se=float(np.std(fluct['shared_two']-fluct['shared_one'],ddof=1)/np.sqrt(71)) if all(v is not None for v in fluct.values()) else None
        contrasts.append(dict(case_id=case['case_id'],design=case['design'],R=case['R'],tau_min=case['tau_min'],tau_max=32,
            contrast='shared_two_minus_shared_one',delta_E0=two['E0']-one['E0'],paired_naive_se=se,
            both_local_checks_pass=two['local_checks_pass'] and one['local_checks_pass']))
    shapes=[]
    for i,n in enumerate(CUTOFFS):
        for j,R in enumerate(RADII):
            for tau,ref_R,ref_data in ((0,1,a[:,i,0,0]),(24,R,a[:,i,0,j])):
                ratio,se=ratio_estimate(a[:,i,tau,j],ref_data)
                shapes.append(dict(Nv=n,R=R,numerator_tau=tau,reference_R=ref_R,reference_tau=0,
                    configurations=71,ratio_of_means=ratio,paired_naive_se=se))
    return {'source_compatibility_71configs.csv':comparison,'paired_model_contrasts_71configs.csv':contrasts,
            'paired_source_energy_contrasts_71configs.csv':source_contrasts,'source_shape_diagnostics_71configs.csv':shapes}


def audit_requested(case,row):
    return ((case['design'] in ('pair_dense','all_four') and case['tau_min']==12 and row['model'] in ('shared_one','shared_two'))
            or (row['model']=='shared_two' and row['local_checks_pass']))


def audit_fit(case,w,row,solutions):
    kind=row['model'];row['jackknife_performed']=True;records=[];energy=[];resolved=0;ei=2*w.m if kind=='shared_two' else w.m
    for i,cfg in enumerate(CONFIGS):
        mean=(71*w.mean-w.samples[i])/70;best,_=jm.solve(w,kind,starts=[solutions[row['fit_id']]['p']],mean=mean)
        record=dict(omitted_configuration=cfg,success=best is not None,p=best['p'] if best else None,q2=best['q2'] if best else None)
        if best:
            d=jm.diagnostics(w,kind,best['p'],best['q2'],mean);record.update(local_checks_pass=d['local_checks_pass'],flags=d['flags'])
            energy.append(best['p'][ei]);resolved+=int(d['local_checks_pass'])
        records.append(record)
    solutions[row['fit_id']]['jackknife']=records;row['jackknife_converged']=len(energy);row['jackknife_resolved']=resolved
    if len(energy)==71:row['E0_jackknife_se_formal']=float(np.sqrt(70/71*np.sum((np.array(energy)-np.mean(energy))**2)))
    row['jackknife_all_local_checks']=len(energy)==71 and resolved==71
    print('JACKKNIFE',row['fit_id'],len(energy),'converged,',resolved,'local checks pass',flush=True)


def main():
    start=time.monotonic();a,seeds=load_data();rows=[];solutions={};bycase={};window_cache={}
    for case in cases():
        w=window_for(case,a);window_cache[case['case_id']]=w;full={};best_models={}
        for kind in jm.KINDS:
            best,candidates=jm.solve(w,kind,shared=best_models.get('shared_one'),independent=best_models.get('independent_one'),
                seed_energies=[seeds[n,case['R'],case['tau_min']] for n in case['cutoffs']])
            if best is None:raise ValueError('No converged start '+case['case_id']+' '+kind)
            if kind!='shared_one' and best['q2']>best_models['shared_one']['q2']+1e-6:raise ValueError('Nested shared-one control failed')
            row,d=fit_row(case,w,kind,best,candidates);rows.append(row);full[kind]=row;best_models[kind]=best
            solutions[row['fit_id']]=dict(p=best['p'],candidates=candidates,local_covariance=d['local']['covariance'].tolist() if d['local'] else None,jackknife=[])
        bycase[case['case_id']]=full
        print('FITS',case['case_id'],'common E',round(full['shared_one']['E0'],6),'delta q2',round(full['shared_one']['q2']-full['independent_one']['q2'],2),flush=True)
    # The covariance and per-source normalization remain fixed in this local-branch audit.
    for case in cases():
        w=window_cache[case['case_id']]
        for kind in ('shared_one','shared_two'):
            row=bycase[case['case_id']][kind];row['jackknife_requested']=audit_requested(case,row)
            if row['jackknife_requested'] and row['local_checks_pass']:audit_fit(case,w,row,solutions)
    indices=np.random.RandomState(20261103).randint(0,71,size=(399,71));write_json(RESULTS/'bootstrap_indices.json',dict(seed=20261103,configuration_labels=CONFIGS,indices=indices.tolist()))
    bootstrap=[]
    for case in cases():
        if case['tau_min']!=12:continue
        w=window_cache[case['case_id']];full=bycase[case['case_id']];one=solutions[full['shared_one']['fit_id']]['p'];ind=solutions[full['independent_one']['fit_id']]['p']
        null=w.samples-w.mean+jm.model(one,w.dt,w.m,'shared_one');records=[]
        for draw,chosen in enumerate(indices):
            try:
                sample=jm.JointWindow(null[chosen].reshape(71,w.m,w.nt),w.source_times,w.cutoffs,scales=np.ones(w.m))
                b1,_=jm.solve(sample,'shared_one',starts=[one,np.r_[np.array(one[:w.m])*.9,one[-1]*1.1]])
                if b1 is None:raise ValueError('Shared null refit failed')
                nested=np.r_[b1['p'][:w.m],np.repeat(b1['p'][w.m],w.m)]
                bi,_=jm.solve(sample,'independent_one',starts=[nested,ind])
                if bi is None or bi['q2']>b1['q2']+1e-6:raise ValueError('Independent null refit/nested control failed')
                records.append(dict(draw=draw,success=True,shared={'p':b1['p'],'q2':b1['q2']},independent={'p':bi['p'],'q2':bi['q2']}))
            except (ValueError,np.linalg.LinAlgError) as exc:records.append(dict(draw=draw,success=False,error=str(exc)))
        solutions[full['shared_one']['fit_id']]['bootstrap']=records
        record=bootstrap_row(case,full,records);bootstrap.append(record)
        print('BOOTSTRAP',case['case_id'],record['successful_draws'],'valid; shared p',record['shared_p_add_one'],'constraint p',record['delta_p_add_one'],flush=True)
    write_json(RESULTS/'solutions.json',solutions);write_csv(RESULTS/'joint_fits_71configs.csv',rows);write_csv(RESULTS/'bootstrap_compatibility_71configs.csv',bootstrap)
    exports=derived(a,rows,solutions)
    for name,records in exports.items():write_csv(RESULTS/name,records)
    audits=[r for r in rows if r['jackknife_requested']];ratios=[r['E0_jackknife_se_formal']/r['E0_naive_se'] for r in audits if r['jackknife_all_local_checks']]
    summary=dict(success=True,configuration_count=71,configuration_labels=CONFIGS,case_count=len(cases()),fit_count=len(rows),
        source_time_covariance_max_dimension=max(r['covariance_dimension'] for r in rows),
        maximum_correlation_condition=max(r['correlation_condition'] for r in rows),
        locally_resolved_shared_two=sum(r['model']=='shared_two' and r['local_checks_pass'] for r in rows),
        shared_two_fit_count=sum(r['model']=='shared_two' for r in rows),
        jackknife_requested=len(audits),jackknife_performed=sum(r['jackknife_performed'] for r in audits),jackknife_all_local_checks=sum(r['jackknife_all_local_checks'] for r in audits),
        jackknife_over_naive_range=[min(ratios),max(ratios)] if ratios else None,
        bootstrap_tests=len(bootstrap),bootstrap_results=bootstrap,derived_exports={k:len(v) for k,v in exports.items()},
        interpretation='Common-energy/source-compatibility diagnostics, not selected ground-state energies or a final potential. All signed configuration measurements retained; local parameter checks differ from fit acceptance.',
        seconds=time.monotonic()-start,versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    write_json(RESULTS/'SUMMARY.json',summary)
    print('PASS joint source analysis:',len(rows),'fits;',summary['locally_resolved_shared_two'],'locally resolved two-state fits',flush=True)


if __name__=='__main__':main()
