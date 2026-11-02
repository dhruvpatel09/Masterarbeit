#!/usr/bin/env python3
"""Fit the existing full71 correlators; no gauge-field measurements required."""
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path
import numpy as np
import scipy
import fit_models as fm

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RESULTS=HERE/'results'
CONFIGS=list(range(100,801,10))
CUTOFFS=[40,60,80,100]
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]
TWO_R=[4,12,16,24,32]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def write_json(path,value):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.pending')
    with tmp.open('w') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def write_csv(path,rows):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.pending')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows);f.flush();os.fsync(f.fileno())
    tmp.replace(path)


def load_data():
    receipt=json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name,expected in receipt['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Published dependency changed: '+name)
    name=next(k for k in receipt['files'] if k.endswith('/correlators.csv'))
    a=np.full((71,4,33,12),np.nan);seen=set()
    with (ROOT/name).open(newline='') as f:
        for r in csv.DictReader(f):
            if r['scheme']=='fixed' and r['axis']=='-1' and r['piter']=='10':
                key=tuple(int(r[k]) for k in ('cfg','Nv','tau','R'))
                if key in seen:raise ValueError('Duplicate primary');
                seen.add(key);c,n,t,R=key
                a[CONFIGS.index(c),CUTOFFS.index(n),t,RADII.index(R)]=float(r['L_RE'])
    if len(seen)!=112464 or not np.all(np.isfinite(a)):raise ValueError('Unmatched full71 primaries')
    return a


def key(n,R,lo,hi,states):return 'Nv%d_R%d_t%d_%d_s%d'%(n,R,lo,hi,states)


def serialize_candidate(c):
    return dict(c,p=c['p'].tolist() if c['p'] is not None else None)


def parameters(row):
    return np.array([row['A0_scaled'],row['E0']]) if row['states']==1 else np.array([
        row['A0_scaled'],row['A1_scaled'],row['E0'],row['gap']])


def make_row(n,R,lo,hi,states,window,best,candidates):
    p=best['p'];d=fm.diagnostics(window,states,p,best['q2'],candidates=candidates)
    local=d['local'];ei=1 if states==1 else 2
    row=dict(fit_id=key(n,R,lo,hi,states),Nv=n,R=R,tau_min=lo,tau_max=hi,states=states,
        configurations=71,scale_L=window.scale,A0_scaled=float(p[0]),A1_scaled=float(p[1]) if states==2 else 0.,
        E0=float(p[ei]),gap=float(p[3]) if states==2 else None,
        E1=float(p[2]+p[3]) if states==2 else None,
        E0_naive_se=float(local['se'][ei]) if local else None,
        gap_naive_se=float(local['se'][3]) if local and states==2 else None,
        ground_fraction_at_tmin=d['ground_fraction'],q2=best['q2'],dof=d['dof'],
        p_gaussian_reference=d['p_gaussian_reference'],p_hotelling_linearized_reference=d['p_hotelling_linearized_reference'],
        identified=d['identified'],flags='|'.join(d['flags']) or 'none',
        scaled_jacobian_condition=local['scaled_jacobian_condition'] if local else None,
        scaled_hessian_condition=local['scaled_hessian_condition'] if local else None,
        near_optimal_start_E0_span=d['near_optimal_E0_span'],correlation_condition=window.condition,
        nonpositive_means_in_window=int(np.sum(window.mean<=0)),minimum_mean_over_se=float(np.min(window.mean/window.se)),
        successful_starts=sum(c['success'] for c in candidates),total_starts=len(candidates),
        best_start=best['start'],jackknife_requested=False,jackknife_performed=False,
        jackknife_converged=0,jackknife_resolved=0,E0_jackknife_se_formal=None,jackknife_reliable=False)
    return row,d


def main():
    start=time.monotonic();a=load_data();rows=[];solutions={};internal={};audit=[]
    jobs=[(n,R,lo,hi) for n in CUTOFFS for R in RADII for lo in (4,8,12,16) for hi in (24,32)]
    for number,(n,R,lo,hi) in enumerate(jobs):
        samples=a[:,CUTOFFS.index(n),lo:hi+1,RADII.index(R)]
        window=fm.Window(samples,np.arange(lo,hi+1))
        single,candidates=fm.solve(window,1)
        if single is None:raise ValueError('No converged single-exponential start')
        for states in (1,2) if R in TWO_R else (1,):
            best,starts=(single,candidates) if states==1 else fm.solve(window,2,single=single['p'])
            if best is None:raise ValueError('No converged two-exponential start')
            if states==2 and best['q2']>single['q2']+1e-6:raise ValueError('Two-state search failed nested one-state control')
            row,d=make_row(n,R,lo,hi,states,window,best,starts)
            rows.append(row);internal[row['fit_id']]=(window,d)
            solutions[row['fit_id']]=dict(candidates=[serialize_candidate(c) for c in starts],
                best_start=best['start'],jackknife=[],local_covariance=d['local']['covariance'].tolist() if d['local'] else None)
        if (number+1)%24==0:print('FITS %d/%d one-state windows; two-state controls included'%(number+1,len(jobs)),flush=True)
    write_json(RESULTS/'solutions.json',solutions);write_csv(RESULTS/'fits_71configs.csv',rows)
    # Audit the nonlinear error propagation without changing the 71-sample unit.
    for row in rows:
        requested=row['Nv'] in (40,100) and row['R'] in (4,12,24,32) and row['tau_min'] in (8,12) and row['tau_max']==32
        row['jackknife_requested']=requested
        if not requested or not row['identified']:continue
        row['jackknife_performed']=True;window,_=internal[row['fit_id']];p=parameters(row);estimates=[];resolved=0
        for i,cfg in enumerate(CONFIGS):
            mean=(71*window.mean-window.samples[i])/70
            best,candidates=fm.solve(window,row['states'],starts=[p],mean=mean)
            success=best is not None
            record=dict(omitted_configuration=cfg,success=success,p=best['p'].tolist() if success else None,q2=best['q2'] if success else None)
            if success:
                d=fm.diagnostics(window,row['states'],best['p'],best['q2'],mean)
                record['identified']=d['identified'];record['flags']=d['flags'];resolved+=int(d['identified'])
                estimates.append(float(best['p'][1 if row['states']==1 else 2]))
            solutions[row['fit_id']]['jackknife'].append(record)
        row['jackknife_converged']=len(estimates);row['jackknife_resolved']=resolved
        if len(estimates)==71:
            x=np.array(estimates);row['E0_jackknife_se_formal']=float(np.sqrt(70/71*np.sum((x-x.mean())**2)))
        row['jackknife_reliable']=len(estimates)==71 and resolved==71
        print('JACKKNIFE',row['fit_id'],len(estimates),'converged;',resolved,'resolved',flush=True)
    # Recenter the empirical configuration distribution under the fitted null.
    rng=np.random.RandomState(20261102);indices=rng.randint(0,71,size=(399,71));bootstrap=[]
    write_json(RESULTS/'bootstrap_indices.json',dict(configuration_labels=CONFIGS,indices=indices.tolist(),seed=20261102))
    for row in rows:
        if row['states']!=1 or row['Nv'] not in (40,100) or row['R'] not in (4,12,24,32) or row['tau_min']!=8 or row['tau_max']!=32:continue
        window,_=internal[row['fit_id']];p=parameters(row);prediction=fm.model(p,window.t,1)
        centred=window.samples-window.mean+prediction;records=[]
        for draw,chosen in enumerate(indices):
            try:
                sample=fm.Window(centred[chosen],window.times,scale=1.)
                best,_=fm.solve(sample,1,starts=[p,np.array([p[0]*.8,p[1]*1.2]),np.array([p[0]*1.2,max(p[1]*.8,1e-6)])])
                if best is None:raise ValueError('No converged bootstrap start')
                records.append(dict(draw=draw,success=True,p=best['p'].tolist(),q2=best['q2']))
            except (ValueError,np.linalg.LinAlgError) as exc:records.append(dict(draw=draw,success=False,p=None,q2=None,error=str(exc)))
        solutions[row['fit_id']]['bootstrap']=records
        valid=[r['q2'] for r in records if r['success']];complete=len(valid)==399
        exceeds=sum(q>=row['q2'] for q in valid)
        value=(exceeds+1)/400 if complete else None
        bootstrap.append(dict(fit_id=row['fit_id'],Nv=row['Nv'],R=row['R'],tau_min=8,tau_max=32,draws=399,
            successful_draws=len(valid),p_recentered_bootstrap=value,
            monte_carlo_se=math.sqrt(value*(1-value)/400) if complete else None,
            p_gaussian_reference=row['p_gaussian_reference'],p_hotelling_linearized_reference=row['p_hotelling_linearized_reference']))
        print('BOOTSTRAP',row['fit_id'],len(valid),'valid; p',value,flush=True)
    write_csv(RESULTS/'fits_71configs.csv',rows);write_json(RESULTS/'solutions.json',solutions)
    write_csv(RESULTS/'bootstrap_goodness_71configs.csv',bootstrap)
    comparison=[];lookup={(r['Nv'],r['R'],r['tau_min'],r['tau_max'],r['states']):r for r in rows}
    for r in rows:
        if r['states']!=2:continue
        one=lookup[r['Nv'],r['R'],r['tau_min'],r['tau_max'],1]
        comparison.append(dict(Nv=r['Nv'],R=r['R'],tau_min=r['tau_min'],tau_max=r['tau_max'],
            single_E0=one['E0'],two_E0=r['E0'],q2_single=one['q2'],q2_two=r['q2'],
            delta_q2=one['q2']-r['q2'],two_state_identified=r['identified'],two_state_flags=r['flags']))
    write_csv(RESULTS/'model_comparison_71configs.csv',comparison)
    summary=dict(success=True,configuration_count=71,configuration_labels=CONFIGS,
        one_state_fits=sum(r['states']==1 for r in rows),two_state_fits=sum(r['states']==2 for r in rows),
        identified_two_state_fits=sum(r['states']==2 and r['identified'] for r in rows),
        nonlinear_jackknife_requested=sum(r['jackknife_requested'] for r in rows),nonlinear_jackknife_performed=sum(r['jackknife_performed'] for r in rows),
        nonlinear_jackknife_all_replicates_resolved=sum(r['jackknife_reliable'] for r in rows),
        bootstrap_null_tests=bootstrap,seconds=time.monotonic()-start,
        versions=dict(numpy=np.__version__,scipy=scipy.__version__),
        interpretation='Forward positive-amplitude one-/two-state models; bulk approximation with open time boundaries. Signed data retained. No automatic state, window or potential selection.',
        p_value_note='Gaussian and linearized Hotelling references are diagnostic. Recentered whole-configuration bootstrap re-estimates covariance in each draw; finite-draw uncertainty is separate from finite-sample bias. Delta q2 is not assigned a chi-square p value for the boundary-nested two-state model.')
    write_json(RESULTS/'SUMMARY.json',summary)
    print('PASS correlated exponential fits:',summary['one_state_fits'],'single,',summary['two_state_fits'],'double;',summary['identified_two_state_fits'],'identified two-state diagnostics',flush=True)


if __name__=='__main__':main()
