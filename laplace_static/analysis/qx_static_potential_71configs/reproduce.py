#!/usr/bin/env python3
"""Replay and analyze all 71 qX static-correlator configurations."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import sys
import tempfile
from decimal import Decimal as D,localcontext
from pathlib import Path
import estimators as stat

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];INPUTS=HERE/'inputs';RESULTS=HERE/'results'
CONFIGS=list(range(100,801,10));AXES=(-1,0,1,2)
EXACT=('correlators_ensemble_71configs.csv','effective_energies_71configs.csv','temporal_covariance_71configs.csv')
FLOAT=('plateau_fits_71configs.csv','window_stability_71configs.csv','potential_covariance_71configs.csv','cornell_fits_71configs.csv','spatial_block_diagnostics_71configs.csv')


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,default=str,allow_nan=False)+'\n')


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def load():
    spec=json.loads((INPUTS/'SOURCES.json').read_text())
    for name,sha in spec['files'].items():
        if digest(HERE/name)!=sha:raise ValueError('Input changed: '+name)
    if digest(HERE/'ANALYSIS_PLAN.json')!=spec['analysis_plan_sha256']:raise ValueError('Analysis settings changed')
    validation=ROOT/'laplace_static/validation/qx_static_potential'
    for r in spec['archive_mapping']:
        if r['preserved_path'].startswith('inputs/workflow/') and digest(ROOT/r['archive_path'])!=r['sha256']:raise ValueError('Executed source differs from committed source')
    sys.path.insert(0,str(validation))
    try:c=module('qx_static71_controller',validation/'campaign.py')
    finally:sys.path.pop(0)
    c.verify_sources()
    full=INPUTS/'job5955993';pilot=INPUTS/'job5955349';plan=json.loads((full/'PREPARED.json').read_text());pp=json.loads((pilot/'PREPARED.json').read_text())
    for target in (full,pilot):
        if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Plan checksum changed')
    if plan['build_manifest_sha256']!=digest(INPUTS/'build_provenance/MANIFEST.json') or plan['pilot']['prepared_sha256']!=digest(pilot/'PREPARED.json'):raise ValueError('Build/reuse association changed')
    refs=c.references();raw=[];blocks=[];summaries=[];launches=[]
    for entry in plan['cases']:
        out=c.complete_case(full,entry)
        if out is None:raise ValueError('Missing production case')
        a,b,d=c.replay_case(out,entry,refs)
        if a!=read_csv(out/'correlators.csv') or b!=read_csv(out/'spatial_blocks.csv'):raise ValueError('Case table differs from stdout')
        raw+=a;blocks+=b;summaries.append(d);launches+=json.loads((out/'LAUNCHES.json').read_text())
    for reuse in plan['pilot']['cases']:
        entry=next(e for e in pp['cases'] if e['configuration']==reuse['configuration'])
        parts=Path(reuse['attempt']).parts;idx=parts.index('campaigns');suffix=parts[idx+2:]
        if len(suffix)!=2 or not re.fullmatch(r'case_\d{3}',suffix[0]) or not re.fullmatch(r'attempt_\d{4}',suffix[1]):raise ValueError('Invalid pilot reuse path')
        out=pilot.joinpath(*suffix);case=json.loads((out/'CASE.json').read_text());summary=json.loads((out/'SUMMARY.json').read_text())
        if case!=dict(entry,prepared_sha256=digest(pilot/'PREPARED.json')) or not summary['success'] or digest(out/'LAUNCHES.json')!=summary['launches_sha256'] or digest(out/'CASE.json')!=summary['case_sha256']:raise ValueError('Pilot provenance changed')
        a,b,d=c.replay_case(out,entry,refs);raw += [r for r in a if r['piter']=='10'];blocks += [r for r in b if r['piter']=='10'];summaries.append(d);launches+=json.loads((out/'LAUNCHES.json').read_text())
    archived=json.loads((full/'SUMMARY.json').read_text())
    if summaries!=archived['cases'] or sorted(s['configuration'] for s in summaries)!=CONFIGS or not archived['success'] or archived['missing_configurations']:raise ValueError('Incomplete cohort')
    if raw!=read_csv(full/'correlators.csv') or blocks!=read_csv(full/'spatial_blocks.csv') or len(raw)!=157052 or len(blocks)!=272640:raise ValueError('Assembled table differs from complete log replay')
    return raw,blocks,summaries,launches


def primary_statistics(raw):
    groups={}
    for r in raw:
        if r['scheme']=='fixed':groups.setdefault(tuple(int(r[k]) for k in ('axis','tau','R')),{})[int(r['cfg'])]=(D(r['L_RE']),D(r['L_IM']))
    if len(groups)!=2176 or any(sorted(g)!=CONFIGS for g in groups.values()):raise ValueError('Unmatched production geometry')
    means=[];energies=[];deltas={};loo={};points={};covrows=[]
    with localcontext() as ctx:
        ctx.prec=50
        for (axis,tau,R),g in sorted(groups.items()):
            x=[g[c][0] for c in CONFIGS];y=[g[c][1] for c in CONFIGS];m,e=stat.mean_se(x);im,ie=stat.mean_se(y)
            means.append(dict(axis=axis,t0=90,tau=tau,R=R,piter=10,configurations=71,L_RE=m,L_RE_naive_se=e,L_IM=im,L_IM_naive_se=ie,
                mean_over_se=m/e if e else '',positive_mean=m>0,nonpositive_individual=sum(v<=0 for v in x)))
        for axis in AXES:
            for R in range(1,33):
                for tau in range(16):
                    x=[groups[axis,tau,R][c][0] for c in CONFIGS];y=[groups[axis,tau+1,R][c][0] for c in CONFIGS]
                    r,d,jk=stat.energy(x,y);row=dict(axis=axis,t0=90,tau=tau,tau_next=tau+1,R=R,piter=10,configurations=71);row.update(r)
                    energies.append(row);points[axis,R,tau]=row;deltas[axis,R,tau]=d;loo[axis,R,tau]=jk
                for t in range(16):
                    for u in range(16):
                        a,b=deltas[axis,R,t],deltas[axis,R,u];v=stat.covariance(a,b) if a is not None and b is not None else ''
                        den=points[axis,R,t]['naive_se']*points[axis,R,u]['naive_se'] if v!='' else 0
                        covrows.append(dict(axis=axis,R=R,tau_i=t,tau_j=u,covariance=v,correlation=v/den if den else ''))
    return means,energies,covrows,points,deltas,loo


def fits(means,points,deltas,loo,settings):
    import numpy as np
    mean_lookup={(r['axis'],r['R'],r['tau']):r for r in means};windows=[tuple(w) for w in settings['plateau_windows']]
    records=[];saved={};stability=[];potential_cov=[];cornell=[]
    for axis in AXES:
        for R in range(1,33):
            for lo,hi in windows:
                row=dict(axis=axis,R=R,tau_min=lo,tau_max=hi,points=hi-lo+1,aV='',naive_se='',jackknife_se='',chi2='',dof=hi-lo,p_gaussian_approx='',condition_inf='',minimum_mean_over_se='',diagnostic_compatible=False,status='invalid_mean')
                taus=list(range(lo,hi+1));rows=[points[axis,R,t] for t in taus]
                if all(r['valid_point'] for r in rows):
                    data=np.array([[float(v) for v in deltas[axis,R,t]] for t in taus]).T
                    fit=stat.gls([float(r['aE_eff']) for r in rows],data,np.ones((len(taus),1)))
                    v=float(fit['params'][0]);err=float(np.sqrt(fit['covariance'][0,0]));snr=min(float(mean_lookup[axis,R,t]['mean_over_se']) for t in range(lo,hi+2))
                    row.update(aV=v,naive_se=err,chi2=fit['chi2'],p_gaussian_approx=fit['p_gaussian_approx'],condition_inf=fit['condition_inf'],minimum_mean_over_se=snr,status='fit')
                    if all(r['valid_jackknife'] for r in rows):
                        jk=np.array([[float(x) for x in loo[axis,R,t]] for t in taus]).T@fit['weights'][0]
                        row['jackknife_se']=float(np.sqrt(70/71*np.sum((jk-jk.mean())**2)))
                    row['diagnostic_compatible']=fit['p_gaussian_approx']>=settings['diagnostic_minimum_gaussian_p'] and snr>=settings['diagnostic_minimum_mean_over_se'] and fit['condition_inf']<=float(settings['condition_inf_limit'])
                    saved[axis,R,lo,hi]=dict(row=row,delta=fit['fluctuations'][:,0],weights=fit['weights'][0])
                records.append(row)
    a,b=settings['reference_window']
    for axis in AXES:
        for R in range(1,33):
            first=saved.get((axis,R,a,b))
            for low,high in settings['temporal_comparison_windows']:
                second=saved.get((axis,R,low,high));r=dict(axis=axis,R=R,reference_min=a,reference_max=b,comparison_min=low,comparison_max=high,delta_aV='',paired_se='',z='')
                if first and second:
                    d=first['delta']-second['delta'];se=float(d.std(ddof=1)/np.sqrt(71));value=first['row']['aV']-second['row']['aV']
                    r.update(delta_aV=value,paired_se=se,z=value/se if se else '')
                stability.append(r)
        for R in range(1,33):
            for S in range(1,33):
                first=saved.get((axis,R,a,b));second=saved.get((axis,S,a,b));v=float(np.cov(first['delta'],second['delta'],ddof=1)[0,1]/71) if first and second else ''
                potential_cov.append(dict(axis=axis,R_i=R,R_j=S,tau_min=a,tau_max=b,covariance=v))
    for low,high in settings['cornell_windows']:
        for rmin,rmax in settings['cornell_ranges']:
            selected=[saved.get((-1,R,low,high)) for R in range(rmin,rmax+1)]
            row=dict(tau_min=low,tau_max=high,R_min=rmin,R_max=rmax,V0='',sigma_a2='',alpha='',V0_se='',sigma_a2_se='',alpha_se='',chi2='',dof=rmax-rmin-2,p_gaussian_approx='',condition_inf='',all_temporal_windows_compatible=False,status='invalid_mean')
            if all(selected):
                R=np.arange(rmin,rmax+1,dtype=float);design=np.column_stack((np.ones(len(R)),R,-1/R))
                fit=stat.gls([s['row']['aV'] for s in selected],np.array([s['delta'] for s in selected]).T,design)
                e=np.sqrt(np.diag(fit['covariance']));p=fit['params']
                row.update(V0=float(p[0]),sigma_a2=float(p[1]),alpha=float(p[2]),V0_se=float(e[0]),sigma_a2_se=float(e[1]),alpha_se=float(e[2]),chi2=fit['chi2'],p_gaussian_approx=fit['p_gaussian_approx'],condition_inf=fit['condition_inf'],all_temporal_windows_compatible=all(s['row']['diagnostic_compatible'] for s in selected),status='fit')
            cornell.append(row)
    return records,stability,potential_cov,cornell


def block_diagnostics(blocks):
    import numpy as np
    groups={}
    for r in blocks:groups.setdefault(tuple(int(r[k]) for k in ('tau','R')),{}).setdefault(int(r['axis']),{}).setdefault(int(r['cfg']),{})[tuple(int(r[k]) for k in ('px','py','pz'))]=float(r['L_RE'])
    out=[]
    for (tau,R),group in sorted(groups.items()):
        positions=[(x,y,z) for x in range(4) for y in range(4) for z in range(4)]
        arrays={a:np.array([[group[a][c][j] for j in positions] for c in CONFIGS]) for a in (0,1,2)}
        arrays[-1]=np.concatenate([arrays[a] for a in (0,1,2)],axis=1)
        for axis in AXES:
            data=arrays[axis];N,M=data.shape;averages=data.mean(axis=1);correct=float(averages.std(ddof=1)/np.sqrt(N))
            diagonal=float(data.var(axis=0,ddof=1).sum()/(N*M*M));full=correct**2
            out.append(dict(axis=axis,tau=tau,R=R,configurations=N,regional_channels=M,mean=float(averages.mean()),configuration_naive_se=correct,ignore_regional_covariance_se=float(np.sqrt(diagonal)),correct_over_ignore_regional_covariance=correct/math.sqrt(diagonal),off_diagonal_variance_contribution=full-diagonal))
    return out


def numerical_exports(output=RESULTS):
    raw,blocks,validation,launches=load();settings=json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    means,energies,covrows,points,deltas,loo=primary_statistics(raw)
    fitrows,stability,potential_cov,cornell=fits(means,points,deltas,loo,settings);spatial=block_diagnostics(blocks)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for name,rows in zip(EXACT+FLOAT,(means,energies,covrows,fitrows,stability,potential_cov,cornell,spatial)):write_csv(output/name,rows)
    pooled=[r for r in means if r['axis']==-1];ref=[r for r in fitrows if r['axis']==-1 and [r['tau_min'],r['tau_max']]==settings['reference_window']]
    late=[r for r in stability if r['axis']==-1 and [r['comparison_min'],r['comparison_max']]==[8,13] and r['z']!='']
    summary=dict(success=True,configuration_count=71,configuration_labels=CONFIGS,piter=10,geometry=dict(t0=90,taus=list(range(17)),R=list(range(1,33)),Nv=100,rho=1),
        raw_correlator_rows=157052,spatial_block_rows=272640,audited_modes=120700,mpi_launches_replayed=145,maximum_relative_mode_residual=max(s['quality']['maximum_relative_residual'] for s in validation),
        maximum_p10_unitarity_defect=max(float(v['unit_max']) for s in validation for v in s['scans'] if v['piter']==10),
        production_mean_rows=len(means),effective_energy_rows=len(energies),temporal_covariance_rows=len(covrows),plateau_fit_rows=len(fitrows),cornell_fit_rows=len(cornell),
        uncertainty='Independent configurations, N=71, with full covariance of ratio/log estimates and fixed-weight paired delete-one controls. Spatial sources and axes are not extra independent samples.',
        pooled_nonpositive_mean_points=sum(not r['positive_mean'] for r in pooled),pooled_nonpositive_individual_measurements=sum(r['nonpositive_individual'] for r in pooled),minimum_pooled_mean_over_se=min(r['mean_over_se'] for r in pooled),
        reference_window=settings['reference_window'],reference_window_compatible_R=[r['R'] for r in ref if r['diagnostic_compatible']],
        late_window_shift_above_two_paired_se_R=[r['R'] for r in late if abs(r['z'])>2],reference_potential=ref,cornell_fits=cornell,
        interpretation='Measured correlators and effective energies cover the full ensemble. Window compatibility is exploratory; later-time drift and fit-range dependence prevent treating every reference-window value as a final ground-state potential. Cornell parameters are diagnostic candidates, not a finalized string tension.',
        p_value_scope=settings['fit_interpretation'])
    write_json(output/'SUMMARY.json',summary)
    return means,energies,fitrows,stability,cornell,spatial,summary


def compare_numeric_csv(a,b):
    aa,bb=read_csv(a),read_csv(b)
    if len(aa)!=len(bb):raise ValueError('Changed fit row count')
    for x,y in zip(aa,bb):
        if list(x)!=list(y):raise ValueError('Changed fit columns')
        for k in x:
            if x[k]==y[k]:continue
            try:good=math.isclose(float(x[k]),float(y[k]),rel_tol=3e-9,abs_tol=1e-14)
            except ValueError:good=False
            if not good:raise ValueError('Fit differs: '+k)


def compare_json(a,b):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise ValueError('Summary keys differ')
        for k in a:compare_json(a[k],b[k])
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise ValueError('Summary list differs')
        for x,y in zip(a,b):compare_json(x,y)
    elif isinstance(a,float):
        if not math.isclose(a,b,rel_tol=3e-9,abs_tol=1e-14):raise ValueError('Summary fit value differs')
    elif a!=b:raise ValueError('Summary value differs')


def verify():
    with tempfile.TemporaryDirectory() as tmp:
        numerical_exports(tmp)
        for name in EXACT:
            if (Path(tmp)/name).read_bytes()!=(RESULTS/name).read_bytes():raise ValueError('Exact numerical export differs: '+name)
        for name in FLOAT:compare_numeric_csv(Path(tmp)/name,RESULTS/name)
        compare_json(json.loads((Path(tmp)/'SUMMARY.json').read_text()),json.loads((RESULTS/'SUMMARY.json').read_text()))
    checked=read_csv(RESULTS/'pyerrors_S0_check_71configs.csv')
    baseline={(r['tau'],r['R']):r for r in read_csv(RESULTS/EXACT[1]) if r['axis']=='-1'}
    if len(checked)!=512 or {(r['tau'],r['R']) for r in checked}!=set(baseline):raise ValueError('Incomplete PyErrors check')
    for r in checked:
        ref=baseline[r['tau'],r['R']]
        if r['analytic_se']!=ref['naive_se'] or not math.isclose(float(r['pyerrors_S0_se']),float(ref['naive_se']),rel_tol=2e-9,abs_tol=2e-15):raise ValueError('Archived PyErrors error differs')
    py=json.loads((RESULTS/'PYERRORS_CHECK.json').read_text())
    if not py['success'] or py['checks']!=512 or py['configuration_labels']!=CONFIGS or py['S']!=0:raise ValueError('PyErrors sample changed')
    print('PASS full71 static analysis: all 71 cases and 145 launches replayed; exact primary/log/covariance exports; correlated window/range fits and regional diagnostics reproduced.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');p.add_argument('--plots',action='store_true');args=p.parse_args()
    if args.verify:verify()
    else:
        data=numerical_exports()
        if args.plots:module('qx_static71_figures',HERE/'render_plots.py').render(*data[:6],output=RESULTS)
        print('PASS full71 static numerical exports')


if __name__=='__main__':main()
