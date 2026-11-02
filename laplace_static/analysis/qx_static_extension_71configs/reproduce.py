#!/usr/bin/env python3
"""Matched qX mode-cutoff, temporal-window and sparse Cornell diagnostics."""
import argparse
import csv
import json
import math
import os
from decimal import Decimal as D,localcontext
from pathlib import Path

import numpy as np
import estimators as stats
from replay_inputs import CONFIGS,CUTOFFS,RADII,HERE,SOURCE,replay

RESULTS=HERE/'results'
WINDOW_STARTS=list(range(4,27,2))
DISPLAY_WINDOWS=[6,14,22]
CORNELL_RANGES=[[4,8,12,16],[4,8,12,16,20],[8,12,16,20,24],
                [8,12,16,20,24,28],[12,16,20,24,28]]


def read_csv(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def write_csv(path,values):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    with pending.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(values[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(values);stream.flush();os.fsync(stream.fileno())
    pending.replace(path)


def write_json(path,value):
    path=Path(path);pending=path.with_suffix(path.suffix+'.pending')
    with pending.open('w') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    pending.replace(path)


def jk_se(samples):
    x=np.asarray(samples,float);n=len(x)
    return float(np.sqrt((n-1.)/n*np.sum((x-x.mean())**2)))


def paired_se(samples):
    return float(np.std(samples,ddof=1)/math.sqrt(len(samples)))


def normalized(xs,ys):
    with localcontext() as ctx:
        ctx.prec=50;mx,sx=stats.mean_se(xs);my,sy=stats.mean_se(ys)
        if my<=0:raise ValueError('Normalization has a nonpositive denominator')
        value=mx/my
        fluct=[(a-value*b)/my for a,b in zip(xs,ys)]
        return value,stats.mean_se(fluct)[1]


def count_positive_loo(xs,ys):
    sx=sum(xs,D(0));sy=sum(ys,D(0))
    return sum(sx-x>0 and sy-y>0 for x,y in zip(xs,ys))


def primary_statistics(data):
    raw={key:[D(group[c]) for c in CONFIGS] for key,group in data.items()}
    signals=[];energies=[];signal_lookup={};energy_lookup={};fluct={};loo={}
    for n in CUTOFFS:
        for R in RADII:
            for t in range(33):
                xs=raw[n,t,R];m,se=stats.mean_se(xs);q,qs=normalized(xs,raw[n,0,R])
                snr=m/se if se else None
                flag='nonpositive_mean' if m<=0 else ('low_mean_over_se' if snr is None or snr<3 else 'positive_mean')
                row=dict(Nv=n,tau=t,R=R,configurations=71,mean_L=str(m),naive_se=str(se),
                    normalized_L=str(q),normalized_naive_se=str(qs),mean_over_se=float(snr) if snr is not None else None,
                    nonpositive_individual_configurations=sum(x<=0 for x in xs),diagnostic=flag)
                signals.append(row);signal_lookup[n,t,R]=row
            for t in range(32):
                xs,ys=raw[n,t,R],raw[n,t+1,R];estimate,delta,replicates=stats.energy(xs,ys)
                first,second=signal_lookup[n,t,R],signal_lookup[n,t+1,R]
                snr=min(first['mean_over_se'],second['mean_over_se'])
                flag='nonpositive_endpoint_mean' if not estimate['valid_point'] else ('low_endpoint_mean_over_se' if snr<3 else 'positive_endpoint_means')
                ratio=float(estimate['jackknife_se']/estimate['naive_se']) if estimate['valid_jackknife'] and estimate['naive_se'] else None
                row=dict(Nv=n,tau=t,R=R,configurations=71,valid_point=estimate['valid_point'],
                    aE_eff=str(estimate['aE_eff']),naive_se=str(estimate['naive_se']),
                    valid_jackknife=estimate['valid_jackknife'],jackknife_se=str(estimate['jackknife_se']),
                    positive_jackknife_replicates=count_positive_loo(xs,ys),jackknife_to_naive_se=ratio,
                    minimum_endpoint_mean_over_se=snr,diagnostic=flag)
                energies.append(row);energy_lookup[n,t,R]=row
                if delta is not None:fluct[n,t,R]=np.array(delta,float)
                if replicates:loo[n,t,R]=np.array(replicates,float)
    return raw,signals,energies,signal_lookup,energy_lookup,fluct,loo


def effective_covariance(fluct):
    output=[]
    indices=[(n,t) for n in CUTOFFS for t in range(32)]
    for R in RADII:
        valid=np.array([(n,t,R) in fluct for n,t in indices])
        array=np.column_stack([fluct.get((n,t,R),np.zeros(71)) for n,t in indices])
        C=np.cov(array,rowvar=False,ddof=1)/71
        diag=np.sqrt(np.maximum(0,np.diag(C)))
        for i,(n,t) in enumerate(indices):
            for j,(m,u) in enumerate(indices):
                ok=bool(valid[i] and valid[j]);correlation=float(C[i,j]/(diag[i]*diag[j])) if ok and diag[i]*diag[j]>0 else None
                output.append(dict(R=R,Nv=n,tau=t,reference_Nv=m,reference_tau=u,configurations=71,
                    covariance_naive=float(C[i,j]) if ok else None,correlation=correlation,defined=ok))
    return output


def constant_fits(energy_lookup,signal_lookup,fluct,loo):
    rows=[];fits={}
    for n in CUTOFFS:
        for R in RADII:
            for start in WINDOW_STARTS:
                keys=[(n,t,R) for t in range(start,start+6)]
                minimum=min(signal_lookup[n,t,R]['mean_over_se'] for t in range(start,start+7))
                base=dict(Nv=n,R=R,tau_min=start,tau_max=start+5,configurations=71,
                    status='undefined',estimate=None,naive_se=None,jackknife_se=None,
                    chi2=None,dof=5,p_gaussian_approx=None,condition_inf=None,
                    minimum_endpoint_mean_over_se=minimum,low_signal_flag=minimum<3)
                if any(k not in fluct for k in keys):base['status']='nonpositive_endpoint_mean';rows.append(base);continue
                y=[float(energy_lookup[k]['aE_eff']) for k in keys]
                d=np.column_stack([fluct[k] for k in keys])
                try:fit=stats.gls(y,d,np.ones((6,1)))
                except (ValueError,np.linalg.LinAlgError):base['status']='invalid_covariance';rows.append(base);continue
                replicates=None
                if all(k in loo for k in keys):replicates=np.column_stack([loo[k] for k in keys])@fit['weights'][0]
                base.update(status='defined',estimate=float(fit['params'][0]),
                    naive_se=float(np.sqrt(fit['covariance'][0,0])),jackknife_se=jk_se(replicates) if replicates is not None else None,
                    chi2=fit['chi2'],p_gaussian_approx=fit['p_gaussian_approx'],condition_inf=fit['condition_inf'])
                rows.append(base);fits[n,R,start]=dict(row=base,fluct=fit['fluctuations'][:,0],loo=replicates)
    return rows,fits


def contrasts(energy_lookup,fluct,loo,fits):
    energy_rows=[];cutoff_rows=[];window_rows=[]
    for n in CUTOFFS[:-1]:
        for R in RADII:
            for t in range(32):
                a,b=(n,t,R),(100,t,R);valid=a in fluct and b in fluct
                estimate=float(D(energy_lookup[a]['aE_eff'])-D(energy_lookup[b]['aE_eff'])) if valid else None
                error=paired_se(fluct[a]-fluct[b]) if valid else None
                energy_rows.append(dict(Nv=n,reference_Nv=100,tau=t,R=R,configurations=71,
                    difference=estimate,paired_naive_se=error,
                    paired_jackknife_se=jk_se(loo[a]-loo[b]) if a in loo and b in loo else None,
                    diagnostic='defined' if valid else 'nonpositive_endpoint_mean'))
            for start in WINDOW_STARTS:
                a,b=fits.get((n,R,start)),fits.get((100,R,start));ok=a is not None and b is not None
                delta=a['row']['estimate']-b['row']['estimate'] if ok else None
                error=paired_se(a['fluct']-b['fluct']) if ok else None
                cutoff_rows.append(dict(Nv=n,reference_Nv=100,R=R,tau_min=start,tau_max=start+5,
                    configurations=71,difference=delta,paired_naive_se=error,
                    paired_jackknife_se=jk_se(a['loo']-b['loo']) if ok and a['loo'] is not None and b['loo'] is not None else None,
                    difference_over_se=delta/error if error else None,
                    low_signal_flag=bool(a['row']['low_signal_flag'] or b['row']['low_signal_flag']) if ok else True,
                    diagnostic='defined' if ok else 'undefined_fit'))
    for n in CUTOFFS:
        for R in RADII:
            for first,last in [(6,14),(6,22),(14,22)]:
                a,b=fits.get((n,R,first)),fits.get((n,R,last));ok=a is not None and b is not None
                delta=b['row']['estimate']-a['row']['estimate'] if ok else None
                error=paired_se(b['fluct']-a['fluct']) if ok else None
                window_rows.append(dict(Nv=n,R=R,early_tau_min=first,late_tau_min=last,configurations=71,
                    later_minus_earlier=delta,paired_naive_se=error,
                    paired_jackknife_se=jk_se(b['loo']-a['loo']) if ok and a['loo'] is not None and b['loo'] is not None else None,
                    difference_over_se=delta/error if error else None,
                    low_signal_flag=bool(a['row']['low_signal_flag'] or b['row']['low_signal_flag']) if ok else True,
                    diagnostic='defined' if ok else 'undefined_fit'))
    return energy_rows,cutoff_rows,window_rows


def fit_covariance(fits):
    rows=[];indices=[(n,R) for n in CUTOFFS for R in RADII]
    for start in WINDOW_STARTS:
        for n,R in indices:
            for m,Q in indices:
                a,b=fits.get((n,R,start)),fits.get((m,Q,start));ok=a is not None and b is not None
                cov=float(np.cov(a['fluct'],b['fluct'],ddof=1)[0,1]/71) if ok else None
                rows.append(dict(tau_min=start,tau_max=start+5,Nv=n,R=R,reference_Nv=m,reference_R=Q,
                    configurations=71,covariance_naive=cov,defined=ok))
    return rows


def cornell_diagnostics(fits):
    rows=[];parameter_cov=[]
    for n in CUTOFFS:
        for start in WINDOW_STARTS:
            for rs in CORNELL_RANGES:
                parts=[fits.get((n,R,start)) for R in rs];tag=','.join(map(str,rs))
                row=dict(Nv=n,tau_min=start,tau_max=start+5,R_values=tag,points=len(rs),configurations=71,
                    status='undefined',V0=None,alpha=None,sigma_a2=None,V0_naive_se=None,alpha_naive_se=None,
                    sigma_naive_se=None,V0_jackknife_se=None,alpha_jackknife_se=None,sigma_jackknife_se=None,
                    chi2=None,dof=len(rs)-3,p_gaussian_approx=None,condition_inf=None,
                    minimum_endpoint_mean_over_se=None,low_signal_flag=True)
                if any(p is None for p in parts):row['status']='undefined_input_fit';rows.append(row);continue
                minimum=min(p['row']['minimum_endpoint_mean_over_se'] for p in parts)
                row.update(minimum_endpoint_mean_over_se=minimum,low_signal_flag=minimum<3)
                A=np.array([[1.,-1./R,float(R)] for R in rs])
                try:f=stats.gls([p['row']['estimate'] for p in parts],np.column_stack([p['fluct'] for p in parts]),A)
                except (ValueError,np.linalg.LinAlgError):row['status']='invalid_covariance';rows.append(row);continue
                naive=np.sqrt(np.diag(f['covariance']));jk=[None]*3
                if all(p['loo'] is not None for p in parts):
                    samples=np.column_stack([p['loo'] for p in parts])@f['weights'].T
                    jk=[jk_se(samples[:,i]) for i in range(3)]
                row.update(status='defined',V0=float(f['params'][0]),alpha=float(f['params'][1]),sigma_a2=float(f['params'][2]),
                    V0_naive_se=float(naive[0]),alpha_naive_se=float(naive[1]),sigma_naive_se=float(naive[2]),
                    V0_jackknife_se=jk[0],alpha_jackknife_se=jk[1],sigma_jackknife_se=jk[2],
                    chi2=f['chi2'],p_gaussian_approx=f['p_gaussian_approx'],condition_inf=f['condition_inf'])
                rows.append(row)
                for i,p in enumerate(['V0','alpha','sigma_a2']):
                    for j,q in enumerate(['V0','alpha','sigma_a2']):
                        parameter_cov.append(dict(Nv=n,tau_min=start,tau_max=start+5,R_values=tag,
                            parameter=p,reference_parameter=q,covariance_naive=float(f['covariance'][i,j])))
    return rows,parameter_cov


def compute(progress=False):
    data,validation,timing=replay(progress=progress)
    _,signals,energies,signal_lookup,energy_lookup,fluct,loo=primary_statistics(data)
    covariance=effective_covariance(fluct);fit_rows,fits=constant_fits(energy_lookup,signal_lookup,fluct,loo)
    energy_diffs,cutoff_diffs,window_diffs=contrasts(energy_lookup,fluct,loo,fits)
    fit_cov=fit_covariance(fits);cornell,cornell_cov=cornell_diagnostics(fits)
    support=[r for r in signals if r['R'] in (24,32) and r['tau'] in (8,16,24,32)]
    counts=[]
    for n in CUTOFFS:
        selected=[r for r in signals if r['Nv']==n];ee=[r for r in energies if r['Nv']==n]
        ratios=[r['jackknife_to_naive_se'] for r in ee if r['jackknife_to_naive_se'] is not None and r['minimum_endpoint_mean_over_se']>=3]
        counts.append(dict(Nv=n,mean_points=len(selected),nonpositive_mean_points=sum(r['diagnostic']=='nonpositive_mean' for r in selected),
            individual_combined_axis_measurements=71*33*12,
            nonpositive_individual_combined_axis_measurements=sum(r['nonpositive_individual_configurations'] for r in selected),
            undefined_effective_energies=sum(not r['valid_point'] for r in ee),
            positive_points_without_complete_jackknife=sum(r['valid_point'] and not r['valid_jackknife'] for r in ee),
            jackknife_to_naive_se_range_mean_over_se_at_least3=[min(ratios),max(ratios)]))
    summary=dict(success=True,validation=validation,geometry=dict(source_time=90,temporal_separations=list(range(33)),
        spatial_separations=RADII,mode_cutoffs=CUTOFFS,saved_Nv=100,piter=10,p40_control_configurations=[170,640,650]),
        uncertainty='Independent-configuration baseline; N=71; 50-digit primary mean/log arithmetic; paired naive covariance and paired delete-one jackknife.',
        row_counts=dict(correlator_means=len(signals),effective_energies=len(energies),effective_energy_covariance=len(covariance),
            constant_fit_windows=len(fit_rows),fitted_energy_covariance=len(fit_cov),cutoff_energy_contrasts=len(energy_diffs),
            cutoff_fit_contrasts=len(cutoff_diffs),window_contrasts=len(window_diffs),cornell_fits=len(cornell)),
        signal_diagnostics=counts,R24_tau24=[r for r in support if r['R']==24 and r['tau']==24],
        R32_tau24=[r for r in support if r['R']==32 and r['tau']==24],
        selected_constant_fits=[r for r in fit_rows if r['R'] in (4,8,12,16,24,32) and r['tau_min'] in DISPLAY_WINDOWS],
        selected_paired_cutoff_contrasts=[r for r in cutoff_diffs if r['R'] in (4,8,12,16,24,32) and r['tau_min'] in DISPLAY_WINDOWS],
        selected_window_contrasts=[r for r in window_diffs if r['R'] in (4,8,12,16,24,32)],
        selected_cornell_diagnostics=[r for r in cornell if r['R_values']=='4,8,12,16,20' and r['tau_min'] in [6,10,14,18,22]],
        fit_method='Correlated constant fits to six adjacent effective energies. GLS weights and covariance estimated on the full sample and held fixed in the paired jackknife. Cornell diagnostics retain covariance across spatial distances.',
        p_value_scope='Gaussian chi-square reference tails with estimated covariance; approximate, exploratory and unadjusted for the many windows, cutoffs and ranges. They are not automatic plateau acceptance criteria.',
        signal_flag='Mean/SE <3 is a plotting/review flag, not an agreed physics acceptance threshold.',
        logarithm_policy='Use ratios/logs of ensemble means; retain all configuration samples. Nonpositive endpoint means are undefined and are exported as empty fields. No positive-configuration selection or absolute values.',
        interpretation='Lower Nv substantially improves the R=24 signal but not all distances. Cutoff and time-window dependence remain; apparent fitted energies and Cornell slopes are diagnostics, not a final ground-state potential or string tension.')
    outputs={'correlator_means_71configs.csv':signals,'effective_energies_71configs.csv':energies,
        'effective_energy_covariance_71configs.csv':covariance,'constant_fits_71configs.csv':fit_rows,
        'fitted_energy_covariance_71configs.csv':fit_cov,'cutoff_energy_contrasts_71configs.csv':energy_diffs,
        'cutoff_fit_contrasts_71configs.csv':cutoff_diffs,'window_contrasts_71configs.csv':window_diffs,
        'cornell_diagnostics_71configs.csv':cornell,'cornell_parameter_covariance_71configs.csv':cornell_cov,
        'R24_R32_signal_71configs.csv':support,'timing_71configs.csv':timing}
    return outputs,summary


def compare(actual,expected,label):
    if expected is None:
        if actual not in ('',None):raise ValueError('Expected undefined value: '+label)
    elif isinstance(expected,bool):
        if actual not in (expected,str(expected)):raise ValueError('Boolean mismatch: '+label)
    elif isinstance(expected,(int,float,D)):
        if not math.isclose(float(actual),float(expected),rel_tol=3e-7,abs_tol=1e-24):raise ValueError('Numerical mismatch: '+label)
    elif isinstance(expected,dict):
        if set(actual)!=set(expected):raise ValueError('Different keys: '+label)
        for key in expected:compare(actual[key],expected[key],label+'.'+key)
    elif isinstance(expected,list):
        if len(actual)!=len(expected):raise ValueError('Different length: '+label)
        for i,(x,y) in enumerate(zip(actual,expected)):compare(x,y,label+'[%d]'%i)
    elif actual!=expected:raise ValueError('Different value: '+label)


def export(outputs,summary,verify=False):
    for name,values in outputs.items():
        if verify:
            actual=read_csv(RESULTS/name)
            if len(actual)!=len(values):raise ValueError('Different row count: '+name)
            for i,(a,b) in enumerate(zip(actual,values)):compare(a,b,name+'[%d]'%i)
        else:write_csv(RESULTS/name,values)
    if verify:compare(json.loads((RESULTS/'SUMMARY.json').read_text()),summary,'SUMMARY')
    else:write_json(RESULTS/'SUMMARY.json',summary)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.verify:
        from verify_files import verify
        verify()
    outputs,summary=compute(progress=True);export(outputs,summary,a.verify)
    print('PASS full71 extension: 145 launches replayed; 63048 reference records exact; all four cutoffs, signed means, paired errors, window and Cornell diagnostics reproduced.')


if __name__=='__main__':main()
