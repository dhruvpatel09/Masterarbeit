"""Distance-by-distance plateau catalogue and paired stability diagnostics."""
import itertools
import json

import numpy as np
from data import HERE, ROOT, SHORT, LONG, FOCUS, LONG_RADII, CUTOFFS, read_csv
from statistics_tools import (assert_close, constant_fit, jackknife_error,
                              standard_error)


def windows(last):
    plan = json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    return sorted(set(tuple(w) for w in plan['constant_windows'] if w[1] <= last))


def fit_catalogue(data):
    rows, saved = [], {}
    series = sorted(set(k[:2] for k in data['energies']))
    for n,r in series:
        last = max(k[2] for k in data['energies'] if k[:2]==(n,r))
        for lo,hi in windows(last):
            keys = [(n,r,t) for t in range(lo,hi+1)]
            minimum = min(data['means'][n,r,t]['mean_over_se'] for t in range(lo,hi+2))
            row = dict(Nv=n,R=r,tau_min=lo,tau_max=hi,points=len(keys),configurations=71,
                status='undefined_endpoint_mean',estimate=None,naive_se=None,jackknife_se=None,
                chi2=None,dof=len(keys)-1,p_gaussian_approx=None,condition_inf=None,
                minimum_endpoint_mean_over_se=minimum,low_signal=minimum<3,
                diagnostic_screen=False)
            if all(k in data['delta'] for k in keys):
                y = np.array([float(data['energies'][k]['aE_eff']) for k in keys])
                d = np.column_stack([data['delta'][k] for k in keys])
                try:
                    fit,w = constant_fit(y,d)
                except (ValueError,np.linalg.LinAlgError):
                    row['status'] = 'invalid_covariance'
                else:
                    j = (np.column_stack([data['loo'][k] for k in keys])@w
                         if all(k in data['loo'] for k in keys) else None)
                    row.update(fit,status='defined',jackknife_se=jackknife_error(j))
                    row['diagnostic_screen'] = bool(minimum>=3 and fit['p_gaussian_approx']>=.05 and fit['condition_inf']<=1e12)
                    saved[n,r,lo,hi] = dict(row=row,fluct=d@w,loo=j,weights=w)
            rows.append(row)
    # Regression against both previously published catalogues, not just selected plots.
    short_count = long_count = 0
    for old in read_csv(ROOT/SHORT/'results/plateau_fits_71configs.csv'):
        if old['axis'] != '-1':
            continue
        key = (100,int(old['R']),int(old['tau_min']),int(old['tau_max']))
        new = saved[key]['row']
        for a,b in (('estimate','aV'),('naive_se','naive_se'),('jackknife_se','jackknife_se'),('chi2','chi2')):
            assert_close(new[a],float(old[b]),'short constant fit '+str(key),rtol=1e-8,atol=5e-12)
        short_count += 1
    lookup = {(r['Nv'],r['R'],r['tau_min'],r['tau_max']):r for r in rows}
    for old in read_csv(ROOT/LONG/'results/constant_fits_71configs.csv'):
        key = tuple(int(old[k]) for k in ('Nv','R','tau_min','tau_max'))
        new = lookup[key]
        if (old['status']=='defined') != (new['status']=='defined'):
            raise ValueError('Published long fit status changed')
        if new['status']=='defined':
            for f in ('estimate','naive_se','jackknife_se','chi2'):
                if old[f]:
                    assert_close(new[f],float(old[f]),'long constant fit '+str(key),rtol=1e-8,atol=5e-12)
        long_count += 1
    return rows,saved,dict(short_fits_verified=short_count,long_fits_verified=long_count)


def contrast(a,b):
    if a is None or b is None:
        return dict(status='undefined_fit',difference=None,paired_naive_se=None,
                    paired_jackknife_se=None,difference_over_se=None,low_signal=True)
    difference = a['row']['estimate']-b['row']['estimate']
    se = standard_error(a['fluct']-b['fluct'])
    js = jackknife_error(a['loo']-b['loo']) if a['loo'] is not None and b['loo'] is not None else None
    return dict(status='defined',difference=difference,paired_naive_se=se,
        paired_jackknife_se=js,difference_over_se=difference/se if se else None,
        low_signal=bool(a['row']['low_signal'] or b['row']['low_signal']))


def contrasts(data,saved):
    winrows,cutrows = [],[]
    plan = json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    for n,r in sorted(set(k[:2] for k in data['energies'])):
        for early,late in plan['window_contrasts']:
            row = dict(Nv=n,R=r,reference_tau_min=early[0],reference_tau_max=early[1],
                comparison_tau_min=late[0],comparison_tau_max=late[1])
            row.update(contrast(saved.get((n,r)+tuple(late)),saved.get((n,r)+tuple(early))))
            winrows.append(row)
    for n in CUTOFFS[:-1]:
        for r in LONG_RADII:
            for lo,hi in windows(31):
                row = dict(Nv=n,reference_Nv=100,R=r,tau_min=lo,tau_max=hi)
                row.update(contrast(saved.get((n,r,lo,hi)),saved.get((100,r,lo,hi))))
                cutrows.append(row)
    return winrows,cutrows


def axis_checks(data):
    """Axes remain paired within each configuration, never extra independent samples."""
    rows = []
    for n in CUTOFFS:
        for r in FOCUS:
            for t in range(32):
                vals,deltas = {},{}
                for axis in (0,1,2):
                    x,y = (np.array(data['raw'][axis,n,r,u],float) for u in (t,t+1))
                    mx,my = x.mean(),y.mean()
                    if mx > 0 and my > 0:
                        vals[axis] = float(np.log(mx/my))
                        deltas[axis] = (x-mx)/mx-(y-my)/my
                for a,b in itertools.combinations((0,1,2),2):
                    ok = a in vals and b in vals
                    se = standard_error(deltas[a]-deltas[b]) if ok else None
                    diff = vals[a]-vals[b] if ok else None
                    rows.append(dict(Nv=n,R=r,tau=t,axis=a,reference_axis=b,
                        axis_energy=vals.get(a),reference_axis_energy=vals.get(b),
                        difference=diff,paired_naive_se=se,difference_over_se=diff/se if se else None,
                        status='defined' if ok else 'undefined_endpoint_mean'))
    return rows


def compute(data):
    fits,saved,audit = fit_catalogue(data)
    winrows,cutrows = contrasts(data,saved)
    means = [data['means'][k] for k in sorted(data['means'])]
    energies = [data['energies'][k] for k in sorted(data['energies'])]
    coverage,spatial = [],[]
    plan = json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    for n,r in sorted(set(k[:2] for k in data['energies'])):
        keys = [k for k in data['energies'] if k[:2]==(n,r)]
        coverage.append(dict(Nv=n,R=r,configurations=71,max_correlator_tau=max(k[2] for k in keys)+1,
            max_effective_lower_tau=max(k[2] for k in keys),effective_points=len(keys),
            undefined_points=sum(not data['energies'][k]['valid_point'] for k in keys),
            low_signal_defined_points=sum(data['energies'][k]['valid_point'] and data['energies'][k]['low_signal'] for k in keys)))
        for lo,hi in plan['spatial_display_windows']:
            record = dict(Nv=n,R=r,tau_min=lo,tau_max=hi)
            fit = saved.get((n,r,lo,hi))
            if fit:
                record.update(fit['row'])
            else:
                exists = (n,r,hi) in data['energies']
                record.update(status='undefined_endpoint_mean' if exists else 'outside_measured_time',
                    estimate=None,naive_se=None,jackknife_se=None,chi2=None,dof=hi-lo,
                    p_gaussian_approx=None,condition_inf=None,minimum_endpoint_mean_over_se=None,
                    low_signal=True,diagnostic_screen=False,points=hi-lo+1,configurations=71)
            spatial.append(record)
    tables = dict(correlator_means=means,effective_energies=energies,data_coverage=coverage,
        constant_fit_catalogue=fits,paired_window_contrasts=winrows,
        paired_cutoff_contrasts=cutrows,spatial_window_estimates=spatial,
        axes_near_R24=axis_checks(data))
    summary = dict(configuration_labels=list(range(100,801,10)),lattice=[192,64,64,64],
        piter=10,source_time=90,statistical_unit='configuration after source and axis averaging',
        half_spatial_extent=32,half_temporal_extent=96,maximum_measured_correlator_tau=32,
        maximum_effective_energy_lower_tau=31,raw_audit=data['audit'],fit_regression=audit,
        table_rows={k:len(v) for k,v in tables.items()},
        constant_fit_status_counts={s:sum(r['status']==s for r in fits) for s in sorted(set(r['status'] for r in fits))},
        diagnostic_screen_count=sum(r['diagnostic_screen'] for r in fits),
        diagnostic_scope='Exploratory correlated window scans. No selected plateau, potential or string tension.',
        signal_at_R24_tau24=[data['means'][n,24,24] for n in CUTOFFS])
    return tables,summary
