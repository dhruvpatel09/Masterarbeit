"""Load pinned published qX measurements and audit the overlapping scans."""
import csv
import hashlib
import json
from decimal import Decimal as D, localcontext
from pathlib import Path

import numpy as np
from statistics_tools import (assert_close, covariance, energy, mean_se,
                              jackknife_error)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONFIGS = list(range(100, 801, 10))
LONG_RADII = [1, 2, 4, 8, 12, 16, 20, 23, 24, 25, 28, 32]
CUTOFFS = [40, 60, 80, 100]
FOCUS = [20, 23, 24, 25, 28]
SHORT = 'laplace_static/analysis/qx_static_potential_71configs/'
LONG = 'laplace_static/analysis/qx_static_extension_71configs/'


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def check_inputs():
    spec = json.loads((HERE/'inputs/SOURCES.json').read_text())
    if spec['configuration_labels'] != CONFIGS:
        raise ValueError('Configuration labels changed')
    for name, receipt in spec['files'].items():
        p = ROOT/name
        if p.stat().st_size != receipt['bytes'] or digest(p) != receipt['sha256']:
            raise ValueError('Published dependency changed: '+name)
    return spec


def load_raw():
    spec = check_inputs()
    groups = {}
    overlap = 0
    for label in ('short', 'long'):
        seen = set()
        with (ROOT/spec['raw_tables'][label]).open(newline='') as stream:
            for row in csv.DictReader(stream):
                if row['scheme'] != 'fixed' or row['piter'] != '10':
                    continue
                axis, nv, r, t, cfg = (int(row[k]) for k in ('axis','Nv','R','tau','cfg'))
                if cfg not in CONFIGS or row['t0'] != '90':
                    raise ValueError('Unexpected cohort or source time')
                if (label == 'short' and (nv != 100 or r not in range(1,33) or t not in range(17))) or (label == 'long' and (nv not in CUTOFFS or r not in LONG_RADII or t not in range(33))):
                    raise ValueError('Unexpected measurement coverage')
                identity = (axis,nv,r,t,cfg)
                if identity in seen:
                    raise ValueError('Duplicate measurement in '+label)
                seen.add(identity)
                if axis != -1 and r not in FOCUS:
                    continue
                key = (axis,nv,r,t)
                bucket = groups.setdefault(key, {})
                if cfg in bucket:
                    if D(bucket[cfg]) != D(row['L_RE']):
                        raise ValueError('Overlapping raw correlators differ')
                    overlap += 1
                bucket[cfg] = row['L_RE']
        expected = 71*4*(32*17 if label == 'short' else 4*12*33)
        if len(seen) != expected:
            raise ValueError('Incomplete {} measurements: {}'.format(label,len(seen)))
    if any(sorted(g) != CONFIGS for g in groups.values()):
        raise ValueError('Unmatched configuration samples')
    if overlap != 71*(204+3*5*17):
        raise ValueError('Wrong overlap inventory')
    return {k:[g[c] for c in CONFIGS] for k,g in groups.items()}, overlap


def load():
    raw, overlap = load_raw()
    means, energies, delta, loo = {}, {}, {}, {}
    # Preserve published primary decimal strings; independently reconstruct their statistics.
    for row in read_csv(ROOT/SHORT/'results/correlators_ensemble_71configs.csv'):
        if row['axis'] == '-1':
            n,r,t = 100,int(row['R']),int(row['tau'])
            means[n,r,t] = dict(Nv=n,R=r,tau=t,configurations=71,
                mean_L=row['L_RE'],naive_se=row['L_RE_naive_se'],mean_over_se=float(row['mean_over_se']),
                nonpositive_individual=int(row['nonpositive_individual']),normalized_L=None,
                normalized_naive_se=None,source='short')
    for row in read_csv(ROOT/LONG/'results/correlator_means_71configs.csv'):
        n,r,t = (int(row[k]) for k in ('Nv','R','tau'))
        if (n,r,t) in means:
            old = means[n,r,t]
            if D(old['mean_L']) != D(row['mean_L']) or D(old['naive_se']) != D(row['naive_se']):
                raise ValueError('Published means disagree on overlap')
        means[n,r,t] = dict(Nv=n,R=r,tau=t,configurations=71,
            mean_L=row['mean_L'],naive_se=row['naive_se'],mean_over_se=float(row['mean_over_se']),
            nonpositive_individual=int(row['nonpositive_individual_configurations']),
            normalized_L=row['normalized_L'],normalized_naive_se=row['normalized_naive_se'],source='long')
    for label,prefix in (('short',SHORT),('long',LONG)):
        for row in read_csv(ROOT/prefix/'results/effective_energies_71configs.csv'):
            if label == 'short' and row['axis'] != '-1':
                continue
            n,r,t = (100 if label=='short' else int(row['Nv'])),int(row['R']),int(row['tau'])
            key = (n,r,t)
            if key in energies:
                for field in ('aE_eff','naive_se','jackknife_se'):
                    if D(energies[key][field]) != D(row[field]):
                        raise ValueError('Published energies disagree on overlap')
            energies[key] = dict(Nv=n,R=r,tau=t,tau_midpoint=t+.5,configurations=71,
                valid_point=row['valid_point']=='True',aE_eff=row['aE_eff'],naive_se=row['naive_se'],
                valid_jackknife=row['valid_jackknife']=='True',jackknife_se=row['jackknife_se'],source=label)
    if len(means) != 1924 or len(energies) != 1856:
        raise ValueError('Unexpected merged geometry')
    with localcontext() as ctx:
        ctx.prec = 50
        for (n,r,t), row in sorted(means.items()):
            x = [D(v) for v in raw[-1,n,r,t]]
            m,se = mean_se(x)
            if str(m) != row['mean_L'] or str(se) != row['naive_se']:
                # Numerically equivalent trailing zeros are immaterial.
                if m != D(row['mean_L']) or se != D(row['naive_se']):
                    raise ValueError('Raw mean/error regression')
            z = [D(v) for v in raw[-1,n,r,0]]
            mz,_ = mean_se(z)
            ratio = m/mz
            norm_delta = [(a-ratio*b)/mz for a,b in zip(x,z)]
            norm_se = mean_se(norm_delta)[1]
            if row['source'] == 'long':
                assert_close(float(ratio),float(row['normalized_L']),'normalized mean',rtol=1e-13,atol=1e-17)
                assert_close(float(norm_se),float(row['normalized_naive_se']),'normalized error',rtol=1e-13,atol=1e-17)
            else:
                row['normalized_L'],row['normalized_naive_se'] = str(ratio),str(norm_se)
            row['positive_mean'] = m > 0
            row['low_signal'] = row['mean_over_se'] < 3
        for (n,r,t), row in sorted(energies.items()):
            x,y = ([D(v) for v in raw[-1,n,r,u]] for u in (t,t+1))
            values,d,j,positive = energy(x,y)
            row['minimum_endpoint_mean_over_se'] = min(means[n,r,u]['mean_over_se'] for u in (t,t+1))
            row['low_signal'] = row['minimum_endpoint_mean_over_se'] < 3
            row['positive_jackknife_replicates'] = positive
            if row['valid_point'] != (values is not None):
                raise ValueError('Defined-energy inventory changed')
            if values is not None:
                if D(values[0]) != D(row['aE_eff']) or D(values[1]) != D(row['naive_se']):
                    raise ValueError('Raw energy/error regression')
                delta[n,r,t] = d
                if j is not None:
                    loo[n,r,t] = j
                    assert_close(jackknife_error(j),float(row['jackknife_se']),'jackknife error',rtol=2e-9,atol=2e-14)
    cov_checks = 0
    # Audit every published combined-axis covariance entry, including cutoff correlations.
    for prefix,name in ((SHORT,'temporal_covariance_71configs.csv'),(LONG,'effective_energy_covariance_71configs.csv')):
        with (ROOT/prefix/'results'/name).open(newline='') as stream:
            cache = {}
            for row in csv.DictReader(stream):
                if prefix == SHORT:
                    if row['axis'] != '-1':
                        continue
                    r = int(row['R']); a=(100,r,int(row['tau_i'])); b=(100,r,int(row['tau_j']))
                    value = row['covariance']
                else:
                    r = int(row['R']); a=(int(row['Nv']),r,int(row['tau'])); b=(int(row['reference_Nv']),r,int(row['reference_tau']))
                    value = row['covariance_naive']
                if (a in delta and b in delta) != bool(value):
                    raise ValueError('Covariance definedness changed')
                if value:
                    if r not in cache:
                        keys = sorted(k for k in delta if k[1]==r)
                        matrix = covariance(np.column_stack([delta[k] for k in keys]))
                        cache[r] = (dict((k,i) for i,k in enumerate(keys)),matrix)
                    ix,C = cache[r]
                    assert_close(C[ix[a],ix[b]],float(value),'published covariance',rtol=3e-12,atol=1e-19)
                cov_checks += 1
    return dict(raw=raw,means=means,energies=energies,delta=delta,loo=loo,
        audit=dict(configuration_count=71,raw_overlap_records=overlap,mean_overlap_points=204,
                   energy_overlap_points=192,covariance_entries_verified=cov_checks,
                   merged_mean_points=len(means),merged_energy_points=len(energies),
                   defined_energy_points=len(delta)))
