"""Conditional A+B exp(-c t) effective-energy diagnostics, not a spectral extraction."""
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from statistics_tools import covariance, whiten


def model(t,p):
    a,b,c = p
    return a+b*np.exp(-c*np.asarray(t,float))


def profile(t,y,solve,c):
    e = np.exp(-c*t)
    X = np.column_stack((np.ones(len(t)),e))
    normal = X.T@solve(X)
    p = np.linalg.solve(normal,X.T@solve(y))
    residual = y-X@p
    return float(residual@solve(residual)),np.array([p[0],p[1],c])


def search(t,y,C,bounds,grid_count=65):
    solve = whiten(C)
    grid = np.geomspace(bounds[0],bounds[1],grid_count)
    candidates = []
    for c in grid:
        q,p = profile(t,y,solve,c)
        candidates.append(dict(params=p.tolist(),chi2=q,kind='grid'))
    # Refine every grid local minimum, not merely one initial guess.
    vals = [r['chi2'] for r in candidates]
    for i in range(1,len(grid)-1):
        if vals[i] <= vals[i-1] and vals[i] <= vals[i+1]:
            result = minimize_scalar(lambda c:profile(t,y,solve,c)[0],
                bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-12})
            q,p = profile(t,y,solve,float(result.x))
            candidates.append(dict(params=p.tolist(),chi2=q,kind='refined'))
    return candidates


def assess(t,y,d,minimum,candidates,bounds):
    C = covariance(d)
    solve = whiten(C)
    # Candidate replay is algebraic; does not repeat optimizer decisions across BLAS versions.
    grid = [row['params'][2] for row in candidates if row.get('kind')=='grid']
    if len(grid)!=65 or not np.allclose(grid,np.geomspace(bounds[0],bounds[1],65),rtol=2e-14,atol=0):
        raise ValueError('Incomplete archived rate grid')
    if any(row.get('kind') not in ('grid','refined') for row in candidates):
        raise ValueError('Unknown approach candidate kind')
    for candidate in candidates:
        p = np.array(candidate['params'],float)
        if p.shape != (3,) or not np.isfinite(p).all() or not bounds[0] <= p[2] <= bounds[1]:
            raise ValueError('Malformed archived approach candidate')
        q,fresh = profile(t,y,solve,p[2])
        if not np.allclose(p,fresh,rtol=2e-8,atol=2e-10) or abs(q-candidate['chi2'])>2e-7*(1+abs(q)):
            raise ValueError('Approach candidate does not solve conditional linear system')
    chosen = min(candidates,key=lambda r:r['chi2'])
    p = np.array(chosen['params'],float)
    a,b,c = p
    residual = y-model(t,p)
    e = np.exp(-c*t)
    J = np.column_stack((np.ones(len(t)),e,-b*t*e))
    H = J.T@solve(J)
    wr = solve(residual)
    H[1,2] -= wr@(-t*e); H[2,1] = H[1,2]
    H[2,2] -= wr@(b*t*t*e)
    flags = []
    if minimum < 3:flags.append('weak_endpoint_signal')
    if c <= bounds[0]*1.001 or c >= bounds[1]/1.001:flags.append('c_at_search_boundary')
    gradient = float(np.max(np.abs(J.T@wr)/np.maximum(1.,np.sqrt(np.diag(J.T@solve(J))))))
    if gradient > 1e-4:flags.append('stationarity_unresolved')
    se = [None,None,None]; influence = None; condition = None
    try:
        np.linalg.cholesky(H)
        condition = float(np.linalg.cond(H))
        response = np.linalg.solve(H,J.T@solve(np.eye(len(t))))
        pcov = response@C@response.T
        se = np.sqrt(np.maximum(0,np.diag(pcov))).tolist()
        influence = d@response.T
        if condition > 1e12:flags.append('ill_conditioned_curvature')
        if abs(b) <= 2*se[1] or c <= 2*se[2]:flags.append('amplitude_or_rate_unresolved')
    except np.linalg.LinAlgError:
        flags.append('nonpositive_curvature')
    q = float(residual@wr); dof = len(y)-3; prob = float(chi2.sf(q,dof))
    if prob < .05:flags.append('poor_approximate_fit')
    row = dict(status='defined',A=a,B=b,c=c,A_conditional_se=se[0],B_conditional_se=se[1],
        c_conditional_se=se[2],chi2=q,dof=dof,p_gaussian_approx=prob,
        minimum_endpoint_mean_over_se=minimum,gradient_scaled=gradient,
        curvature_condition=condition,resolution_flags=';'.join(flags),
        diagnostic_screen=not flags,candidate_count=len(candidates))
    return row,influence


def compute(data,plan,archived=None):
    rows,solutions = [],{}
    for n,r in sorted(set(k[:2] for k in data['energies'])):
        last = max(k[2] for k in data['energies'] if k[:2]==(n,r))
        starts = plan['approach_starts_long'] if last==31 else plan['approach_starts_short']
        for start in starts:
            lo,hi = start,last
            key = '{}:{}:{}:{}'.format(n,r,lo,hi)
            keys = [(n,r,t) for t in range(lo,hi+1)]
            minimum = min(data['means'][n,r,t]['mean_over_se'] for t in range(lo,hi+2))
            row = dict(Nv=n,R=r,tau_min=lo,tau_max=hi,status='undefined_endpoint_mean',A=None,B=None,c=None,
                A_conditional_se=None,B_conditional_se=None,c_conditional_se=None,chi2=None,dof=len(keys)-3,
                p_gaussian_approx=None,minimum_endpoint_mean_over_se=minimum,gradient_scaled=None,
                curvature_condition=None,resolution_flags='undefined_endpoint_mean',diagnostic_screen=False,
                candidate_count=0)
            if all(k in data['delta'] for k in keys):
                t = np.arange(len(keys),dtype=float)
                y = np.array([float(data['energies'][k]['aE_eff']) for k in keys])
                d = np.column_stack([data['delta'][k] for k in keys])
                try:
                    candidates = search(t,y,covariance(d),plan['approach_c_bounds']) if archived is None else archived[key]
                    result,_ = assess(t,y,d,minimum,candidates,plan['approach_c_bounds'])
                except np.linalg.LinAlgError:
                    row.update(status='invalid_covariance',resolution_flags='invalid_covariance')
                else:
                    row.update(result);solutions[key]=candidates
            rows.append(row)
    if archived is not None and set(archived) != set(solutions):
        raise ValueError('Approach solution inventory mismatch')
    return rows,solutions
