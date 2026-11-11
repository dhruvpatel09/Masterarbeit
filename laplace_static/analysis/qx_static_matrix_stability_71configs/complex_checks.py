"""Paired complex components and a joint recentered max-|t| bootstrap."""
import numpy as np
from scipy.stats import beta
from common import RADII,CUTOFFS,number


def antihermitian_components(matrices):
    A=(matrices-matrices.conj().swapaxes(-1,-2))/2
    values=[];names=[]
    for i in range(4):
        names.append(('diag_im',i,i));values.append(A[:,i,i].imag)
    for i in range(4):
        for j in range(i+1,4):
            names.extend([('offdiag_re',i,j),('offdiag_im',i,j)])
            values.extend([A[:,i,j].real,A[:,i,j].imag])
    return names,np.array(values).T


def max_t_bootstrap(X,indices):
    """Resample whole centered vectors; never invert their covariance."""
    X=np.asarray(X,float);n=len(X)
    mean=X.mean(0);se=X.std(0,ddof=1)/np.sqrt(n)
    if not np.isfinite(X).all() or np.any(se<=0):raise ValueError('Nonfinite/zero-variance tested component')
    centered=X-mean
    counts=np.array([np.bincount(i,minlength=n) for i in indices],float)
    if counts.shape[1]!=n or not np.all(counts.sum(1)==n):raise ValueError('Bootstrap draw inventory')
    bm=counts.dot(centered)/n
    variance=(counts.dot(centered**2)-n*bm**2)/(n-1)
    if np.any(variance<=0) or not np.isfinite(variance).all():
        raise ValueError('Degenerate bootstrap studentizer; draws are not discarded')
    statistics=np.abs(bm/np.sqrt(variance/n));maxima=statistics.max(1)
    observed=np.abs(mean/se)
    adjusted=(1+np.sum(maxima[:,None]>=observed[None,:],axis=0))/(len(indices)+1.)
    global_tail=float((1+np.sum(maxima>=observed.max()))/(len(indices)+1.))
    critical=float(np.quantile(maxima,.95))
    return mean,se,observed,adjusted,maxima,global_tail,critical


def analyze_complex(C,plan,draws):
    entries=[];tested=[];tests=[];imaginary=[]
    for ri,r in enumerate(RADII):
        scale=np.sqrt(np.diag(C[:,ri,0].mean(0)).real)
        for t in range(33):
            names,X=antihermitian_components(C[:,ri,t])
            for j,(kind,i,k) in enumerate(names):
                x=X[:,j]/(scale[i]*scale[k]);se=float(x.std(ddof=1)/np.sqrt(len(x)));mean=float(x.mean())
                structural=(r==32 and kind!='offdiag_re')
                if structural and np.max(abs(x))>64*np.finfo(float).eps:
                    raise ValueError('Half-box imaginary component exceeds declared roundoff bound')
                test=t in plan['component_test_times'] and not structural
                row=dict(R=r,tau=t,component=kind,Nv_row=CUTOFFS[i],Nv_column=CUTOFFS[k],
                    scaled_mean=mean,independent_se=se,pointwise_t=number(mean/se) if se and not structural else '',
                    max_t_adjusted_tail='',family_tested=bool(test),
                    status='roundoff_only_half_box' if structural else ('family_test' if test else 'descriptive_time'))
                entries.append(row)
                if test:
                    if se<=0:raise ValueError('Zero variance in the prescribed test family')
                    tested.append(x);tests.append(len(entries)-1)
            for i in range(4):
                for k in range(4):
                    x=C[:,ri,t,i,k].imag/(scale[i]*scale[k])
                    se=float(x.std(ddof=1)/np.sqrt(len(x)))
                    imaginary.append(dict(R=r,tau=t,Nv_row=CUTOFFS[i],Nv_column=CUTOFFS[k],
                        scaled_mean=float(x.mean()),independent_se=se,
                        pointwise_t=number(x.mean()/se) if se and r!=32 else '',
                        interpretation='roundoff-only imaginary part at half box; no significance test' if r==32 else
                            'descriptive; an imaginary off-diagonal entry can belong to a Hermitian matrix'))
    X=np.array(tested).T
    if X.shape!=(71,324):raise ValueError('Prescribed anti-Hermitian family changed')
    mean,se,observed,adjusted,maxima,tail,critical=max_t_bootstrap(X,draws)
    for index,value in zip(tests,adjusted):entries[index]['max_t_adjusted_tail']=float(value)
    selected=[]
    for r in RADII:
        for t in plan['component_test_times']:
            rs=[entries[j] for j in tests if entries[j]['R']==r and entries[j]['tau']==t]
            selected.append(dict(R=r,tau=t,components=len(rs),maximum_absolute_pointwise_t=max(abs(d['pointwise_t']) for d in rs),
                minimum_global_adjusted_tail=min(d['max_t_adjusted_tail'] for d in rs),global_95percent_max_t=critical))
    summary=dict(tested_components=324,bootstrap_draws=len(draws),maximum_observed_absolute_t=float(observed.max()),
        joint_null_tail=tail,bootstrap_95percent_max_t=critical,
        family_adjusted_components_below_005=int(np.sum(adjusted<.05)),
        interpretation='Approximate recentered empirical global-null calibration across the fixed 324-component family; non-rejection does not prove Hermiticity or boundary independence.')
    exceedances=int(np.sum(maxima>=observed.max()));B=len(draws)
    summary.update(global_exceedances=exceedances,
        finite_draw_tail_interval_95=[float(beta.ppf(.025,exceedances,B-exceedances+1)) if exceedances else 0.,
            float(beta.ppf(.975,exceedances+1,B-exceedances)) if exceedances<B else 1.],
        interval_scope='Binomial finite-draw uncertainty conditional on the empirical centered distribution; excludes finite-ensemble calibration bias.')
    return dict(antihermitian_components=entries,matrix_imaginary_components=imaginary,complex_family_summary=selected),maxima,summary
