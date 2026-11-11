"""Full71 configuration statistics, conditioning and fixed-source diagnostics."""
import csv
from decimal import Decimal, localcontext
import numpy as np
from scipy import linalg
from scipy.stats import chi2
from replay import HERE, ROOT, CONFIGS, RADII, CUTOFFS
from statistics_tools import (energies, finite_value as f, independent_se, jk_cov, jk_se,
                              metric_transform, assert_close)
from projections import BASES, all_projections, submatrix

WINDOWS = [(6, 11), (10, 15), (14, 19), (18, 23), (22, 27)]
PHASE_LIMIT = .05


def baseline(C, text):
    rows, energy_rows, means_array, reps_array = [], [], [], []
    with localcontext() as context:
        context.prec = 50
        n = Decimal(71)
        for ri, r in enumerate(RADII):
            for j, nv in enumerate(CUTOFFS):
                samples = [[Decimal(text[i, ri, t, j]) for i in range(71)] for t in range(33)]
                totals = [sum(s,Decimal(0)) for s in samples]
                means = [total/n for total in totals]
                ses = [(sum((v-m)**2 for v in s)/(n*(n-1))).sqrt() for s,m in zip(samples,means)]
                loo = [[(total-v)/(n-1) for v in s] for s,total in zip(samples,totals)]
                means_array.append([float(m) for m in means])
                reps_array.append(np.array([[float(v) for v in a] for a in loo]).T)
                for t in range(33):
                    rows.append(dict(R=r, Nv=nv, tau=t, configurations=71, mean_re=str(means[t]),
                        independent_se=str(ses[t]), mean_im=float(C[:,ri,t,j,j].imag.mean()),
                        nonpositive_configurations=sum(v <= 0 for v in samples[t])))
                for t in range(32):
                    valid = means[t] > 0 and means[t+1] > 0
                    influence = [(a-means[t])/means[t]-(b-means[t+1])/means[t+1] for a,b in zip(samples[t],samples[t+1])] if valid else []
                    avg_delta = sum(influence,Decimal(0))/n if valid else Decimal(0)
                    naive = (sum(((v-avg_delta)**2 for v in influence),Decimal(0))/(n*(n-1))).sqrt() if valid else None
                    reps = [((totals[t]-a)/(totals[t+1]-b)).ln() for a,b in zip(samples[t],samples[t+1]) if totals[t]-a>0 and totals[t+1]-b>0]
                    jk = (((n-1)/n)*sum((v-sum(reps)/n)**2 for v in reps)).sqrt() if len(reps)==71 else None
                    energy_rows.append(dict(R=r, Nv=nv, tau=t, midpoint=t+.5, valid_point=valid,
                        aE_eff=str((means[t]/means[t+1]).ln()) if valid else '',
                        independent_se=str(naive) if naive is not None else '',
                        jackknife_se=str(jk) if jk is not None else '', valid_jackknife=len(reps)==71,
                        positive_jackknife_replicates=len(reps),
                        minimum_endpoint_signal=float(min(means[t]/ses[t], means[t+1]/ses[t+1]))))
    # Exact decimal regression against previously published diagonal results.
    old = ROOT/'laplace_static/analysis/qx_static_extension_71configs/results'
    for filename, values, mapping in (
            ('correlator_means_71configs.csv', rows, {'mean_re':'mean_L','independent_se':'naive_se'}),
            ('effective_energies_71configs.csv', energy_rows, {'aE_eff':'aE_eff','independent_se':'naive_se','jackknife_se':'jackknife_se'})):
        with (old/filename).open(newline='') as stream:
            reference = {(int(a['R']),int(a['Nv']),int(a['tau'])):a for a in csv.DictReader(stream)}
        for row in values:
            ref = reference[row['R'],row['Nv'],row['tau']]
            for newkey, oldkey in mapping.items():
                a,b = row[newkey],ref[oldkey]
                if bool(a) != bool(b) or (a and Decimal(a) != Decimal(b)):
                    raise ValueError('Exact diagonal regression: {} R{} Nv{} t{} {}'.format(filename,row['R'],row['Nv'],row['tau'],newkey))
    return rows, energy_rows, np.asarray(means_array), np.asarray(reps_array)


def matrix_diagnostics(C):
    mean = C.mean(axis=0)
    metrics, defects, elements = [], [], []
    for ri,r in enumerate(RADII):
        for basis, idx in BASES.items():
            raw = submatrix(C[:,ri,0],idx)
            spectra = []
            for label, metric in [('ensemble',raw.mean(axis=0))]+[('delete_n'+str(c),(raw.sum(axis=0)-raw[i])/70) for i,c in enumerate(CONFIGS)]:
                w,e,error=metric_transform(metric)
                spectra.append(e)
                for k,value in enumerate(e):
                    metrics.append(dict(R=r,basis=basis,sample=label,mode=k,correlation_eigenvalue=float(value),
                        correlation_condition=float(e[-1]/e[0]),whitening_residual=float(error)))
        d = np.sqrt(np.diag(mean[ri,0]).real)
        for t in range(33):
            A = mean[ri,t]/d[:,None]/d[None,:]
            H = (A+A.conj().T)/2
            norm = np.linalg.norm(A)
            defects.append(dict(R=r,tau=t,relative_antihermitian_norm=float(np.linalg.norm(A-A.conj().T)/(2*norm)),
                relative_imaginary_norm=float(np.linalg.norm(A.imag)/norm),
                hermitian_part_minimum_eigenvalue=float(np.linalg.eigvalsh(H)[0]),
                metric_scaled_frobenius_norm=float(norm)))
            for sink in range(4):
                for src in range(4):
                    x = C[:,ri,t,sink,src]
                    anti = (x-C[:,ri,t,src,sink].conj())/2
                    elements.append(dict(R=r,tau=t,Nv_sink=CUTOFFS[sink],Nv_src=CUTOFFS[src],
                        mean_re=float(x.real.mean()),mean_im=float(x.imag.mean()),
                        independent_se_re=float(independent_se(x.real)),independent_se_im=float(independent_se(x.imag)),
                        antihermitian_mean_re=float(anti.real.mean()),antihermitian_mean_im=float(anti.imag.mean()),
                        antihermitian_se_re=float(independent_se(anti.real)),antihermitian_se_im=float(independent_se(anti.imag))))
    return metrics, defects, elements


def complex_pencils(C):
    # Central algebraic diagnostics only: no eigenvalue energy ordering or level claim.
    rows=[]
    mean=C.mean(axis=0)
    for ri,r in enumerate(RADII):
        w,_,_=metric_transform(mean[ri,0])
        A=np.einsum('ai,tab,bj->tij',w.conj(),mean[ri],w,optimize=True)
        for ref in (0,2,4):
            for t in range(ref+1,33):
                values,vectors=linalg.eig(A[t],A[ref])
                order=np.lexsort((values.imag,-values.real))
                for k,j in enumerate(order):
                    v=vectors[:,j]; value=values[j]
                    residual=np.linalg.norm(A[t].dot(v)-value*A[ref].dot(v))/((np.linalg.norm(A[t])+abs(value)*np.linalg.norm(A[ref]))*np.linalg.norm(v))
                    if not np.isfinite(value) or residual>1e-10:
                        raise ValueError('General complex pencil residual')
                    rows.append(dict(R=r,reference_tau=ref,tau=t,algebraic_mode=k,
                        eigenvalue_re=float(value.real),eigenvalue_im=float(value.imag),phase=float(np.angle(value)),
                        backward_residual=float(residual),reference_condition=float(np.linalg.cond(A[ref]))))
    return rows


def projection_tables(keys, means, reps):
    mean_rows, energy_rows = [], []
    E, J = energies(means), energies(reps)
    for i,key in enumerate(keys):
        se_re,se_im=jk_se(reps[i].real),jk_se(reps[i].imag)
        for t in range(33):
            mean_rows.append(dict(key,tau=t,mean_re=float(means[i,t].real),mean_im=float(means[i,t].imag),
                jackknife_se_re=float(se_re[t]),jackknife_se_im=float(se_im[t]),phase=float(np.angle(means[i,t])),
                real_mean_over_se=float(means[i,t].real/se_re[t]),
                uncertainty='delete-one configuration with complete source retraining'))
        for t in range(32):
            valid=np.isfinite(E[i,t]); validjk=np.isfinite(J[i,:,t]).all()
            phase=max(abs(np.angle(means[i,t])),abs(np.angle(means[i,t+1])))
            signal=min(means[i,t].real/se_re[t],means[i,t+1].real/se_re[t+1])
            energy_rows.append(dict(key,tau=t,midpoint=t+.5,aE_eff=f(E[i,t]),
                jackknife_se=f(jk_se(J[i,:,t])) if validjk else '',valid_point=bool(valid),
                valid_jackknife=bool(validjk),positive_jackknife_replicates=int(np.isfinite(J[i,:,t]).sum()),
                maximum_endpoint_phase=float(phase),phase_within_005=bool(phase<=PHASE_LIMIT),
                minimum_endpoint_signal=float(signal),signal_at_least_3=bool(signal>=3)))
    return mean_rows,energy_rows,E,J


def contrasts(keys,E,J,base_mean,base_reps):
    rows=[]; lookup={tuple(k.values()):i for i,k in enumerate(keys)}
    B,BJ=energies(base_mean),energies(base_reps)
    for i,key in enumerate(keys):
        ri=RADII.index(key['R'])
        comparisons=[('minus_Nv40',B[ri*4],BJ[ri*4]),('minus_Nv100',B[ri*4+3],BJ[ri*4+3])]
        if key['procedure']=='cross_fit':
            j=lookup[(key['R'],key['basis'],key['anchor'],'same_sample')]
            comparisons.append(('minus_same_sample',E[j],J[j]))
        if key['anchor']!=4:
            j=lookup[(key['R'],key['basis'],4,key['procedure'])]
            comparisons.append(('minus_anchor4',E[j],J[j]))
        if key['basis']!='pair':
            j=lookup[(key['R'],'pair',key['anchor'],key['procedure'])]
            comparisons.append(('minus_pair_basis',E[j],J[j]))
        for label,other,other_j in comparisons:
            diff=E[i]-other; reps=J[i]-other_j
            for t in range(32):
                valid=np.isfinite(reps[:,t]).all()
                se=float(jk_se(reps[:,t])) if valid else np.nan
                rows.append(dict(key,contrast=label,tau=t,midpoint=t+.5,difference=f(diff[t]),
                    paired_jackknife_se=f(se),valid_jackknife=bool(valid),
                    pointwise_difference_over_se=f(diff[t]/se) if se>0 else ''))
    return rows


def windows(keys, E, J, energy_rows):
    rows=[]; fits={}
    for i,key in enumerate(keys):
        for start,end in WINDOWS:
            ids=np.arange(start,end+1); y=E[i,ids]; samples=J[i][:,ids]
            row=dict(key,start_tau=start,end_tau=end,points=len(ids),aE='',jackknife_se='',
                chi2_dof='',gaussian_p_approx='',covariance_condition='',minimum_endpoint_signal='',
                maximum_endpoint_phase='',phase_and_signal_screen=False,status='invalid_positive_mean_or_resample')
            if not np.isfinite(y).all() or not np.isfinite(samples).all():
                rows.append(row);continue
            cov=jk_cov(samples)
            ev=np.linalg.eigvalsh(cov)
            if ev[0]<=1e-12*ev[-1]:
                row['status']='singular_covariance';rows.append(row);continue
            ones=np.ones(len(ids)); solved=np.linalg.solve(cov,ones); weights=solved/solved.sum()
            center=float(weights.dot(y)); jkmeans=samples.dot(weights)
            residual=y-center
            statistic=float(residual.dot(np.linalg.solve(cov,residual)))
            detail=energy_rows[i*32+start:i*32+end+1]
            signal=min(d['minimum_endpoint_signal'] for d in detail); phase=max(d['maximum_endpoint_phase'] for d in detail)
            row.update(aE=center,jackknife_se=float(jk_se(jkmeans)),chi2_dof=statistic/(len(ids)-1),
                gaussian_p_approx=float(chi2.sf(statistic,len(ids)-1)),covariance_condition=float(ev[-1]/ev[0]),
                minimum_endpoint_signal=signal,maximum_endpoint_phase=phase,
                phase_and_signal_screen=bool(signal>=3 and phase<=PHASE_LIMIT),status='defined_conditional_diagnostic')
            rows.append(row)
            fits[(tuple(key.values()),start)]=(center,jkmeans)
    return rows,fits


def fit_contrasts(keys,fits):
    rows=[]
    for key in keys:
        if key['procedure']=='diagonal':continue
        for start,end in WINDOWS:
            own=fits.get((tuple(key.values()),start))
            other_keys=[('minus_Nv40',(key['R'],'Nv40',0,'diagonal'),start),
                        ('minus_Nv100',(key['R'],'Nv100',0,'diagonal'),start)]
            if key['procedure']=='cross_fit':
                other_keys.append(('minus_same_sample',(key['R'],key['basis'],key['anchor'],'same_sample'),start))
            if key['basis']!='pair':
                other_keys.append(('minus_pair_basis',(key['R'],'pair',key['anchor'],key['procedure']),start))
            if key['anchor']!=4:
                other_keys.append(('minus_anchor4',(key['R'],key['basis'],4,key['procedure']),start))
            if start!=6:
                other_keys.append(('minus_window6_11',tuple(key.values()),6))
            for contrast,other_key,other_start in other_keys:
                other=fits.get((other_key,other_start));valid=own is not None and other is not None
                diff=own[0]-other[0] if valid else np.nan
                se=float(jk_se(own[1]-other[1])) if valid else np.nan
                rows.append(dict(key,start_tau=start,end_tau=end,contrast=contrast,
                    difference=f(diff),paired_jackknife_se=f(se),defined=valid,
                    pointwise_difference_over_se=f(diff/se) if se>0 else ''))
    return rows


def analyze(C,text,progress=False):
    b,e,bm,bj=baseline(C,text)
    metrics,defects,elements=matrix_diagnostics(C)
    pencils=complex_pencils(C)
    keys,means,reps,vectors=all_projections(C,RADII,progress)
    m,pe,E,J=projection_tables(keys,means,reps)
    paired=contrasts(keys,E,J,bm,bj)
    fits,fit_samples=windows(keys,E,J,pe)
    base_keys=[];base_quality=[]
    for ri,r in enumerate(RADII):
        for j,nv in enumerate(CUTOFFS):
            base_keys.append(dict(R=r,basis='Nv'+str(nv),anchor=0,procedure='diagonal'))
            for t in range(32):
                detail=e[(ri*4+j)*32+t]
                phase=max(abs(np.angle(C[:,ri,t,j,j].mean())),abs(np.angle(C[:,ri,t+1,j,j].mean())))
                base_quality.append(dict(minimum_endpoint_signal=detail['minimum_endpoint_signal'],maximum_endpoint_phase=float(phase)))
    base_fits,base_fit_samples=windows(base_keys,energies(bm),energies(bj),base_quality)
    fit_samples.update(base_fit_samples)
    paired_fits=fit_contrasts(keys,fit_samples)
    tables=dict(diagonal_correlators=b,diagonal_effective_energies=e,metric_spectra=metrics,
        matrix_defects=defects,matrix_elements=elements,complex_pencil_diagnostics=pencils,
        fixed_source_vectors=vectors,projected_correlators=m,projected_effective_energies=pe,
        paired_energy_contrasts=paired,constant_window_diagnostics=fits,
        diagonal_window_diagnostics=base_fits,paired_window_contrasts=paired_fits)
    arrays=dict(projected_means=means,projected_jackknife_means=reps,diagonal_means=bm,
        diagonal_jackknife_means=bj)
    summary=dict(configurations=71,radii=RADII,cutoffs=CUTOFFS,piter=10,t0=90,
        temporal_extent=192,spatial_extent=64,correlator_taus=list(range(33)),effective_energy_taus=list(range(32)),
        exact_diagonal_correlator_points=len(b),exact_diagonal_effective_energy_points=len(e),
        source_candidates=len(keys)//2,source_procedures=2,outer_jackknife_resamples=71,
        projected_effective_energy_rows=len(pe),undefined_projected_energy_points=sum(not d['valid_point'] for d in pe),
        projected_energy_points_with_invalid_jackknife=sum(not d['valid_jackknife'] for d in pe),
        constant_window_attempts=len(fits),defined_constant_windows=sum(d['status']=='defined_conditional_diagnostic' for d in fits),
        metric_conditions_full_basis={str(r):next(d['correlation_condition'] for d in metrics if d['R']==r and d['basis']=='full' and d['sample']=='ensemble') for r in RADII},
        analysis_interpretation='Source-selection and finite-sample sensitivity diagnostics. No final ground-state static potential, positive spectral decomposition, string tension or Cornell fit.',
        statistical_unit='one configuration; axes, sites, sources, folds and reused pilot cases are not additional independent configurations',
        cross_fit='Alternating original configuration indices; fit on one fold and evaluate the other. Outer delete-one retrains both folds. Cross-fit combines two trained sources, not a single fixed operator.',
        constant_windows='Generalized least-squares constants with full-sample outer-jackknife covariance, held fixed while applying delete-one fits. Gaussian p values are approximate diagnostics; no selection rule or multiplicity-adjusted discovery claim.',
        covariance_limit='71 independent configurations: empirical covariance rank at most 70. No inverse of the full 16 x 33 feature covariance.')
    return tables,arrays,summary
