"""Joint complex-component checks, held-out stability and direct-correlator fits."""
import numpy as np
from scipy.stats import chi2
from common import HERE,RADII,read,load_data,number,jk_cov,jk_se,energies,check_equal_time_reality
from complex_checks import analyze_complex
from source_splits import compute_splits,make_assignments
from correlator_fits import ExponentialFit,whitening,SingularCovariance


def curve_tables(keys,means,deletions):
    corrs=[];ee=[];E=energies(means);J=energies(deletions)
    for ci,key in enumerate(keys):
        check_equal_time_reality(means[ci,0],deletions[ci,:,0],str(key))
        re,im=jk_se(deletions[ci].real),jk_se(deletions[ci].imag)
        circular=np.angle(deletions[ci]/means[ci][None,:]);phase_se=jk_se(circular)
        for t in range(33):
            normalization=key['variant']=='same_sample' and t==0
            corrs.append(dict(key,tau=t,mean_re=float(means[ci,t].real),mean_im=float(means[ci,t].imag),
                jackknife_se_re=float(re[t]),jackknife_se_im=float(im[t]),
                real_mean_over_se=number(means[ci,t].real/re[t]) if re[t]>0 and not normalization else '',
                phase=float(np.angle(means[ci,t])),phase_jackknife_se=float(phase_se[t]),
                imaginary_pointwise_t=number(means[ci,t].imag/im[t]) if im[t]>0 and key['R']!=32 and t!=0 else '',
                diagnostic='training_normalization' if normalization else ('roundoff_imaginary_half_box' if key['R']==32 else
                    ('roundoff_imaginary_equal_time' if t==0 else 'physical_time_diagnostic'))))
        for t in range(32):
            valid=np.isfinite(J[ci,:,t]).all()
            ee.append(dict(key,tau=t,midpoint=t+.5,aE_eff=number(E[ci,t]),
                jackknife_se=number(jk_se(J[ci,:,t])) if valid else '',
                valid_point=bool(np.isfinite(E[ci,t])),valid_jackknife=bool(valid)))
    return corrs,ee,E,J


def constant_fit(y,samples):
    if not np.isfinite(y).all() or not np.isfinite(samples).all():return None
    try:W,condition=whitening(jk_cov(samples))
    except SingularCovariance:return None
    one=W.dot(np.ones(len(y)));weights=one.dot(W)/np.dot(one,one)
    value=float(weights.dot(y));rep=samples.dot(weights)
    return value,rep,condition


def fit_curves(keys,means,deletions,E,J,plan,progress=False):
    rows=[];replicates=[];amplitudes=[]
    for ci,key in enumerate(keys):
        for start,end in plan['correlator_fit_windows']:
            times=np.arange(start,end+1);y=means[ci,times].real;samples=deletions[ci][:,times].real
            row=dict(key,start_tau=start,end_tau=end,points=len(times),status='invalid_covariance',
                aE='',energy_jackknife_se='',amplitude='',amplitude_jackknife_se='',
                chi2_dof='',gaussian_p_approx='',covariance_condition='',energy_at_boundary='',
                jackknife_interior_fits='',positive_profile_curvature='',stationary_energy_error='',
                minimum_endpoint_signal='',maximum_endpoint_phase='',constant_aE='',constant_jackknife_se='',
                direct_minus_constant='',paired_direct_minus_constant_se='')
            try:fitter=ExponentialFit(jk_cov(samples),times,plan['energy_search_interval'],plan['energy_grid_points'])
            except SingularCovariance:
                rows.append(row);replicates.append(np.full(71,np.nan));amplitudes.append(np.full(71,np.nan));continue
            fits=fitter.fit_many(np.vstack([y,samples]));fit=fits[0]
            rE=np.array([d['energy'] for d in fits[1:]]);rA=np.array([d['amplitude'] for d in fits[1:]])
            minimum_signal=float(np.min(y/jk_se(samples)))
            maximum_phase=float(np.max(abs(np.angle(means[ci,times]))))
            const=constant_fit(E[ci,start:end],J[ci][:,start:end])
            row.update(status='defined_conditional_fit',aE=fit['energy'],energy_jackknife_se=float(jk_se(rE)),
                amplitude=fit['amplitude'],amplitude_jackknife_se=float(jk_se(rA)),
                chi2_dof=fit['q']/(len(times)-2),gaussian_p_approx=float(chi2.sf(fit['q'],len(times)-2)),
                covariance_condition=fitter.condition,energy_at_boundary=fit['boundary'],
                jackknife_interior_fits=sum(not d['boundary'] for d in fits[1:]),
                positive_profile_curvature=bool(fit['curvature']>0),
                stationary_energy_error=number(fit['energy_residual']) if not fit['boundary'] else '',
                minimum_endpoint_signal=minimum_signal,maximum_endpoint_phase=maximum_phase)
            if const is not None:
                value,reps,condition=const
                row.update(constant_aE=value,constant_jackknife_se=float(jk_se(reps)),
                    direct_minus_constant=fit['energy']-value,paired_direct_minus_constant_se=float(jk_se(rE-reps)))
            rows.append(row);replicates.append(rE);amplitudes.append(rA)
        if progress and (ci+1)%54==0:print('Direct-correlator and paired deletion fits: {}/216 curves'.format(ci+1),flush=True)
    return rows,np.asarray(replicates),np.asarray(amplitudes)


def fit_contrasts(rows,replicates):
    fields=('R','basis','anchor','variant','start_tau');lookup={tuple(d[k] for k in fields):i for i,d in enumerate(rows)}
    output=[]
    for i,row in enumerate(rows):
        key=tuple(row[k] for k in fields);r,b,a,v,t=key;references=[]
        if v!='alternating':references.append(('minus_alternating',(r,b,a,'alternating',t)))
        if b!='pair':references.append(('minus_pair_basis',(r,'pair',a,v,t)))
        if a!=4:references.append(('minus_anchor4',(r,b,4,v,t)))
        if t!=6:references.append(('minus_window6_12',(r,b,a,v,6)))
        for label,reference in references:
            j=lookup[reference];other=rows[j]
            valid=row['status']=='defined_conditional_fit' and other['status']=='defined_conditional_fit'
            difference=float(row['aE']-other['aE']) if valid else np.nan
            error=float(jk_se(replicates[i]-replicates[j])) if valid else np.nan
            output.append(dict(R=r,basis=b,anchor=a,variant=v,start_tau=t,end_tau=row['end_tau'],contrast=label,
                difference=number(difference),paired_jackknife_se=number(error),defined=bool(valid),
                pointwise_difference_over_se=number(difference/error) if error>0 else ''))
    return output


def split_summary(fits,contrasts):
    result=[]
    for r in RADII:
        for basis in ('pair','triple','full'):
            for anchor in (4,6,8):
                for start in (6,10,14,18,22):
                    group=[d for d in fits if d['R']==r and d['basis']==basis and d['anchor']==anchor
                           and d['start_tau']==start and d['variant']!='same_sample' and d['status']=='defined_conditional_fit']
                    comparisons=[d for d in contrasts if d['R']==r and d['basis']==basis and d['anchor']==anchor
                                 and d['start_tau']==start and d['variant'].startswith('seed') and d['contrast']=='minus_alternating'
                                 and d['pointwise_difference_over_se']!='']
                    values=[d['aE'] for d in group]
                    result.append(dict(R=r,basis=basis,anchor=anchor,start_tau=start,end_tau=start+6,
                        defined_holdout_splits=len(values),minimum_energy=min(values) if values else '',
                        maximum_energy=max(values) if values else '',
                        maximum_absolute_paired_shift_over_se=max(abs(d['pointwise_difference_over_se']) for d in comparisons) if comparisons else '',
                        interpretation='five split choices share the same 71 configurations; spread is not an independent-sample error'))
    return result


def compute(progress=False):
    plan=read(HERE/'ANALYSIS_PLAN.json');C=load_data()
    resampling=np.load(str(HERE/'inputs/resampling_plan.npz'),allow_pickle=False)
    if list(resampling['configurations'])!=plan['configuration_labels']:raise ValueError('Resampling cohort')
    expected=np.random.RandomState(plan['bootstrap_seed']).randint(0,71,size=(plan['bootstrap_draws'],71)).astype(np.uint8)
    if (not np.array_equal(resampling['bootstrap_indices'],expected)
            or not np.array_equal(resampling['split_assignments'],make_assignments(plan))):
        raise ValueError('Resampling indices differ from prescribed seeds or folds')
    complex_tables,maxima,complex_summary=analyze_complex(C,plan,resampling['bootstrap_indices'])
    if progress:print('Complex-family bootstrap:',complex_summary,flush=True)
    keys,means,deletions,membership=compute_splits(C,plan,resampling['split_assignments'],progress)
    corrs,ee,E,J=curve_tables(keys,means,deletions)
    fits,fit_reps,fit_amplitudes=fit_curves(keys,means,deletions,E,J,plan,progress)
    contrasts=fit_contrasts(fits,fit_reps);split_ranges=split_summary(fits,contrasts)
    tables=dict(complex_tables,split_membership=membership,projected_correlators=corrs,
        projected_effective_energies=ee,correlator_fits=fits,paired_fit_contrasts=contrasts,holdout_split_summary=split_ranges)
    arrays=dict(projected_means=means,projected_jackknife_means=deletions,fit_energy_deletions=fit_reps,
        fit_amplitude_deletions=fit_amplitudes,complex_bootstrap_max_t=maxima)
    summary=dict(base_commit='2b9830b3e5de49b8234b2257c24a0175f247cc61',configurations=71,
        complex_component_tests=complex_summary,holdout_splits=5,additional_random_splits=4,
        source_candidates=36,total_curves=len(keys),outer_retraining_deletions=71,
        correlator_fit_attempts=len(fits),defined_correlator_fits=sum(d['status']=='defined_conditional_fit' for d in fits),
        boundary_central_fits=sum(d['energy_at_boundary'] is True for d in fits),
        fits_with_incomplete_interior_jackknife=sum(d['jackknife_interior_fits']!='' and d['jackknife_interior_fits']<71 for d in fits),
        tables={name:len(values) for name,values in tables.items()},
        interpretation='Resampling and fit stability diagnostics. No final ground-state potential, proof of Hermiticity, positive-spectrum GEVP, string tension or new MPI measurements.',
        goodness_scope='Correlator-fit Gaussian tails use estimated full-sample fixed covariance; uncalibrated and exploratory. Component max-t bootstrap is a separate null test, not fit goodness calibration.',
        resampling_scope='All variants share the same 71 configurations. Pair contrasts by identical deleted configuration; never treat folds or split choices as extra independent samples.')
    return tables,arrays,summary
