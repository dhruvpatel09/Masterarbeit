"""Complex matrix conditioning and explicitly exploratory generalized pencils."""
import json
import math
from pathlib import Path
import numpy as np
from scipy import linalg
from scipy.optimize import linear_sum_assignment

HERE=Path(__file__).resolve().parent


def assemble(rows,plan):
    configs=plan['configuration_labels'];radii=plan['radii'];cuts=plan['cutoffs']
    data=np.empty((len(configs),len(radii),33,4,4),dtype=complex);control={};seen=set()
    for row in rows:
        if int(row['axis'])!=plan['axis']:continue
        cfg,r,t,s,k=(int(row[name]) for name in ('cfg','R','tau','Nv_src','Nv_sink'))
        z=complex(float(row['L_RE']),float(row['L_IM']))
        if row['piter']=='40':control[r,t,k,s]=z;continue
        key=(cfg,r,t,k,s)
        if key in seen:raise ValueError('Duplicate analysis matrix entry')
        seen.add(key);data[configs.index(cfg),radii.index(r),t,cuts.index(k),cuts.index(s)]=z
    expected={(c,r,t,k,s) for c in configs for r in radii for t in range(33) for k in cuts for s in cuts}
    if seen!=expected or not np.all(np.isfinite(data)):raise ValueError('Incomplete complex analysis matrices')
    if len(control)!=4*33*16:raise ValueError('Incomplete p40 matrix control')
    return data,control


def subset(series,indices):return series[:,indices,:][:,:,indices]


def metric_transform(c0):
    scale=np.linalg.norm(c0)
    if not scale or np.linalg.norm(c0-c0.conj().T)>1e-10*scale:raise ValueError('Equal-time metric must be Hermitian')
    diagonal=c0.diagonal().real
    if np.any(diagonal<=0):raise ValueError('Equal-time diagonal must be positive')
    d=np.sqrt(diagonal);normalized=c0/d[:,None]/d[None,:]
    eigenvalues,u=linalg.eigh((normalized+normalized.conj().T)/2)
    if eigenvalues[0]<=1e-12*eigenvalues[-1]:raise ValueError('Equal-time metric numerically singular; no automatic mode deletion')
    whiten=(u/np.sqrt(eigenvalues))@u.conj().T
    transform=whiten/d[:,None]
    if not np.allclose(transform.conj().T@c0@transform,np.eye(len(d)),rtol=1e-10,atol=1e-10):raise ValueError('Equal-time whitening failed')
    return transform,d,eigenvalues


def overlaps(left,right,metric):
    a=np.sqrt(np.maximum(np.real(np.einsum('ij,ij->j',left.conj(),metric@left)),0))
    b=np.sqrt(np.maximum(np.real(np.einsum('ij,ij->j',right.conj(),metric@right)),0))
    if np.any(a==0) or np.any(b==0):raise ValueError('Zero eigenvector norm')
    return np.minimum(np.abs(left.conj().T@metric@right)/a[:,None]/b[None,:],1.)


def match(left,right,metric):
    score=overlaps(left,right,metric);rows,columns=linear_sum_assignment(-score)
    if not np.array_equal(rows,np.arange(score.shape[0])):raise ValueError('Incomplete mode assignment')
    return columns,score[rows,columns]


def match_spectrum(labels,eigenvalues):
    """Choose archive labels after independently solving the complete pencil."""
    labels=np.asarray(labels,dtype=complex)
    if labels.shape!=eigenvalues.shape or not np.all(np.isfinite(labels)):raise ValueError('Incomplete replay eigenvalue labels')
    distance=np.abs(labels[:,None]-eigenvalues[None,:])
    rows,columns=linear_sum_assignment(distance)
    if not np.array_equal(rows,np.arange(len(labels))):raise ValueError('Incomplete spectrum assignment')
    wanted=eigenvalues[columns]
    for first,second in ((labels.real,wanted.real),(labels.imag,wanted.imag)):
        if np.any(np.abs(first-second)>2e-10+2e-8*np.maximum(np.abs(first),np.abs(second))):raise ValueError('Recomputed spectrum differs from archived labels')
    return columns


def solve_series(series,tref,central=None,central_metric=None,labels=None):
    transform,_,_=metric_transform(series[0]);n=series.shape[1]
    if labels is not None and np.shape(labels)!=(33,n):raise ValueError('Wrong replay label geometry')
    white=np.array([transform.conj().T@a@transform for a in series]);reference=white[tref]
    values=np.full((33,n),complex(np.nan,np.nan));vectors=np.full((33,n,n),complex(np.nan,np.nan))
    residuals=np.full((33,n),np.nan);temporal=np.full((33,n),np.nan);cross=np.ones((33,n));previous=None
    for t in range(tref+1,33):
        eigenvalues,u=linalg.eig(white[t],reference)
        if not np.all(np.isfinite(eigenvalues)):raise ValueError('Singular generalized pencil')
        v=transform@u
        if labels is not None:order=match_spectrum(labels[t],eigenvalues)
        elif central is not None:
            order,quality=match(central['vectors'][t],v,central_metric)
        elif previous is not None:order,quality=match(previous,v,series[0])
        else:order=np.lexsort((eigenvalues.imag,-eigenvalues.real));quality=np.ones(n)
        eigenvalues=eigenvalues[order];v=v[:,order];u=u[:,order]
        norms=np.sqrt(np.real(np.einsum('ij,ij->j',v.conj(),series[0]@v)))
        v=v/norms;u=u/norms
        for j in range(n):
            pivot=int(np.argmax(np.abs(v[:,j])));rotation=np.exp(-1j*np.angle(v[pivot,j]));v[:,j]*=rotation;u[:,j]*=rotation
        if central is not None:cross[t]=np.diag(overlaps(central['vectors'][t],v,central_metric)) if labels is not None else quality
        else:cross[t]=np.ones(n)
        temporal[t]=np.diag(overlaps(previous,v,series[0])) if previous is not None else np.ones(n)
        for j in range(n):
            denominator=(linalg.norm(white[t])+abs(eigenvalues[j])*linalg.norm(reference))*linalg.norm(u[:,j])
            residuals[t,j]=linalg.norm(white[t]@u[:,j]-eigenvalues[j]*reference@u[:,j])/denominator if denominator else 0.
        values[t]=eigenvalues;vectors[t]=v;previous=v
    return dict(values=values,vectors=vectors,residuals=residuals,temporal_overlap=temporal,central_overlap=cross)


def energy(first,second,phase_bound):
    if not np.isfinite(first) or not np.isfinite(second) or first==0 or second==0:return None,False
    result=np.log(first/second)
    valid=first.real>0 and second.real>0 and abs(np.angle(first))<=phase_bound and abs(np.angle(second))<=phase_bound
    return complex(result),bool(valid)


def compute(rows,validation,mode_labels=None):
    plan=json.loads((HERE/'ANALYSIS_PLAN.json').read_text());data,control=assemble(rows,plan)
    configs=plan['configuration_labels'];cuts=plan['cutoffs'];radii=plan['radii'];phase=plan['display_phase_bound_radians']
    indices={name:[cuts.index(n) for n in selected] for name,selected in plan['bases'].items()}
    mean=data.mean(axis=0)
    samples={'mean':mean};labels={'mean':configs}
    for i,cfg in enumerate(configs):
        samples['cfg%d'%cfg]=data[i];labels['cfg%d'%cfg]=[cfg]
        samples['drop%d'%cfg]=np.delete(data,i,axis=0).mean(axis=0);labels['drop%d'%cfg]=[c for c in configs if c!=cfg]
    diagnostics=[];metric_spectrum=[];modes=[];fixed=[];coefficients=[];hermitian=[];controls=[];leading=[];solutions={}
    for sample,matrices in samples.items():
        for ri,radius in enumerate(radii):
            for basis,selected in indices.items():
                series=subset(matrices[ri],selected);transform,d,eig=metric_transform(series[0]);normal=series/d[None,:,None]/d[None,None,:]
                for j,value in enumerate(eig):metric_spectrum.append(dict(sample=sample,R=radius,basis=basis,mode=j,eigenvalue=float(value),condition=float(eig[-1]/eig[0]),configurations=';'.join(map(str,labels[sample]))))
                for t,a in enumerate(normal):
                    singular=linalg.svdvals(a);norm=linalg.norm(a);part=(a+a.conj().T)/2
                    diagnostics.append(dict(sample=sample,R=radius,basis=basis,tau=t,dimension=len(selected),
                        singular_min=float(singular[-1]),singular_max=float(singular[0]),condition=float(singular[0]/singular[-1]),
                        relative_antihermitian_norm=float(linalg.norm(a-a.conj().T)/2/norm),relative_imaginary_norm=float(linalg.norm(a.imag)/norm),
                        hermitian_part_min_eigenvalue=float(linalg.eigvalsh(part)[0])))
    for ri,radius in enumerate(radii):
        for basis,selected in indices.items():
            central_series=subset(mean[ri],selected)
            for tref in plan['reference_times']:
                replay_labels=mode_labels['mean',radius,basis,tref] if mode_labels is not None else None
                central=solve_series(central_series,tref,labels=replay_labels);solutions['mean',radius,basis,tref]=central
                for sample in ['mean']+(['drop%d'%c for c in configs] if tref==0 else []):
                    if sample=='mean':result=central
                    else:
                        replay_labels=mode_labels[sample,radius,basis,tref] if mode_labels is not None else None
                        result=solve_series(subset(samples[sample][ri],selected),tref,central,central_series[0],labels=replay_labels);solutions[sample,radius,basis,tref]=result
                    for t in range(tref+1,33):
                        for j,z in enumerate(result['values'][t]):
                            effective,valid=energy(z,result['values'][t+1,j],phase) if t<32 else (None,False)
                            modes.append(dict(sample=sample,R=radius,basis=basis,t_ref=tref,tau=t,mode=j,lambda_re=float(z.real),lambda_im=float(z.imag),phase=float(np.angle(z)),
                                relative_residual=float(result['residuals'][t,j]),temporal_overlap=float(result['temporal_overlap'][t,j]),central_overlap=float(result['central_overlap'][t,j]),
                                effective_energy_re=effective.real if effective is not None else '',effective_energy_im=effective.imag if effective is not None else '',energy_display_valid=int(valid)))
                    failures=[t for t in range(tref+1,33) if result['values'][t,0].real<=0 or abs(np.angle(result['values'][t,0]))>phase]
                    leading.append(dict(sample=sample,R=radius,basis=basis,t_ref=tref,first_leading_phase_failure=min(failures) if failures else '',
                        maximum_relative_residual=float(np.nanmax(result['residuals'])),minimum_early_temporal_overlap=float(np.nanmin(result['temporal_overlap'][tref+1:9,0])),
                        minimum_early_central_overlap=float(np.nanmin(result['central_overlap'][tref+1:9,0]))))
            transform,_,_=metric_transform(central_series[0])
            for t in range(1,33):
                a=transform.conj().T@central_series[t]@transform
                hvalues=linalg.eigvalsh((a+a.conj().T)/2)[::-1]
                for j,z in enumerate(hvalues):hermitian.append(dict(R=radius,basis=basis,tau=t,ordered_mode=j,hermitian_part_lambda=float(z),interpretation='diagnostic_only'))
        for t in range(33):
            p10=data[configs.index(650),ri,t];p40=np.array([[control[radius,t,k,s] for s in cuts] for k in cuts]);difference=p40-p10
            controls.append(dict(R=radius,tau=t,maximum_component_difference=float(max(np.max(abs(difference.real)),np.max(abs(difference.imag)))),relative_frobenius_difference=float(linalg.norm(difference)/linalg.norm(p10))))
        projections=[]
        for n in (40,100):
            vector=np.zeros(4,complex);vector[cuts.index(n)]=1/np.sqrt(mean[ri,0,cuts.index(n),cuts.index(n)].real)
            projections.append(('Nv%d'%n,'mean',vector,list(range(3)),None))
        for basis,selected in indices.items():
            t=plan['fixed_vector_optimization_time']
            for sample in ['mean']+['drop%d'%c for c in configs]:
                small=solutions[sample,radius,basis,0]['vectors'][t,:,0];vector=np.zeros(4,complex);vector[selected]=small
                training=[configs.index(c) for c in labels[sample]];heldout=next((i for i in range(3) if i not in training),None)
                projections.append(('fixed_'+basis,sample,vector,training,heldout))
                for i,z in enumerate(vector):coefficients.append(dict(sample=sample,R=radius,basis=basis,optimization_tau=t,Nv=cuts[i],coefficient_re=float(z.real),coefficient_im=float(z.imag)))
        for model,sample,vector,training,heldout in projections:
            curve=np.einsum('i,ntij,j->nt',vector.conj(),data[:,ri],vector)
            for evaluation,selected in [('training',training)]+([('heldout',[heldout])] if heldout is not None else []):
                q=curve[selected].mean(axis=0)
                for t,z in enumerate(q):
                    effective,valid=energy(z,q[t+1],phase) if t<32 else (None,False)
                    fixed.append(dict(sample=sample,R=radius,model=model,evaluation=evaluation,training_configurations=';'.join(str(configs[i]) for i in training),
                        evaluation_configurations=';'.join(str(configs[i]) for i in selected),tau=t,correlator_re=float(z.real),correlator_im=float(z.imag),phase=float(np.angle(z)),
                        effective_energy_re=effective.real if effective is not None else '',effective_energy_im=effective.imag if effective is not None else '',energy_display_valid=int(valid)))
    maximum_residual=max(row['relative_residual'] for row in modes)
    if maximum_residual>plan['numerical_relative_residual_limit']:raise ValueError('Generalized eigenpair residual exceeds numerical limit')
    features=data.reshape(3,-1);features=np.concatenate([features.real,features.imag],axis=1);centered=features-features.mean(axis=0)
    rank=int(np.linalg.matrix_rank(centered));conditions={str(r):{name:next(row['condition'] for row in metric_spectrum if row['sample']=='mean' and row['R']==r and row['basis']==name) for name in indices} for r in radii}
    summary=dict(technical_pilot=True,production_submitted=False,configuration_labels=configs,statistical_sample_size=3,independent_unit='configuration',
        input_matrix_rows=len(rows),analysis_axis=-1,bases=plan['bases'],equal_time_diagonally_normalized_conditions=conditions,
        empirical_covariance_rank=rank,empirical_covariance_rank_bound=2,maximum_generalized_eigenpair_relative_residual=maximum_residual,
        maximum_p10_p40_relative_frobenius_difference=max(r['relative_frobenius_difference'] for r in controls),
        raw_complex_matrices_retained=True,nonzero_time_hermitian_symmetrization_used_for_primary_results=False,
        reference_times=plan['reference_times'],fixed_vector_optimization_time=plan['fixed_vector_optimization_time'],display_phase_bound_radians=phase,
        validation=validation,output_rows=dict(matrix_diagnostics=len(diagnostics),metric_spectrum=len(metric_spectrum),complex_modes=len(modes),fixed_projections=len(fixed),fixed_vectors=len(coefficients),projection_controls=len(controls),hermitian_part_diagnostics=len(hermitian),leading_mode_diagnostics=len(leading)),
        conclusions=['All three baseline diagonal regressions are exact; the n650 full-matrix projection control agrees within roundoff.',
            'The equal-time basis is full rank; nested cutoffs have strong short-distance correlations. Smaller bases provide conditioning controls.',
            'Finite-time matrices and generalized eigenvalues retain their complex nonsymmetric form. Long-distance and late-time sensitivity requires full-ensemble data.',
            'Three configurations support technical and leave-one-out sensitivity diagnostics. No final energies, confidence intervals, correlated fits or string tension are inferred.',
            'Next: implement full71 production with validated pilot reuse; submissions are deferred.'])
    tables=dict(matrix_diagnostics=diagnostics,metric_spectrum=metric_spectrum,complex_modes=modes,fixed_projections=fixed,fixed_vectors=coefficients,projection_controls=controls,hermitian_part_diagnostics=hermitian,leading_mode_diagnostics=leading)
    return tables,summary
