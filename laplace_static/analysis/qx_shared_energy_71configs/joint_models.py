"""Joint source-operator fits retaining the complete source/time covariance."""
import importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'qx_correlated_exponential_71configs/fit_models.py'
spec=importlib.util.spec_from_file_location('qx_published_fit_models',str(BASE))
published=importlib.util.module_from_spec(spec);spec.loader.exec_module(published)
KINDS=('shared_one','independent_one','shared_two')


def model(p,t,m,kind):
    p=np.asarray(p);t=np.asarray(t)
    if kind=='shared_one':return (p[:m,None]*np.exp(-p[m]*t)).reshape(-1)
    if kind=='independent_one':return (p[:m,None]*np.exp(-p[m:,None]*t)).reshape(-1)
    if kind=='shared_two':return (p[:m,None]*np.exp(-p[2*m]*t)+p[m:2*m,None]*np.exp(-(p[2*m]+p[2*m+1])*t)).reshape(-1)
    raise ValueError('Unknown model')


def derivatives(p,t,m,kind):
    p=np.asarray(p);t=np.asarray(t)
    nt=len(t);d=m+1 if kind=='shared_one' else 2*m if kind=='independent_one' else 2*m+2
    J=np.zeros((m*nt,d));H=np.zeros((m*nt,d,d))
    for i in range(m):
        sl=slice(i*nt,(i+1)*nt)
        if kind!='shared_two':
            ei=m if kind=='shared_one' else m+i;e=np.exp(-p[ei]*t);A=p[i]
            J[sl,i]=e;J[sl,ei]=-A*t*e
            H[sl,i,ei]=H[sl,ei,i]=-t*e;H[sl,ei,ei]=A*t*t*e
        else:
            ei=2*m;gi=ei+1;ai=i;bi=m+i;A,B=p[ai],p[bi]
            e=np.exp(-p[ei]*t);b=np.exp(-(p[ei]+p[gi])*t)
            J[sl,ai]=e;J[sl,bi]=b;J[sl,ei]=-t*(A*e+B*b);J[sl,gi]=-t*B*b
            H[sl,ai,ei]=H[sl,ei,ai]=-t*e
            H[sl,bi,ei]=H[sl,ei,bi]=-t*b;H[sl,bi,gi]=H[sl,gi,bi]=-t*b
            H[sl,ei,ei]=t*t*(A*e+B*b)
            H[sl,ei,gi]=H[sl,gi,ei]=t*t*B*b;H[sl,gi,gi]=t*t*B*b
    return J,H


def bounds(m,kind):
    n=m+1 if kind=='shared_one' else 2*m if kind=='independent_one' else 2*m+2
    upper=np.full(n,np.inf)
    upper[m:]=10. if kind!='shared_two' else np.inf
    if kind=='shared_two':upper[-2:]=10.
    return np.zeros(n),upper


class JointWindow(published.Window):
    def __init__(self,samples,times,cutoffs,scales=None):
        raw=np.asarray(samples,float);self.source_times=np.asarray(times,float);self.cutoffs=list(cutoffs)
        self.m=len(cutoffs);self.nt=len(times)
        if raw.ndim!=3 or raw.shape[1:]!=(self.m,self.nt):raise ValueError('Source/time shape mismatch')
        self.source_scales=np.max(np.abs(raw.mean(0)),axis=1) if scales is None else np.asarray(scales,float)
        if not np.all(self.source_scales>0):raise ValueError('Zero source normalization')
        normalized=raw/self.source_scales[None,:,None]
        super().__init__(normalized.reshape(raw.shape[0],-1),np.tile(times,self.m),scale=1.)
        self.dt=self.source_times-self.source_times[0]

    def residual(self,p,kind,mean=None):
        return self.W@(model(p,self.dt,self.m,kind)-(self.mean if mean is None else mean))

    def jacobian(self,p,kind):return self.W@derivatives(p,self.dt,self.m,kind)[0]

    def local_errors(self,p,kind,mean=None):
        J,H=derivatives(p,self.dt,self.m,kind);Jw=self.W@J;r=self.residual(p,kind,mean)
        curvature=Jw.T@Jw+np.einsum('t,tij->ij',self.W.T@r,H)
        norms=np.linalg.norm(Jw,axis=0)
        if np.any(norms<=0):raise ValueError('Unidentified model derivative')
        singular=np.linalg.svd(Jw/norms,compute_uv=False)
        eigen=np.linalg.eigvalsh(curvature/norms[:,None]/norms[None,:])
        if singular[-1]<=1e-10*singular[0] or eigen[0]<=1e-10*eigen[-1]:raise ValueError('Unidentified or non-positive local curvature')
        G=np.linalg.solve(curvature,Jw.T@self.W);cov=G@self.covariance@G.T
        if np.any(np.diag(cov)<0):raise ValueError('Invalid parameter covariance')
        return dict(se=np.sqrt(np.diag(cov)),covariance=cov,gradient=G,fluctuations=(self.samples-self.mean)@G.T,
            jacobian_condition=float(singular[0]/singular[-1]),curvature_condition=float(eigen[-1]/eigen[0]))


def initial_guesses(w,kind,shared=None,independent=None,seed_energies=None):
    m=w.m
    if kind=='shared_one':return [np.r_[np.ones(m),e] for e in (.06,.12,.22,.35,.55)]
    if kind=='independent_one':
        p=shared['p'];starts=[np.r_[p[:m],np.repeat(p[m],m)]]
        if seed_energies is not None:starts.append(np.r_[np.ones(m),seed_energies])
        for e in (.06,.15,.3,.55):starts.append(np.r_[np.ones(m),np.repeat(e,m)])
        return starts
    p=np.asarray(shared['p']);A,E=p[:m],p[m];starts=[np.r_[A,np.repeat(1e-7,m),max(E,1e-5),.2]]
    for ratio in (.5,.9):
        for fraction in (.2,.8):
            for gap in (.05,.2,.6):starts.append(np.r_[fraction*A,(1-fraction)*A,max(E*ratio,1e-5),gap])
    if independent is not None:
        e=np.asarray(independent['p'][m:]);lo=max(float(np.min(e))*.9,1e-5);gap=max(float(np.max(e)-lo),.05)
        excited=np.clip((e-lo)/gap,.1,.9);amps=np.asarray(independent['p'][:m])
        starts.append(np.r_[amps*(1-excited),amps*excited,lo,gap])
    return starts


def solve(w,kind,starts=None,mean=None,shared=None,independent=None,seed_energies=None):
    if starts is None:starts=initial_guesses(w,kind,shared,independent,seed_energies)
    candidates=[]
    for i,p in enumerate(starts):
        try:
            r=least_squares(lambda v:w.residual(v,kind,mean),p,jac=lambda v:w.jacobian(v,kind),
                bounds=bounds(w.m,kind),method='trf',x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=1800)
            candidates.append(dict(start=i,p=r.x.tolist(),success=bool(r.success),q2=float(r.fun@r.fun),nfev=int(r.nfev),optimality=float(r.optimality)))
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as exc:
            candidates.append(dict(start=i,p=None,success=False,q2=None,nfev=0,error=str(exc)))
    valid=[c for c in candidates if c['success'] and c['p'] is not None and np.isfinite(c['q2'])]
    return (min(valid,key=lambda c:c['q2']) if valid else None),candidates


def diagnostics(w,kind,p,q2,mean=None,candidates=()):
    p=np.asarray(p);lower,upper=bounds(w.m,kind);flags=[]
    if np.any(p<=1e-7) or np.any(np.isfinite(upper)&(p>=upper-1e-5)):flags.append('parameter_boundary')
    try:local=w.local_errors(p,kind,mean)
    except (ValueError,np.linalg.LinAlgError):local=None;flags.append('unidentified_curvature')
    ei=2*w.m if kind=='shared_two' else w.m
    if kind=='shared_two':
        fractions=p[:w.m]/(p[:w.m]+p[w.m:2*w.m])
        if np.any(np.minimum(fractions,1-fractions)<1e-4):flags.append('vanishing_source_amplitude')
        if p[-1]<1e-4:flags.append('nearly_equal_energies')
        if local is not None:
            if np.any(local['se'][:2*w.m]>=p[:2*w.m]):flags.append('source_amplitude_unresolved')
            if local['se'][-1]>=p[-1]:flags.append('gap_unresolved')
    energy_indices=list(range(w.m,2*w.m)) if kind=='independent_one' else [ei]
    near=[np.array(c['p'])[energy_indices] for c in candidates if c['success'] and c['p'] is not None and c['q2']<=q2+1.]
    span=float(np.max(np.ptp(np.array(near),axis=0))) if near else 0.
    reference_se=float(np.max(local['se'][energy_indices])) if local is not None else None
    if local is not None and span>max(.02,5*reference_se):flags.append('different_near_optimal_starts')
    dof=len(w.mean)-len(p)
    return dict(local=local,flags=flags,local_checks_pass=not flags,near_E_span=span,dof=dof,
        p_gaussian_reference=float(chi2.sf(q2,dof)))
