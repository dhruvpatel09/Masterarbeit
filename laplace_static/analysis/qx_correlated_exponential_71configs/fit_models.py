"""Bounded forward exponentials with a full fixed covariance of the mean."""
import math
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2,f

ENERGY_CEILING=10.


def model(p,t,states):
    if states==1:return p[0]*np.exp(-p[1]*t)
    return p[0]*np.exp(-p[2]*t)+p[1]*np.exp(-(p[2]+p[3])*t)


def derivatives(p,t,states):
    """Analytic model Jacobian and Hessian at each time."""
    m=len(t)
    if states==1:
        A,E=p;e=np.exp(-E*t);J=np.column_stack([e,-A*t*e]);H=np.zeros((m,2,2))
        H[:,0,1]=H[:,1,0]=-t*e;H[:,1,1]=A*t*t*e
    else:
        A,B,E,gap=p;e=np.exp(-E*t);d=np.exp(-(E+gap)*t)
        J=np.column_stack([e,d,-t*(A*e+B*d),-t*B*d]);H=np.zeros((m,4,4))
        H[:,0,2]=H[:,2,0]=-t*e
        H[:,1,2]=H[:,2,1]=-t*d;H[:,1,3]=H[:,3,1]=-t*d
        H[:,2,2]=t*t*(A*e+B*d)
        H[:,2,3]=H[:,3,2]=t*t*B*d;H[:,3,3]=t*t*B*d
    return J,H


class Window:
    def __init__(self,samples,times,scale=None):
        raw=np.asarray(samples,float);self.times=np.asarray(times,float)
        if raw.ndim!=2 or raw.shape[1]!=len(self.times) or raw.shape[0]<2 or not np.all(np.isfinite(raw)):
            raise ValueError('Invalid matched configuration samples')
        self.n=raw.shape[0];self.t=self.times-self.times[0]
        self.scale=float(np.max(np.abs(raw.mean(0)))) if scale is None else float(scale)
        if not math.isfinite(self.scale) or self.scale<=0:raise ValueError('No finite correlator normalization')
        self.samples=raw/self.scale;self.mean=self.samples.mean(0)
        self.covariance=np.cov(self.samples,rowvar=False,ddof=1)/self.n
        self.se=np.sqrt(np.diag(self.covariance))
        if not np.all(self.se>0):raise ValueError('Zero or invalid covariance diagonal')
        self.correlation=self.covariance/self.se[:,None]/self.se[None,:]
        eigen=np.linalg.eigvalsh(self.correlation)
        if eigen[0]<=0:raise ValueError('Covariance is not positive definite; no modes removed')
        self.condition=float(eigen[-1]/eigen[0]);L=np.linalg.cholesky(self.correlation)
        self.W=solve_triangular(L,np.diag(1/self.se),lower=True)

    def residual(self,p,states,mean=None):
        return self.W@(model(p,self.t,states)-(self.mean if mean is None else mean))

    def jacobian(self,p,states):return self.W@derivatives(p,self.t,states)[0]

    def local_errors(self,p,states,mean=None):
        """Implicit derivative of the actual optimum, including residual curvature."""
        J,H=derivatives(p,self.t,states);Jw=self.W@J;r=self.residual(p,states,mean)
        hessian=Jw.T@Jw+np.einsum('t,tij->ij',self.W.T@r,H)
        norms=np.linalg.norm(Jw,axis=0)
        if np.any(norms<=0):raise ValueError('Unidentified model Jacobian')
        singular=np.linalg.svd(Jw/norms,compute_uv=False)
        scaled=hessian/norms[:,None]/norms[None,:]
        eigen=np.linalg.eigvalsh(scaled)
        if singular[-1]<=1e-10*singular[0] or eigen[0]<=1e-10*eigen[-1]:
            raise ValueError('Unidentified or non-positive local curvature')
        G=np.linalg.solve(hessian,Jw.T@self.W)
        covariance=G@self.covariance@G.T
        if not np.all(np.diag(covariance)>=0):raise ValueError('Invalid local covariance')
        return dict(se=np.sqrt(np.diag(covariance)),covariance=covariance,gradient=G,
            fluctuations=(self.samples-self.mean)@G.T,
            scaled_jacobian_condition=float(singular[0]/singular[-1]),
            scaled_hessian_condition=float(eigen[-1]/eigen[0]))


def bounds(states):
    return (np.zeros(2),np.array([np.inf,ENERGY_CEILING])) if states==1 else (
        np.zeros(4),np.array([np.inf,np.inf,ENERGY_CEILING,ENERGY_CEILING]))


def initial_guesses(states,single=None):
    if states==1:return [np.array([1.,E]) for E in (.02,.1,.25,.5,1.)]
    A,E=single
    result=[np.array([A,1e-7,max(E,1e-5),.2])]
    for ratio in (.5,.9):
        for fraction in (.2,.8):
            for gap in (.05,.2,.6):result.append(np.array([fraction*A,(1-fraction)*A,max(E*ratio,1e-5),gap]))
    return result


def solve(window,states,starts=None,mean=None,single=None):
    starts=initial_guesses(states,single) if starts is None else starts
    candidates=[]
    for i,start in enumerate(starts):
        try:
            result=least_squares(lambda p:window.residual(p,states,mean),start,
                jac=lambda p:window.jacobian(p,states),bounds=bounds(states),method='trf',
                x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=1500)
            candidates.append(dict(start=i,p=result.x,success=bool(result.success),
                q2=float(result.fun@result.fun),nfev=int(result.nfev),optimality=float(result.optimality)))
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as exc:
            candidates.append(dict(start=i,p=None,success=False,q2=None,nfev=0,error=str(exc)))
    valid=[c for c in candidates if c['success'] and c['p'] is not None and math.isfinite(c['q2'])]
    if not valid:return None,candidates
    return min(valid,key=lambda c:c['q2']),candidates


def diagnostics(window,states,p,q2,mean=None,candidates=()):
    flags=[];lower,upper=bounds(states);p=np.asarray(p,float)
    if np.any(p<=1e-7) or np.any(np.isfinite(upper)&(p>=upper-1e-5)):
        flags.append('parameter_boundary')
    local=None
    try:local=window.local_errors(p,states,mean)
    except (ValueError,np.linalg.LinAlgError):flags.append('unidentified_curvature')
    if states==2:
        fraction=float(p[0]/(p[0]+p[1])) if p[0]+p[1] else 0.
        if min(fraction,1-fraction)<1e-4:flags.append('vanishing_state_amplitude')
        if p[3]<1e-4:flags.append('nearly_equal_energies')
        if local is not None:
            if local['se'][0]>=p[0]:flags.append('ground_amplitude_unresolved')
            if local['se'][1]>=p[1]:flags.append('excited_amplitude_unresolved')
            if local['se'][3]>=p[3]:flags.append('gap_unresolved')
    else:fraction=1.
    near=[float(c['p'][1 if states==1 else 2]) for c in candidates
          if c['success'] and c['p'] is not None and c['q2']<=q2+1.]
    span=max(near)-min(near) if near else 0.
    eindex=1 if states==1 else 2
    if local is not None and span>max(.02,5*local['se'][eindex]):flags.append('different_near_optimal_starts')
    npar=2 if states==1 else 4;dof=len(window.t)-npar
    # Hotelling reference uses the linearized interior-model residual dimension.
    hotelling=None
    if not flags and window.n>dof:
        statistic=(window.n-dof)*q2/(dof*(window.n-1))
        hotelling=float(f.sf(statistic,dof,window.n-dof))
    return dict(flags=flags,local=local,identified=not flags,ground_fraction=fraction,
        near_optimal_E0_span=span,p_gaussian_reference=float(chi2.sf(q2,dof)),
        p_hotelling_linearized_reference=hotelling,dof=dof)
