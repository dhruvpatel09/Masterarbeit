"""Independent-configuration statistics and small correlated linear fits."""
import math
from decimal import Decimal as D,localcontext


def mean_se(values):
    with localcontext() as ctx:
        ctx.prec=50;n=len(values)
        if n<2:raise ValueError('At least two matched configurations required')
        mean=sum(values,D(0))/n
        return mean,(sum(((v-mean)**2 for v in values),D(0))/(n*(n-1))).sqrt()


def covariance(x,y):
    if len(x)!=len(y):raise ValueError('Unmatched configurations')
    with localcontext() as ctx:
        ctx.prec=50;n=len(x);mx=sum(x,D(0))/n;my=sum(y,D(0))/n
        return sum(((a-mx)*(b-my) for a,b in zip(x,y)),D(0))/(n*(n-1))


def energy(x,y):
    if len(x)!=len(y):raise ValueError('Unmatched configurations')
    with localcontext() as ctx:
        ctx.prec=50;n=len(x);mx,sx=mean_se(x);my,sy=mean_se(y)
        out=dict(mean_L_tau=mx,mean_L_next=my,mean_L_tau_se=sx,mean_L_next_se=sy,
            valid_point=mx>0 and my>0,aE_eff='',naive_se='',valid_jackknife=False,jackknife_se='')
        if not out['valid_point']:return out,None,None
        out['aE_eff']=(mx/my).ln()
        delta=[(a-mx)/mx-(b-my)/my for a,b in zip(x,y)];out['naive_se']=mean_se(delta)[1]
        ax,ay=sum(x,D(0)),sum(y,D(0));loo=[]
        if all(ax-a>0 and ay-b>0 for a,b in zip(x,y)):
            loo=[((ax-a)/(ay-b)).ln() for a,b in zip(x,y)]
            avg=sum(loo,D(0))/n
            out['valid_jackknife']=True;out['jackknife_se']=(D(n-1)/n*sum(((v-avg)**2 for v in loo),D(0))).sqrt()
        return out,delta,loo


def chi2_sf(q,df):
    """Upper chi-square tail for positive integer degrees of freedom."""
    if q<0 or df<1 or not math.isfinite(q):raise ValueError('Invalid chi-square arguments')
    z=q/2
    if z==0:return 1.
    if df%2==0:
        logs=[-z+k*math.log(z)-math.lgamma(k+1) for k in range(df//2)]
        maximum=max(logs);result=math.exp(maximum)*sum(math.exp(v-maximum) for v in logs)
    else:
        result=math.erfc(math.sqrt(z))
        for k in range((df-1)//2):
            shape=k+.5;result+=math.exp(shape*math.log(z)-z-math.lgamma(shape+1))
    return min(1.,max(0.,result))


def gls(values,fluctuations,design):
    """Fixed full-sample generalized least squares, preserving paired samples."""
    import numpy as np
    y=np.asarray(values,float);d=np.asarray(fluctuations,float);A=np.asarray(design,float)
    n,m=d.shape
    if y.shape!=(m,) or A.shape[0]!=m or n<2:raise ValueError('Invalid fit geometry')
    C=np.cov(d,rowvar=False,ddof=1)/n
    eigen=np.linalg.eigvalsh(C)
    if eigen[0]<=0:raise ValueError('Covariance is not positive definite')
    np.linalg.cholesky(C)
    invA=np.linalg.solve(C,A);normal=A.T@invA
    W=np.linalg.solve(normal,invA.T);params=W@y
    pd=d@W.T;cov=np.cov(pd,rowvar=False,ddof=1)/n
    cov=np.atleast_2d(cov);res=y-A@params;q=float(res@np.linalg.solve(C,res));dof=m-A.shape[1]
    return dict(params=params,weights=W,fluctuations=pd,covariance=cov,chi2=q,dof=dof,
        p_gaussian_approx=chi2_sf(q,dof),condition_inf=float(np.linalg.cond(C,np.inf)),minimum_eigenvalue=float(eigen[0]))
