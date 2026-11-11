"""Signed-amplitude one-exponential fits with fixed full-sample covariance."""
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import brentq


class SingularCovariance(ValueError):pass


def whitening(covariance):
    cov=np.asarray(covariance,float)
    if not np.isfinite(cov).all() or np.any(np.diag(cov)<=0):raise SingularCovariance('Nonpositive variance')
    d=np.sqrt(np.diag(cov));corr=cov/d[:,None]/d[None,:]
    eig=np.linalg.eigvalsh((corr+corr.T)/2)
    if eig[0]<=1e-12*eig[-1]:raise SingularCovariance('Covariance below rank threshold')
    try:L=np.linalg.cholesky(corr)
    except np.linalg.LinAlgError:raise SingularCovariance('Covariance Cholesky failed')
    W=solve_triangular(L,np.diag(1/d),lower=True)
    return W,float(eig[-1]/eig[0])


class ExponentialFit:
    def __init__(self,covariance,times,limits=(-.5,1.5),grid_points=161):
        self.W,self.condition=whitening(covariance)
        self.t=np.asarray(times,float)-times[0]
        self.limits=limits;self.grid=np.linspace(limits[0],limits[1],grid_points)
        self.F=np.exp(-self.t[:,None]*self.grid[None,:]);self.G=self.W.dot(self.F)
        self.Gp=self.W.dot(-self.t[:,None]*self.F)
        self.b=np.sum(self.G*self.G,axis=0);self.bp=2*np.sum(self.Gp*self.G,axis=0)

    def evaluate(self,E,y):
        f=np.exp(-E*self.t);g=self.W.dot(f);gp=self.W.dot(-self.t*f);gpp=self.W.dot(self.t*self.t*f)
        a=np.dot(g,y);ap=np.dot(gp,y);app=np.dot(gpp,y)
        b=np.dot(g,g);bp=2*np.dot(gp,g);bpp=2*(np.dot(gpp,g)+np.dot(gp,gp))
        amplitude=a/b;res=y-amplitude*g;q=float(np.dot(res,res))
        grad=-2*amplitude*np.dot(gp,res)
        curvature=-2*(ap*ap+a*app)/b+4*a*ap*bp/(b*b)+a*a*bpp/(b*b)-2*a*a*bp*bp/(b*b*b)
        return q,float(amplitude),float(grad),float(curvature)

    def fit_many(self,observations):
        data=np.atleast_2d(np.asarray(observations,float))
        if data.shape[1]!=len(self.t) or not np.isfinite(data).all():raise ValueError('Fit observation inventory')
        Y=data.dot(self.W.T);a=Y.dot(self.G);ap=Y.dot(self.Gp)
        grad=-2*a*ap/self.b+a*a*self.bp/(self.b*self.b)
        results=[]
        for y,g in zip(Y,grad):
            candidates=[self.limits[0],self.limits[1]]
            for i in range(len(self.grid)-1):
                if g[i]<=0 and g[i+1]>=0:
                    if g[i]==0:candidates.append(float(self.grid[i]))
                    elif g[i+1]==0:candidates.append(float(self.grid[i+1]))
                    else:
                        E=float(brentq(lambda E:self.evaluate(E,y)[2],self.grid[i],self.grid[i+1],xtol=5e-13,rtol=1e-14))
                        for _ in range(2):
                            _,_,derivative,curvature=self.evaluate(E,y)
                            trial=E-derivative/curvature if curvature>0 else E
                            if self.grid[i]<=trial<=self.grid[i+1]:E=float(trial)
                        candidates.append(E)
            values=[self.evaluate(E,y) for E in candidates];which=int(np.argmin([v[0] for v in values]))
            E=candidates[which];q,A,gradient,curvature=values[which]
            boundary=E==self.limits[0] or E==self.limits[1]
            scaled_gradient=abs(gradient)/(1+abs(q))
            energy_residual=abs(gradient/curvature) if curvature>0 else float('inf')
            if not np.isfinite([E,q,A,gradient,curvature]).all():raise ValueError('Nonfinite exponential fit')
            if not boundary and energy_residual>1e-10:raise ValueError('Profile stationarity failed')
            results.append(dict(energy=E,amplitude=A,q=q,boundary=boundary,curvature=curvature,
                scaled_gradient=scaled_gradient,energy_residual=energy_residual,candidates=len(candidates)))
        return results
