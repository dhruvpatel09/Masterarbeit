"""Paired configuration statistics and unregularized correlated fits."""
import math
from decimal import Decimal as D, localcontext

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2


def mean_se(x):
    with localcontext() as ctx:
        ctx.prec = 50
        n = len(x)
        m = sum(x, D(0)) / n
        return m, (sum(((a-m)**2 for a in x), D(0))/(n*(n-1))).sqrt()


def energy(x, y, jackknife=True):
    """Log of ensemble means; retain all signed configuration observations."""
    with localcontext() as ctx:
        ctx.prec = 50
        n = len(x)
        mx, _ = mean_se(x)
        my, _ = mean_se(y)
        if mx <= 0 or my <= 0:
            return None, None, None, 0
        value = (mx/my).ln()
        delta = [(a-mx)/mx-(b-my)/my for a, b in zip(x, y)]
        error = mean_se(delta)[1]
        sx, sy = sum(x, D(0)), sum(y, D(0))
        positive = sum(sx-a > 0 and sy-b > 0 for a, b in zip(x, y))
        loo = None
        if jackknife and positive == n:
            loo = np.array([float(((sx-a)/(sy-b)).ln()) for a, b in zip(x, y)])
        return (str(value), str(error)), np.array(delta, float), loo, positive


def covariance(delta):
    d = np.asarray(delta, float)
    if d.ndim != 2 or not np.isfinite(d).all():
        raise ValueError('Nonfinite or malformed configuration fluctuations')
    return np.atleast_2d(np.cov(d, rowvar=False, ddof=1)/d.shape[0])


def standard_error(delta):
    return float(np.std(delta, ddof=1)/math.sqrt(len(delta)))


def jackknife_error(samples):
    if samples is None:
        return None
    x = np.asarray(samples, float)
    return float(np.sqrt((len(x)-1.)/len(x)*np.sum((x-x.mean())**2)))


def whiten(C):
    """Strictly positive definite covariance. No pseudoinverse or ridge."""
    C = np.asarray(C, float)
    if not np.isfinite(C).all() or not np.allclose(C, C.T, rtol=1e-12, atol=0):
        raise ValueError('Invalid covariance')
    scale = np.sqrt(np.diag(C))
    if np.any(scale <= 0):
        raise ValueError('Nonpositive variance')
    corr = C/np.outer(scale, scale)
    factor = cho_factor(corr, lower=True, check_finite=True)
    def solve(b):
        b = np.asarray(b, float)
        return (cho_solve(factor, b/scale)/scale if b.ndim == 1 else
                cho_solve(factor, b/scale[:, None])/scale[:, None])
    return solve


def constant_fit(y, delta):
    y = np.asarray(y, float)
    if not np.isfinite(y).all():
        raise ValueError('Nonfinite fit values')
    C = covariance(delta)
    solve = whiten(C)
    one = np.ones(len(y))
    weight = solve(one)
    weight /= one@weight
    estimate = float(weight@y)
    residual = y-estimate
    q = float(residual@solve(residual))
    return dict(estimate=estimate, naive_se=standard_error(delta@weight),
                chi2=q, dof=len(y)-1,
                p_gaussian_approx=float(chi2.sf(q, len(y)-1)),
                condition_inf=float(np.linalg.cond(C, np.inf))), weight


def assert_close(a, b, name, rtol=2e-9, atol=2e-12):
    if not np.isfinite([a, b]).all() or abs(a-b) > atol+rtol*max(abs(a), abs(b)):
        raise ValueError('Numerical mismatch: {} ({!r}, {!r})'.format(name, a, b))
