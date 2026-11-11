"""Configuration-level independent errors and paired delete-one statistics."""
import csv
import json
from pathlib import Path
import numpy as np


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def independent_se(values):
    return np.std(values, axis=0, ddof=1)/np.sqrt(len(values))


def jk_cov(replicates):
    r = np.asarray(replicates, dtype=float)
    if r.ndim == 1:
        r = r[:, None]
    d = r-r.mean(axis=0)
    return (len(r)-1.)/len(r)*d.T.dot(d)


def jk_se(replicates):
    r = np.asarray(replicates, dtype=float)
    return np.sqrt((len(r)-1.)/len(r)*np.sum((r-r.mean(axis=0))**2, axis=0))


def energies(mean):
    mean = np.asarray(mean).real
    result = np.full(mean.shape[:-1]+(mean.shape[-1]-1,), np.nan)
    valid = (mean[..., :-1] > 0) & (mean[..., 1:] > 0)
    result[valid] = np.log(mean[..., :-1][valid]/mean[..., 1:][valid])
    return result


def finite_value(value):
    value = float(value)
    return value if np.isfinite(value) else ''


def bool_value(value):
    return bool(value)


def metric_transform(metric):
    # D^{-1} U Lambda^{-1/2} U^dagger; no eigenvalue truncation.
    diagonal = np.diag(metric).real
    if np.any(diagonal <= 0) or not np.isfinite(metric).all():
        raise ValueError('Nonpositive or nonfinite equal-time metric')
    d = np.sqrt(diagonal)
    corr = metric/d[:, None]/d[None, :]
    if np.linalg.norm(corr-corr.conj().T) > 1e-11*np.linalg.norm(corr):
        raise ValueError('Equal-time metric is not Hermitian')
    eigenvalues, u = np.linalg.eigh((corr+corr.conj().T)/2)
    if eigenvalues[0] <= 1e-12*eigenvalues[-1]:
        raise ValueError('Equal-time metric below admitted rank threshold')
    w = ((u/np.sqrt(eigenvalues)).dot(u.conj().T))/d[:, None]
    error = np.linalg.norm(w.conj().T.dot(metric).dot(w)-np.eye(len(d)))
    if error > 1e-9:
        raise ValueError('Metric whitening residual')
    return w, eigenvalues, error


def assert_close(a, b, name, rtol=3e-7, atol=3e-11):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or not np.allclose(a, b, rtol=rtol, atol=atol):
        raise ValueError('Numerical mismatch: '+name)
