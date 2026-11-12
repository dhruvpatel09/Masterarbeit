"""Fixed-source tests with training uncertainty and honest two-fold holdout.

The Hermitian part at ONE early training time only defines a real-Rayleigh
source-selection objective. All measured projections use the original complex,
possibly non-Hermitian matrix. This is not a positive-spectrum assumption.
"""
import numpy as np
from statistics_tools import metric_transform

BASES = {'pair': [0, 3], 'triple': [0, 1, 3], 'full': [0, 1, 2, 3]}
ANCHORS = [4, 6, 8]
PROCEDURES = ['same_sample', 'cross_fit']


def submatrix(C, indices):
    return C[..., indices, :][..., :, indices]


def train(mean, anchor):
    metric = mean[0]
    w, spectrum, whitening_error = metric_transform(metric)
    training = (mean[anchor]+mean[anchor].conj().T)/2
    h = w.conj().T.dot(training).dot(w)
    values, vectors = np.linalg.eigh((h+h.conj().T)/2)
    v = w.dot(vectors[:, -1])
    v /= np.sqrt(np.vdot(v, metric.dot(v)).real)
    pivot = np.argmax(np.abs(v))
    v *= np.exp(-1j*np.angle(v[pivot]))
    residual = np.linalg.norm(training.dot(v)-values[-1]*metric.dot(v))
    scale = (np.linalg.norm(training)+abs(values[-1])*np.linalg.norm(metric))*np.linalg.norm(v)
    residual /= scale
    if residual > 1e-10 or abs(np.vdot(v, metric.dot(v))-1) > 1e-9:
        raise ValueError('Source-training eigenvector residual or normalization')
    return v, dict(top_eigenvalue=float(values[-1]), second_eigenvalue=float(values[-2]),
        relative_eigengap=float((values[-1]-values[-2])/max(abs(values[-1]), 1e-300)),
        metric_condition=float(spectrum[-1]/spectrum[0]), rayleigh_residual=float(residual),
        whitening_error=float(whitening_error))


def quadratic(C, v):
    return np.einsum('i,...ij,j->...', v.conj(), C, v, optimize=True)


def estimate(C, identities, anchor, procedure, return_vectors=False):
    records = []
    if procedure == 'same_sample':
        v, info = train(C.mean(axis=0), anchor)
        result = quadratic(C, v).mean(axis=0)
        records.append(('all', len(C), len(C), v, info))
    elif procedure == 'cross_fit':
        total = np.zeros(C.shape[1], dtype=complex)
        for fold in (0, 1):
            test = identities % 2 == fold
            training = ~test
            if not test.any() or not training.any():
                raise ValueError('Empty holdout or training fold')
            v, info = train(C[training].mean(axis=0), anchor)
            total += quadratic(C[test], v).sum(axis=0)
            records.append(('test_fold_'+str(fold), int(training.sum()), int(test.sum()), v, info))
        result = total/len(C)
    else:
        raise ValueError('Unknown source procedure')
    return (result, records) if return_vectors else result


def all_projections(C, radii, progress=False):
    identities = np.arange(len(C))
    means, resamples, keys, vectors = [], [], [], []
    for ri, r in enumerate(radii):
        for basis, indices in BASES.items():
            selected = submatrix(C[:, ri], indices)
            for anchor in ANCHORS:
                for procedure in PROCEDURES:
                    key = dict(R=r, basis=basis, anchor=anchor, procedure=procedure)
                    mean, training = estimate(selected, identities, anchor, procedure, True)
                    reps = []
                    for deleted in identities:
                        keep = identities != deleted
                        reps.append(estimate(selected[keep], identities[keep], anchor, procedure))
                    means.append(mean)
                    resamples.append(reps)
                    keys.append(key)
                    for split, ntrain, ntest, v, info in training:
                        for j, index in enumerate(indices):
                            vectors.append(dict(key, split=split, training_configurations=ntrain,
                                evaluated_configurations=ntest, Nv=[40,60,80,100][index],
                                coefficient_re=float(v[j].real), coefficient_im=float(v[j].imag), **info))
        if progress:
            print('Source projections and full retraining jackknife: R={}'.format(r), flush=True)
    return keys, np.asarray(means), np.asarray(resamples), vectors
