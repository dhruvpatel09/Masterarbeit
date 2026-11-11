"""Pinned data access and configuration-level numerical utilities."""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PREVIOUS=ROOT/'laplace_static/analysis/qx_static_matrix_71configs'
RADII=[4,12,24,32]
CUTOFFS=[40,60,80,100]
CONFIGS=list(range(100,801,10))


def read(path):return json.loads(Path(path).read_text())


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(1048576),b''):h.update(b)
    return h.hexdigest()


def verify_dependencies():
    source=read(HERE/'inputs/SOURCES.json')
    for name,r in source['dependencies'].items():
        p=ROOT/name
        if p.stat().st_size!=r['bytes'] or digest(p)!=r['sha256']:
            raise ValueError('Published dependency changed: '+name)
    b=(PREVIOUS/'FILE_MANIFEST.json').read_bytes()
    if hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()!='5a99731faf0e0e097a670865dfa4d6653337266b':
        raise ValueError('Previous package manifest differs from published 2b9830b')


def load_data():
    verify_dependencies()
    a=np.load(str(PREVIOUS/'inputs/primary_matrices.npz'),allow_pickle=False)
    if (a['C'].shape!=(71,4,33,4,4) or not np.isfinite(a['C']).all()
            or list(a['configurations'])!=CONFIGS or list(a['radii'])!=RADII
            or list(a['cutoffs'])!=CUTOFFS):raise ValueError('Primary matrix cohort or orientation')
    return a['C']


def previous_projections():
    """Load the pinned source algorithm without polluting generic module names."""
    old=sys.modules.get('statistics_tools');previous=sys.path[:]
    sys.modules.pop('statistics_tools',None);sys.path.insert(0,str(PREVIOUS))
    try:
        spec=importlib.util.spec_from_file_location('qx_published_source_projection',str(PREVIOUS/'projections.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module
    finally:
        sys.path[:]=previous;sys.modules.pop('statistics_tools',None)
        if old is not None:sys.modules['statistics_tools']=old


def csv_rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def number(value):
    value=float(value)
    return value if np.isfinite(value) else ''


def jk_cov(samples):
    x=np.asarray(samples,float)
    if x.ndim==1:x=x[:,None]
    d=x-x.mean(axis=0)
    return (len(x)-1.)/len(x)*d.T.dot(d)


def jk_se(samples):
    x=np.asarray(samples,float)
    return np.sqrt((len(x)-1.)/len(x)*np.sum((x-x.mean(axis=0))**2,axis=0))


def check_equal_time_reality(mean,deletions,label):
    """A Hermitian equal-time quadratic form is real, within bounded roundoff.

    Keep the raw imaginary numbers, but never studentize their roundoff.
    Check central and deletion projections so a substantive defect cannot
    pass merely because the corresponding significance ratio is undefined.
    """
    values=np.concatenate([np.atleast_1d(mean),np.asarray(deletions).ravel()])
    if not np.isfinite(values).all():raise ValueError('Nonfinite equal-time projection '+label)
    scale=float(np.max(abs(values.real)))
    if scale<=0:raise ValueError('Nonpositive equal-time projection scale '+label)
    bound=64*np.finfo(float).eps*scale
    maximum=float(np.max(abs(values.imag)))
    if maximum>bound:
        raise ValueError('Equal-time imaginary projection exceeds roundoff bound: {} ({:.6g} > {:.6g})'.format(label,maximum,bound))
    return bound


def energies(mean):
    a=np.asarray(mean).real;out=np.full(a.shape[:-1]+(a.shape[-1]-1,),np.nan)
    valid=(a[...,:-1]>0)&(a[...,1:]>0)
    out[valid]=np.log(a[...,:-1][valid]/a[...,1:][valid])
    return out


def close(a,b,label,rtol=5e-7,atol=5e-11):
    a,b=np.asarray(a),np.asarray(b)
    if (a.shape!=b.shape or np.isinf(a).any() or np.isinf(b).any()
            or not np.array_equal(np.isnan(a),np.isnan(b))):
        raise ValueError('Numerical geometry/validity mismatch '+label)
    good=np.isfinite(a)&np.isfinite(b)
    actual,replayed=a[good],b[good]
    allowed=atol+rtol*abs(replayed)
    failed=abs(actual-replayed)>allowed
    if np.any(failed):
        i=int(np.flatnonzero(failed)[0])
        raise ValueError('Numerical mismatch {}: archived={!r}, replayed={!r}, absolute difference={:.6g}, tolerance={:.6g}'.format(
            label,float(actual[i]),float(replayed[i]),float(abs(actual[i]-replayed[i])),float(allowed[i])))
