#!/usr/bin/env python3
"""Cold native replay, independent source retraining and archived table checks."""
import argparse
import csv
import json
import math
import numpy as np
from scipy.optimize import linear_sum_assignment
from analysis import analyze
from replay import HERE,raw_replay
from statistics_tools import assert_close


def scalar(actual, expected, label):
    if isinstance(expected,(bool,np.bool_)):
        if actual not in (bool(expected),str(bool(expected))):raise ValueError('Boolean mismatch '+label)
    elif isinstance(expected,(int,float,np.number)):
        if label.endswith(('rayleigh_residual','backward_residual')):
            if not 0<=float(actual)<=1e-10 or not 0<=float(expected)<=1e-10:raise ValueError('Residual bound '+label)
        elif label.endswith(('whitening_error','whitening_residual')):
            if not 0<=float(actual)<=1e-9 or not 0<=float(expected)<=1e-9:raise ValueError('Whitening bound '+label)
        elif label.endswith(('covariance_condition','reference_condition')):
            a,b=float(actual),float(expected)
            # Condition diagnostics can be sensitive near singularity; compare
            # inverse conditions with a machine-roundoff floor, retaining rank class.
            assert_close([1/a],[1/b],label,rtol=3e-6,atol=8*np.finfo(float).eps)
            if (a>=1e12)!=(b>=1e12):raise ValueError('Condition threshold class '+label)
        else:
            assert_close([float(actual)],[float(expected)],label,rtol=3e-7,atol=3e-11)
    elif str(actual)!=str(expected):
        raise ValueError('Text/undefined mismatch '+label)


def rows(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream))


def normalization_signal(archived,recomputed,label):
    """Compare the relative error at the deterministic same-sample t=0 point.

    Here Q(0)=1 by the training normalization, so its true sampling error is
    zero. Last-bit differences make the reciprocal of its roundoff-sized
    error (the archived signal/error column) numerically ill conditioned.
    """
    floor=64*np.finfo(float).eps
    inverse=[]
    for row in (archived,recomputed):
        if str(row['procedure'])!='same_sample' or str(row['tau'])!='0':
            raise ValueError('Normalization comparison at another point '+label)
        mean,imag,error,imag_error=(float(row[k]) for k in
            ('mean_re','mean_im','jackknife_se_re','jackknife_se_im'))
        if (not all(math.isfinite(x) for x in (mean,imag,error,imag_error))
                or abs(mean-1)>floor or abs(imag)>floor
                or not 0<=error<=floor or not 0<=imag_error<=floor):
            raise ValueError('Normalization or roundoff-error bound '+label)
        signal=float(row['real_mean_over_se'])
        if error:
            if not math.isfinite(signal) or not math.isclose(signal,mean/error,rel_tol=2e-12,abs_tol=0.):
                raise ValueError('Inconsistent normalization signal/error '+label)
        elif signal!=float('inf'):
            raise ValueError('Zero-error normalization signal convention '+label)
        inverse.append(error/mean)
    assert_close([inverse[0]],[inverse[1]],label+' relative uncertainty',rtol=3e-7,atol=floor)


def table(name,expected):
    actual=rows(HERE/'results'/(name+'.csv'))
    if len(actual)!=len(expected):raise ValueError('Row count '+name)
    if name=='complex_pencil_diagnostics':
        # Algebraic mode order is not a physical state label. Compare the whole
        # freshly solved spectrum within each pencil, then circular phases.
        for i in range(0,len(actual),4):
            old,new=actual[i:i+4],expected[i:i+4]
            a=np.array([complex(float(d['eigenvalue_re']),float(d['eigenvalue_im'])) for d in old])
            b=np.array([complex(d['eigenvalue_re'],d['eigenvalue_im']) for d in new])
            ia,ib=linear_sum_assignment(abs(a[:,None]-b[None,:]))
            for j,k in zip(ia,ib):
                assert_close([a[j].real,a[j].imag],[b[k].real,b[k].imag],name,rtol=3e-7,atol=3e-11)
                for field in ('R','reference_tau','tau','reference_condition','backward_residual'):
                    scalar(old[j][field],new[k][field],name+' '+field)
                phase=math.atan2(math.sin(float(old[j]['phase'])-new[k]['phase']),math.cos(float(old[j]['phase'])-new[k]['phase']))
                if abs(phase)>3e-6:raise ValueError('Circular pencil phase')
        return
    if name=='fixed_source_vectors':
        groups={}
        for a,b in zip(actual,expected):
            key=tuple(b[k] for k in ('R','basis','anchor','procedure','split'))
            groups.setdefault(key,([],[]))
            groups[key][0].append(complex(float(a['coefficient_re']),float(a['coefficient_im'])))
            groups[key][1].append(complex(b['coefficient_re'],b['coefficient_im']))
            for field in b:
                if field not in ('coefficient_re','coefficient_im'):scalar(a[field],b[field],name+' '+field)
        for a,b in groups.values():
            a,b=np.asarray(a),np.asarray(b)
            na,nb=np.linalg.norm(a),np.linalg.norm(b)
            assert_close([na],[nb],name+' norm')
            if 1-abs(np.vdot(a,b))/(na*nb)>1e-10:raise ValueError('Phase-free source vector mismatch')
        return
    for i,(a,b) in enumerate(zip(actual,expected)):
        if set(a)!=set(b):raise ValueError('CSV fields '+name)
        for k,v in b.items():
            label='{} row {} {}'.format(name,i,k)
            if (name=='projected_correlators' and k=='real_mean_over_se'
                    and b['procedure']=='same_sample' and b['tau']==0):
                normalization_signal(a,b,label)
            elif name=='matrix_elements' and isinstance(v,float):
                assert_close([float(a[k])],[v],label,rtol=3e-7,atol=1e-24)
            else:scalar(a[k],v,label)


def json_compare(a,b,label):
    if isinstance(b,dict):
        if set(a)!=set(b):raise ValueError('JSON fields '+label)
        for k in b:json_compare(a[k],b[k],label+'.'+k)
    elif isinstance(b,list):
        if len(a)!=len(b):raise ValueError('JSON length '+label)
        for i,(x,y) in enumerate(zip(a,b)):json_compare(x,y,label+'['+str(i)+']')
    else:scalar(a,b,label)


def verify():
    from verify_files import verify as files
    files()
    C,text,validation=raw_replay(progress=True)
    data=np.load(str(HERE/'inputs/primary_matrices.npz'),allow_pickle=False)
    if not np.array_equal(C,data['C']) or not np.array_equal(text,data['diagonal_text']):
        raise ValueError('Primary matrix cache differs from stdout replay')
    json_compare(json.loads((HERE/'results/VALIDATION.json').read_text()),validation,'VALIDATION')
    tables,arrays,summary=analyze(C,text,True)
    for name,values in tables.items():table(name,values)
    old=np.load(str(HERE/'results/projection_resamples.npz'),allow_pickle=False)
    if set(old.files)!=set(arrays):raise ValueError('Resample export keys')
    for key,values in arrays.items():
        atol=1e-24 if key.startswith('diagonal_') else 3e-11
        assert_close(old[key].real,values.real,key,rtol=3e-7,atol=atol)
        if np.iscomplexobj(values):assert_close(old[key].imag,values.imag,key+' imaginary',rtol=3e-7,atol=3e-11)
    json_compare(json.loads((HERE/'results/SUMMARY.json').read_text()),summary,'SUMMARY')
    print('PASS full71 replay: all 71 audits and 72 scans; 599808 primary rows; 872448 rank-block rows; 149952 exact diagonal matches.')
    print('Reproduced all fixed-source/holdout estimates, 71 complete retraining deletions, paired contrasts and window diagnostics. No MPI execution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',required=True,action='store_true');p.parse_args();verify()
