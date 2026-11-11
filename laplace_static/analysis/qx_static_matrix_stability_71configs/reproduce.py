#!/usr/bin/env python3
"""Recompute every split, deletion fit and complex-null bootstrap statistic."""
import argparse,math
import numpy as np
from common import HERE,read,csv_rows,close
from analysis import compute


def scalar(a,b,label):
    if isinstance(b,(bool,np.bool_)):
        if a not in (bool(b),str(bool(b))):raise ValueError('Boolean mismatch '+label)
    elif isinstance(b,(int,np.integer)):
        if str(a)!=str(b):raise ValueError('Count/label mismatch '+label)
    elif isinstance(b,(float,np.floating)):
        if not math.isfinite(float(a)) or not math.isfinite(float(b)):raise ValueError('Nonfinite scalar '+label)
        if label.endswith('stationary_energy_error'):
            if not 0<=float(a)<=1e-10 or not 0<=float(b)<=1e-10:raise ValueError('Fit stationarity bound '+label)
        elif label.endswith('covariance_condition'):
            close([1/float(a)],[1/float(b)],label,rtol=5e-7,atol=8*np.finfo(float).eps)
            if (float(a)>=1e12)!=(float(b)>=1e12):raise ValueError('Condition class '+label)
        else:close([float(a)],[float(b)],label)
    elif str(a)!=str(b):raise ValueError('Text/undefined mismatch '+label)


def json_compare(a,b,label):
    if isinstance(b,dict):
        if set(a)!=set(b):raise ValueError('JSON field inventory '+label)
        for k in b:json_compare(a[k],b[k],label+'.'+k)
    elif isinstance(b,list):
        if len(a)!=len(b):raise ValueError('JSON length '+label)
        for i,(x,y) in enumerate(zip(a,b)):json_compare(x,y,label+'['+str(i)+']')
    else:scalar(a,b,label)


def verify():
    from verify_files import verify as files
    files();tables,arrays,summary=compute(True)
    for name,values in tables.items():
        old=csv_rows(HERE/'results'/(name+'.csv'))
        if len(old)!=len(values):raise ValueError('CSV row count '+name)
        for i,(a,b) in enumerate(zip(old,values)):
            if set(a)!=set(b):raise ValueError('CSV fields '+name)
            for key,value in b.items():scalar(a[key],value,'{} row {} {}'.format(name,i,key))
    old=np.load(str(HERE/'results/stability_resamples.npz'),allow_pickle=False)
    if set(old.files)!=set(arrays):raise ValueError('Resample array inventory')
    for name,value in arrays.items():
        close(old[name].real,value.real,name)
        if np.iscomplexobj(value):close(old[name].imag,value.imag,name+' imaginary')
    json_compare(read(HERE/'results/SUMMARY.json'),summary,'SUMMARY')
    print('PASS complete stability replay: 324 component tests; 1999 bootstrap draws; five holdout splits; 216 curves; 1080 fits and every 71-deletion parameter refit.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',required=True,action='store_true');p.parse_args();verify()
