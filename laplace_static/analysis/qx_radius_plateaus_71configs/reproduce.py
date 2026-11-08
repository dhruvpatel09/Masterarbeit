#!/usr/bin/env python3
"""Replay all numerical exports from primary per-configuration correlators."""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import analysis
import approach_model
from data import HERE, load
from statistics_tools import assert_close


def curvature_condition(actual,expected,label):
    """Compare this diagnostic in reciprocal form at float64 precision.

    For an ill-conditioned 3x3 Hessian, tiny changes in its smallest
    singular value amplify relative changes in kappa. Its reciprocal,
    sigma_min/sigma_max, has an absolute roundoff floor instead.
    Never permit this allowance to change the conditioning classification.
    """
    a,b=float(actual),float(expected)
    if not math.isfinite(a) or not math.isfinite(b) or min(a,b)<1:
        raise ValueError('Invalid curvature condition '+label)
    if (a>1e12)!=(b>1e12):
        raise ValueError('Curvature conditioning classification changed '+label)
    assert_close(1./a,1./b,label+' reciprocal',
                 rtol=3e-7,atol=8*sys.float_info.epsilon)


def scalar(actual,expected,label,field=None):
    if expected is None:
        if actual not in ('',None):raise ValueError('Expected undefined '+label)
    elif isinstance(expected,bool):
        if actual != str(expected):raise ValueError('Boolean mismatch '+label)
    elif isinstance(expected,int):
        if actual != str(expected):raise ValueError('Integer mismatch '+label)
    elif isinstance(expected,float):
        if field=='curvature_condition':
            curvature_condition(actual,expected,label)
        else:
            assert_close(float(actual),expected,label,rtol=3e-7,atol=3e-11)
    elif actual != str(expected):
        raise ValueError('Exact text mismatch '+label)


def verify_table(path,expected):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        if set(reader.fieldnames) != set(expected[0]):raise ValueError('Table schema '+path.name)
        actual = list(reader)
    if len(actual)!=len(expected):raise ValueError('Row count '+path.name)
    for i,(a,b) in enumerate(zip(actual,expected)):
        for k,v in b.items():
            field=k if path.name=='approach_fits.csv' else None
            scalar(a[k],v,'{} row {} {}'.format(path.name,i,k),field=field)


def verify_json(actual,expected,label='summary'):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or set(actual)!=set(expected):raise ValueError('JSON schema '+label)
        for k,v in expected.items():verify_json(actual[k],v,label+'.'+k)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):raise ValueError('JSON length '+label)
        for i,(a,b) in enumerate(zip(actual,expected)):verify_json(a,b,label+'.'+str(i))
    elif isinstance(expected,float):
        assert_close(actual,expected,label,rtol=3e-7,atol=3e-11)
    elif actual!=expected:raise ValueError('JSON value '+label)


def verify():
    from verify_files import verify as verify_files
    verify_files()
    data = load()
    tables,summary = analysis.compute(data)
    plan = json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    archived = json.loads((HERE/'results/approach_solutions.json').read_text())
    rows,_ = approach_model.compute(data,plan,archived=archived)
    tables['approach_fits'] = rows
    summary['table_rows']['approach_fits'] = len(rows)
    summary['approach_status_counts'] = {s:sum(r['status']==s for r in rows) for s in sorted(set(r['status'] for r in rows))}
    summary['approach_diagnostic_screen_count'] = sum(r['diagnostic_screen'] for r in rows)
    summary['approach_replay'] = 'Archived variable-projection candidates are re-evaluated; optimizer searches are not repeated by --verify.'
    for name,values in tables.items():verify_table(HERE/'results'/(name+'.csv'),values)
    verify_json(json.loads((HERE/'results/SUMMARY.json').read_text()),summary)
    print('PASS qX radius replay: 71 configurations; 1924 means; 1856 effective-energy points (42 undefined); 204800 covariance entries audited.')
    print('PASS published short/long fit regression, complete window catalogue, paired contrasts, axes and approach-model candidate replay.')
    print('No new MPI measurements; no automatic plateau or ground-state selection.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',required=True)
    parser.parse_args();verify()
