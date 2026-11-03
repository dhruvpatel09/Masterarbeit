#!/usr/bin/env python3
"""Recompute every archived numerical result with explicit portable tolerances."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from analysis import compute
from replay import compare,replay

HERE=Path(__file__).resolve().parent


def read_mode_labels(path,plan):
    """Retain the published mode convention without supplying fitted values."""
    labels={};inventory=set();expected=set()
    for radius in plan['radii']:
        for basis,cutoffs in plan['bases'].items():
            n=len(cutoffs)
            for tref in plan['reference_times']:
                samples=['mean']+(['drop%d'%cfg for cfg in plan['configuration_labels']] if tref==0 else [])
                for sample in samples:
                    key=(sample,radius,basis,tref);labels[key]=np.full((33,n),complex(np.nan,np.nan))
                    expected.update(key+(tau,mode) for tau in range(tref+1,33) for mode in range(n))
    with Path(path).open(newline='') as f:
        for row in csv.DictReader(f):
            key=(row['sample'],int(row['R']),row['basis'],int(row['t_ref']));tau,mode=int(row['tau']),int(row['mode']);entry=key+(tau,mode)
            if entry not in expected or entry in inventory:raise ValueError('Invalid or duplicate archived mode label')
            value=complex(float(row['lambda_re']),float(row['lambda_im']))
            if not np.isfinite(value):raise ValueError('Nonfinite archived mode label')
            labels[key][tau,mode]=value;inventory.add(entry)
    if inventory!=expected:raise ValueError('Incomplete archived mode-label inventory')
    return labels


def numerical_difference(first,second,angular=False):
    if not math.isfinite(first) or not math.isfinite(second):return math.inf
    if angular:
        if max(abs(first),abs(second))>math.pi+2e-10:return math.inf
        return abs(math.atan2(math.sin(first-second),math.cos(first-second)))
    return abs(first-second)


def verify_table(path,expected):
    with path.open(newline='') as f:recorded=list(csv.DictReader(f))
    if len(recorded)!=len(expected):raise ValueError('Table row inventory changed: '+path.name)
    for i,(actual,wanted) in enumerate(zip(recorded,expected)):
        if set(actual)!=set(wanted):raise ValueError('Table fields changed: '+path.name)
        for key,b in wanted.items():
            a=actual[key]
            if isinstance(b,float):
                a=float(a)
                if not math.isfinite(a) or not math.isfinite(b):raise ValueError('Nonfinite numerical replay value')
                if key in ('relative_residual','maximum_relative_residual'):
                    if not math.isfinite(a) or not math.isfinite(b) or not 0<=a<=1e-10 or not 0<=b<=1e-10:raise ValueError('Eigenpair residual failed')
                elif numerical_difference(a,b,angular=key in ('phase','effective_energy_im'))>2e-10+2e-8*max(abs(a),abs(b)):
                    raise ValueError('Numerical mismatch %s row %d field %s: archived %.17g; recomputed %.17g'%(path.name,i,key,a,b))
            elif a!=str(b):raise ValueError('Inventory mismatch %s row %d field %s'%(path.name,i,key))


def verify_exports(tables,summary,validation):
    compare(json.loads((HERE/'results/VALIDATION.json').read_text()),validation,'validation')
    for name,values in tables.items():verify_table(HERE/'results'/(name+'.csv'),values)
    recorded=json.loads((HERE/'results/SUMMARY.json').read_text());summary=dict(summary)
    for obj in (recorded,summary):
        residual=obj.pop('maximum_generalized_eigenpair_relative_residual')
        if not 0<=residual<=1e-10:raise ValueError('Generalized eigenpair residual failed')
    compare(recorded,summary,'summary',rtol=2e-8,atol=2e-10)


def verify():
    from verify_files import verify as verify_files
    verify_files()
    plan=json.loads((HERE/'ANALYSIS_PLAN.json').read_text());labels=read_mode_labels(HERE/'results/complex_modes.csv',plan)
    rows,_,validation=replay();tables,summary=compute(rows,validation,mode_labels=labels)
    verify_exports(tables,summary,validation)
    print('PASS matrix pilot replay: 3 audits; 4 scans; 6336 exact p10 diagonal matches; all complex matrix, conditioning, mode-tracking and fixed-projection exports reproduced.')
    print('Fresh eigenpairs matched to archived mode labels; principal phases compared modulo 2*pi; all spectra and residual bounds verified.')
    print('Three-configuration sensitivity diagnostics only; no new MPI measurements or production submissions.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',required=True,action='store_true');parser.parse_args();verify()
