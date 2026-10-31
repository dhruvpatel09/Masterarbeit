#!/usr/bin/env python3
"""Optional PyErrors S=0 check of the complete-sample covariance propagation."""
import argparse
import importlib.metadata
import json
import math
import sys
from pathlib import Path

import numpy as np
import pyerrors as pe

import reproduce as analysis

HERE=Path(__file__).resolve().parent


def compare(value,reference):
    if not math.isclose(value,float(reference),rel_tol=2e-8,abs_tol=2e-15):
        raise ValueError('PyErrors and analytic propagation disagree: %r / %s'%(value,reference))


def run(output):
    _,groups,_,_=analysis.load()
    calculated={};report=[];largest=0.
    for key,group in sorted(groups.items()):
        samples=np.array([[float(v) for v in group[c,40]] for c in analysis.CONFIGS])
        primaries=[pe.Obs([samples[:,j]],['qX'],idl=[analysis.CONFIGS]) for j in range(5)]
        L,LS,S,T,LT=primaries
        channels={'RAW':LS/L-S,'TL':(LS-LT)/L-(S-T)}
        row=dict(zip(analysis.KEYS,key))
        for kind,observable in channels.items():
            observable.gamma_method(S=0)
            point,delta=analysis.channel([group[c,40] for c in analysis.CONFIGS],kind)
            reference=analysis.se(delta)
            compare(observable.value,point);compare(observable.dvalue,reference)
            largest=max(largest,abs(observable.dvalue/float(reference)-1))
            row[kind+'_analytic_se']=str(reference)
            row[kind+'_pyerrors_S0_se']=float(observable.dvalue)
            calculated[key,kind]=observable
        report.append(row)
    contrasts=0
    for geometry in sorted({key[1:] for key in groups}):
        for kind in ('RAW','TL'):
            for high,low in ((4,2),(6,4),(6,2)):
                observable=calculated[(high,)+geometry,kind]-calculated[(low,)+geometry,kind]
                observable.gamma_method(S=0)
                _,a=analysis.channel([groups[(high,)+geometry][c,40] for c in analysis.CONFIGS],kind)
                _,b=analysis.channel([groups[(low,)+geometry][c,40] for c in analysis.CONFIGS],kind)
                reference=analysis.se([x-y for x,y in zip(a,b)])
                compare(observable.dvalue,reference)
                largest=max(largest,abs(observable.dvalue/float(reference)-1));contrasts+=1
    versions={name:importlib.metadata.version(name) for name in ('pyerrors','numpy','pandas','scipy','autograd')}
    summary=dict(success=True,configuration_count=71,configuration_labels=analysis.CONFIGS,piter=40,
        profile_geometry_rows=len(report),profile_error_checks=2*len(report),temporal_contrast_checks=contrasts,
        setting='gamma_method(S=0)',ensemble_label='qX',idl='100,110,...,800; the same labels and samples for all primaries and temporal separations',
        maximum_relative_error_difference=largest,versions=versions,python=sys.version.split()[0],
        interpretation='Independent-configuration covariance propagation. This check does not infer a Monte Carlo autocorrelation time or multiply the sample count by spatial sources.')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    analysis.write_csv(output/'pyerrors_S0_checks_71configs.csv',report)
    analysis.write_json(output/'PYERRORS_CHECK.json',summary)
    print('PASS PyErrors S=0: %d profile errors and %d paired temporal errors; all 71 configuration labels retained.'%(2*len(report),contrasts))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=HERE/'results')
    run(p.parse_args().output)
