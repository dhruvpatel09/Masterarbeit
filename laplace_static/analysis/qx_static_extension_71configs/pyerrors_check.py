#!/usr/bin/env python3
"""Optional modern-Python PyErrors S=0 check, preserving matched samples."""
import importlib.metadata
import math
import numpy as np
import pyerrors as pe

import reproduce as a
from verify_files import verify


def main():
    verify()
    groups={}
    with (a.SOURCE/'campaigns/late_modes_full71/correlators.csv').open(newline='') as stream:
        import csv
        for row in csv.DictReader(stream):
            if row['scheme']=='fixed' and row['axis']=='-1':
                key=tuple(int(row[k]) for k in ('Nv','tau','R'))
                groups.setdefault(key,{})[int(row['cfg'])]=float(row['L_RE'])
    obs={key:pe.Obs([[g[c] for c in a.CONFIGS]],['qX'],idl=[a.CONFIGS]) for key,g in groups.items()}
    derived={};rows=[];maximum=0.;invalid=0
    for ref in a.read_csv(a.RESULTS/'effective_energies_71configs.csv'):
        n,t,R=(int(ref[k]) for k in ('Nv','tau','R'))
        if ref['valid_point']!='True':invalid+=1;continue
        value=np.log(obs[n,t,R]/obs[n,t+1,R]);value.gamma_method(S=0)
        expected=float(ref['naive_se'])
        if not math.isclose(value.value,float(ref['aE_eff']),rel_tol=2e-10,abs_tol=2e-13) or not math.isclose(value.dvalue,expected,rel_tol=2e-9,abs_tol=2e-15):
            raise ValueError('PyErrors S=0 disagrees with matched naive propagation')
        maximum=max(maximum,abs(value.dvalue/expected-1));derived[n,t,R]=value
        rows.append(dict(Nv=n,tau=t,R=R,configurations=71,aE_eff=float(value.value),
            pyerrors_S0_se=float(value.dvalue),analytic_naive_se=ref['naive_se']))
    comparisons=[];maximum_contrast=0.
    for ref in a.read_csv(a.RESULTS/'cutoff_energy_contrasts_71configs.csv'):
        n,t,R=(int(ref[k]) for k in ('Nv','tau','R'))
        if ref['diagnostic']!='defined':continue
        value=derived[n,t,R]-derived[100,t,R];value.gamma_method(S=0)
        expected=float(ref['paired_naive_se'])
        if not math.isclose(value.value,float(ref['difference']),rel_tol=2e-9,abs_tol=2e-13) or not math.isclose(value.dvalue,expected,rel_tol=2e-9,abs_tol=2e-15):
            raise ValueError('PyErrors S=0 disagrees with the paired cutoff contrast')
        maximum_contrast=max(maximum_contrast,abs(value.dvalue/expected-1))
        comparisons.append(dict(Nv=n,reference_Nv=100,tau=t,R=R,configurations=71,
            difference=float(value.value),pyerrors_S0_se=float(value.dvalue),analytic_paired_se=ref['paired_naive_se']))
    a.write_csv(a.RESULTS/'pyerrors_S0_energies_71configs.csv',rows)
    a.write_csv(a.RESULTS/'pyerrors_S0_cutoff_contrasts_71configs.csv',comparisons)
    a.write_json(a.RESULTS/'PYERRORS_CHECK.json',dict(success=True,S=0,configuration_labels=a.CONFIGS,
        defined_energy_checks=len(rows),undefined_energy_points_retained=invalid,paired_cutoff_checks=len(comparisons),
        maximum_relative_energy_error_difference=maximum,maximum_relative_contrast_error_difference=maximum_contrast,
        versions={k:importlib.metadata.version(k) for k in ('pyerrors','numpy','autograd')},
        interpretation='S=0 validates the independent-configuration covariance algebra; it does not establish Gaussian errors for weak signals or a ground-state plateau.'))
    print('PASS PyErrors S=0: %d defined energies and %d paired cutoff contrasts.'%(len(rows),len(comparisons)))


if __name__=='__main__':main()
