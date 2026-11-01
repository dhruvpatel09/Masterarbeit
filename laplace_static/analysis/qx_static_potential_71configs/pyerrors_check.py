#!/usr/bin/env python3
"""Optional modern-Python S=0 automatic-differentiation check for all 512 pooled energies."""
import importlib.metadata
import math
import numpy as np
import pyerrors as pe
import reproduce as a


def main():
    raw,_,_,_=a.load();groups={}
    for r in raw:
        if r['scheme']=='fixed' and r['axis']=='-1':groups.setdefault((int(r['tau']),int(r['R'])),{})[int(r['cfg'])]=float(r['L_RE'])
    obs={key:pe.Obs([[g[c] for c in a.CONFIGS]],['qX'],idl=[a.CONFIGS]) for key,g in groups.items()}
    refs={(int(r['tau']),int(r['R'])):r for r in a.read_csv(a.RESULTS/'effective_energies_71configs.csv') if r['axis']=='-1'}
    rows=[];maximum=0.
    for key,ref in sorted(refs.items()):
        t,R=key;value=np.log(obs[t,R]/obs[t+1,R]);value.gamma_method(S=0)
        expected=float(ref['naive_se'])
        if not math.isclose(value.value,float(ref['aE_eff']),rel_tol=2e-10,abs_tol=2e-14) or not math.isclose(value.dvalue,expected,rel_tol=2e-9,abs_tol=2e-15):raise ValueError('PyErrors S=0 disagrees with paired propagation')
        maximum=max(maximum,abs(value.dvalue/expected-1));rows.append(dict(tau=t,R=R,configurations=71,aE_eff=float(value.value),pyerrors_S0_se=float(value.dvalue),analytic_se=ref['naive_se']))
    a.write_csv(a.RESULTS/'pyerrors_S0_check_71configs.csv',rows)
    a.write_json(a.RESULTS/'PYERRORS_CHECK.json',dict(success=True,checks=len(rows),configuration_labels=a.CONFIGS,ensemble='qX',S=0,maximum_relative_error_difference=maximum,versions={k:importlib.metadata.version(k) for k in ('pyerrors','numpy','autograd')}))
    print('PASS PyErrors S=0: all 512 combined-axis effective-energy errors agree with analytic covariance propagation.')


if __name__=='__main__':main()
