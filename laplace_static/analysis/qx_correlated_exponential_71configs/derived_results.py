#!/usr/bin/env python3
"""Paired local contrasts and explicit audits of nonlinear fit errors."""
import csv
import json
import numpy as np
import fit_models as fm
from run_analysis import HERE,ROOT,RESULTS,CUTOFFS,RADII,load_data,parameters,write_csv,write_json


def load_rows():
    rows=[]
    integers=('Nv','R','tau_min','tau_max','states','configurations','dof','nonpositive_means_in_window','successful_starts','total_starts','best_start','jackknife_converged','jackknife_resolved')
    booleans=('identified','jackknife_requested','jackknife_performed','jackknife_reliable')
    strings=('fit_id','flags')
    with (RESULTS/'fits_71configs.csv').open(newline='') as f:
        for raw in csv.DictReader(f):
            row={}
            for k,v in raw.items():
                row[k]=int(v) if k in integers else (v=='True') if k in booleans else v if k in strings else float(v) if v else None
            rows.append(row)
    return rows


def window_for(row,a):
    lo,hi=row['tau_min'],row['tau_max']
    return fm.Window(a[:,CUTOFFS.index(row['Nv']),lo:hi+1,RADII.index(row['R'])],np.arange(lo,hi+1))


def derive(a,rows,solutions):
    cache={}
    for r in rows:
        w=window_for(r,a);p=parameters(r);ei=1 if r['states']==1 else 2
        try:local=w.local_errors(p,r['states'])
        except (ValueError,np.linalg.LinAlgError):local=None
        cache[r['fit_id']]=(w,local,ei)
    lookup={(r['Nv'],r['R'],r['tau_min'],r['tau_max'],r['states']):r for r in rows}
    pairs=[]
    for r in rows:
        if r['states']==1 and r['Nv']!=100:
            other=lookup[100,r['R'],r['tau_min'],r['tau_max'],1]
            kind='cutoff_minus_Nv100'
        elif r['states']==2:
            other=lookup[r['Nv'],r['R'],r['tau_min'],r['tau_max'],1]
            kind='two_minus_one_state'
        else:continue
        _,local,ei=cache[r['fit_id']];_,other_local,ej=cache[other['fit_id']]
        se=None
        if local is not None and other_local is not None:
            delta=local['fluctuations'][:,ei]-other_local['fluctuations'][:,ej]
            se=float(np.std(delta,ddof=1)/np.sqrt(71))
        ja=solutions[r['fit_id']]['jackknife'];jb=solutions[other['fit_id']]['jackknife'];jk=None
        if len(ja)==len(jb)==71 and all(v['success'] for v in ja+jb):
            if [v['omitted_configuration'] for v in ja]!=[v['omitted_configuration'] for v in jb]:raise ValueError('Unmatched jackknife contrast')
            d=np.array([v['p'][ei]-u['p'][ej] for v,u in zip(ja,jb)])
            jk=float(np.sqrt(70/71*np.sum((d-d.mean())**2)))
        pairs.append(dict(contrast=kind,fit_id=r['fit_id'],reference_fit_id=other['fit_id'],Nv=r['Nv'],R=r['R'],
            tau_min=r['tau_min'],tau_max=r['tau_max'],delta_E0=r['E0']-other['E0'],paired_naive_se=se,
            both_local_checks_pass=r['identified'] and other['identified'],paired_jackknife_se_formal=jk,
            paired_jackknife_reliable=r['jackknife_reliable'] and other['jackknife_reliable']))
    audits=[]
    for r in rows:
        if not r['jackknife_requested']:continue
        naive=r['E0_naive_se'];jk=r['E0_jackknife_se_formal']
        audits.append(dict(fit_id=r['fit_id'],Nv=r['Nv'],R=r['R'],tau_min=r['tau_min'],tau_max=r['tau_max'],states=r['states'],
            full_local_checks_pass=r['identified'],requested=True,performed=r['jackknife_performed'],
            converged=r['jackknife_converged'],locally_resolved=r['jackknife_resolved'],
            E0_naive_se=naive,E0_jackknife_se_formal=jk,jackknife_over_naive=jk/naive if jk is not None and naive else None,
            reliable=r['jackknife_reliable'],full_flags=r['flags']))
    source=json.loads((HERE/'inputs/SOURCES.json').read_text())
    path=next(k for k in source['files'] if k.endswith('/constant_fits_71configs.csv'))
    with (ROOT/path).open(newline='') as f:
        old={(int(v['Nv']),int(v['R']),int(v['tau_min'])):v for v in csv.DictReader(f)}
    references=[]
    for r in rows:
        if r['states']!=1:continue
        v=old[r['Nv'],r['R'],r['tau_min']]
        references.append(dict(fit_id=r['fit_id'],Nv=r['Nv'],R=r['R'],correlator_tau_min=r['tau_min'],correlator_tau_max=r['tau_max'],
            correlator_energy=r['E0'],correlator_energy_naive_se=r['E0_naive_se'],
            effective_energy_tau_min=int(v['tau_min']),effective_energy_tau_max=int(v['tau_max']),
            effective_energy_estimate=float(v['estimate']) if v['estimate'] else None,
            effective_energy_naive_se=float(v['naive_se']) if v['naive_se'] else None,
            effective_energy_status=v['status'],effective_energy_low_signal=v['low_signal_flag'],
            note='Different time support; no significance assigned to this difference.'))
    return {'paired_energy_contrasts_71configs.csv':pairs,'parameter_error_audit_71configs.csv':audits,
            'effective_energy_reference_comparison.csv':references}


def main():
    a=load_data();rows=load_rows();solutions=json.loads((RESULTS/'solutions.json').read_text())
    exports=derive(a,rows,solutions)
    for name,records in exports.items():write_csv(RESULTS/name,records)
    summary=json.loads((RESULTS/'SUMMARY.json').read_text())
    ratios=[r['jackknife_over_naive'] for r in exports['parameter_error_audit_71configs.csv'] if r['reliable']]
    summary['derived_exports']={name:len(records) for name,records in exports.items()}
    summary['reliable_audit_jackknife_over_naive_range']=[min(ratios),max(ratios)]
    summary['identifiability_note']='The identified field reports local parameter diagnostics, not goodness of fit, a plateau, or physical ground-state identification.'
    write_json(RESULTS/'SUMMARY.json',summary)
    print('PASS paired exponential contrasts and explicit nonlinear-error audits')


if __name__=='__main__':main()
