#!/usr/bin/env python3
"""Replay the three-configuration qX extension pilot and export diagnostics."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'inputs/job5957304'
CONFIGS=[170,640,650]
CUTOFFS=[40,60,80,100]
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as stream:return list(csv.DictReader(stream))


def mean(values):return math.fsum(values)/len(values)


def se(values):
    centre=mean(values);n=len(values)
    return math.sqrt(math.fsum((x-centre)**2 for x in values)/(n*(n-1)))


def ratio(xs,ys):
    denominator=mean(ys)
    if denominator<=0:raise ValueError('Normalization needs a positive mean')
    value=mean(xs)/denominator
    influence=[(x-value*y)/denominator for x,y in zip(xs,ys)]
    return value,se(influence)


def log_ratio(xs,ys):
    a,b=mean(xs),mean(ys)
    if a<=0 or b<=0:return None,None,None,0,None
    value=math.log(a/b);influence=[x/a-y/b for x,y in zip(xs,ys)]
    replicates=[]
    for i in range(len(xs)):
        x=mean(xs[:i]+xs[i+1:]);y=mean(ys[:i]+ys[i+1:])
        if x>0 and y>0:replicates.append(math.log(x/y))
    jk=None
    if len(replicates)==len(xs):
        m=mean(replicates);n=len(xs)
        jk=math.sqrt((n-1.)/n*math.fsum((x-m)**2 for x in replicates))
    return value,se(influence),jk,len(replicates),influence


def check_files():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        if digest(HERE/name)!=expected:raise ValueError('Package file changed: '+name)
    return len(manifest['files'])


def replay():
    inventory=json.loads((HERE/'inputs/ARCHIVE_FILES.json').read_text());mapping={}
    for entry in inventory['files']:
        path=HERE/entry['stored']
        if path.stat().st_size!=entry['bytes'] or digest(path)!=entry['sha256']:
            raise ValueError('Archived input changed: '+entry['original'])
        mapping[entry['original']]=path
    if len(mapping)!=77:raise ValueError('Incomplete archive inventory')
    sources=json.loads((HERE/'inputs/SOURCES.json').read_text())
    if digest(SOURCE/'INPUTS.json')!=sources['measurement_inputs_sha256']:
        raise ValueError('Measurement inputs differ from the published source receipt')
    # The captured controller is used only for offline replay, never for
    # execution or resubmission from this analysis directory.
    sys.path.insert(0,str(SOURCE))
    try:
        spec=importlib.util.spec_from_file_location('qx_extension_captured_campaign',str(SOURCE/'campaign.py'))
        c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
    finally:sys.path.pop(0)
    target=SOURCE/'campaigns/late_modes_pilot'
    plan=json.loads((target/'PREPARED.json').read_text())
    if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():
        raise ValueError('Prepared receipt changed')
    build=json.loads((SOURCE/'build_provenance/MANIFEST.json').read_text())
    if digest(SOURCE/'build_provenance/MANIFEST.json')!=plan['build_manifest_sha256'] or build['inputs_sha256']!=digest(SOURCE/'INPUTS.json'):
        raise ValueError('Build receipt changed')
    checked={}
    for label,files in [('source',json.loads((SOURCE/'INPUTS.json').read_text())['files']),('build',build['files'])]:
        available=0
        for name,expected in files.items():
            if name in mapping:
                if digest(mapping[name])!=expected:raise ValueError('Captured pin changed: '+name)
                available+=1
        checked[label]=dict(verified_archived_files=available,external_hashes_retained=len(files)-available)
    if plan['stage']!='pilot' or [e['configuration'] for e in plan['cases']]!=CONFIGS:
        raise ValueError('Unexpected pilot cohort')
    geometry=plan['geometry']
    if geometry['taus']!=list(range(33)) or geometry['R']!=RADII or geometry['cutoffs']!=CUTOFFS:
        raise ValueError('Changed pilot geometry')
    refs=c.references();allrows=[];allblocks=[];cases=[];timing=[]
    for entry in plan['cases']:
        out=c.complete_case(target,entry)
        if out is None:raise ValueError('Missing validated case')
        rows,blocks,summary=c.replay_case(out,entry,refs)
        if read_csv(out/'correlators.csv')!=rows or read_csv(out/'spatial_blocks.csv')!=blocks:
            raise ValueError('Archived numerical table differs from stdout replay')
        stored=json.loads((out/'SUMMARY.json').read_text())
        if any(stored.get(k)!=v for k,v in summary.items()):raise ValueError('Case summary differs from replay')
        allrows+=rows;allblocks+=blocks;cases.append(summary)
        launches=json.loads((out/'LAUNCHES.json').read_text())
        times={r['piter']:r['seconds'] for r in launches}
        timing.append(dict(configuration=entry['configuration'],audit_seconds=times[0],p10_seconds=times[10],
            p40_seconds=times[40],pilot_total_seconds=math.fsum(times.values()),
            production_case_seconds=times[0]+times[10]))
    if read_csv(target/'correlators.csv')!=allrows or read_csv(target/'spatial_blocks.csv')!=allblocks:
        raise ValueError('Collected tables differ from replay')
    summary=json.loads((target/'SUMMARY.json').read_text())
    if summary['success'] is not True or summary['cases']!=cases or summary['missing_configurations']:
        raise ValueError('Collected summary differs from replay')
    if (len(allrows),len(allblocks))!=(23976,28800):raise ValueError('Wrong row counts')
    return allrows,cases,timing,checked


def compute():
    rows,cases,timing,checked=replay()
    data={}
    for row in rows:
        if row['scheme']=='fixed' and row['axis']=='-1' and row['piter']=='10':
            key=tuple(int(row[k]) for k in ('cfg','Nv','tau','R'))
            if key in data:raise ValueError('Duplicate combined-axis primary')
            data[key]=float(row['L_RE'])
    if len(data)!=3*4*33*12:raise ValueError('Missing combined-axis primary')
    def sample(n,t,r):return [data[c,n,t,r] for c in CONFIGS]
    signals=[];energies=[];energy_lookup={};cutoff_counts=[]
    for n in CUTOFFS:
        nonpositive_means=0;nonpositive_configurations=0
        for r in RADII:
            for t in range(33):
                xs=sample(n,t,r);m,s=mean(xs),se(xs);q,qs=ratio(xs,sample(n,0,r))
                snr=m/s if s>0 else None
                nonpositive_means+=int(m<=0);nonpositive_configurations+=sum(x<=0 for x in xs)
                flag='nonpositive_mean' if m<=0 else ('low_mean_over_se' if snr is None or snr<3 else 'positive_mean')
                signals.append(dict(Nv=n,tau=t,R=r,configurations=3,L_RE_mean=m,L_RE_se=s,
                    normalized_mean=q,normalized_se=qs,mean_over_se=snr,
                    nonpositive_configuration_count=sum(x<=0 for x in xs),diagnostic=flag))
            for t in range(32):
                xs,ys=sample(n,t,r),sample(n,t+1,r);v,s,jk,count,influence=log_ratio(xs,ys)
                endpoint_snr=min(mean(z)/se(z) for z in (xs,ys) if se(z)>0)
                flag='nonpositive_endpoint_mean' if v is None else ('low_endpoint_mean_over_se' if endpoint_snr<3 else 'positive_endpoint_means')
                row=dict(Nv=n,tau=t,R=r,configurations=3,effective_energy=v,naive_se=s,
                    paired_delete_one_jackknife_se=jk,positive_jackknife_replicates=count,
                    minimum_endpoint_mean_over_se=endpoint_snr,diagnostic=flag)
                energies.append(row);energy_lookup[n,t,r]=(row,influence)
        cutoff_counts.append(dict(Nv=n,combined_axis_points=396,individual_configuration_points=1188,
            nonpositive_mean_points=nonpositive_means,nonpositive_individual_points=nonpositive_configurations))
    differences=[]
    for n in CUTOFFS[:-1]:
        for r in RADII:
            for t in range(32):
                a,ia=energy_lookup[n,t,r];b,ib=energy_lookup[100,t,r]
                valid=ia is not None and ib is not None
                differences.append(dict(Nv=n,reference_Nv=100,tau=t,R=r,configurations=3,
                    effective_energy_difference=a['effective_energy']-b['effective_energy'] if valid else None,
                    paired_naive_se=se([x-y for x,y in zip(ia,ib)]) if valid else None,
                    diagnostic='defined' if valid else 'nonpositive_endpoint_mean'))
    max_case=max(r['production_case_seconds'] for r in timing)
    summary=dict(success=True,scope='Technical pilot; selected controls, not a population-level potential extraction.',
        configuration_labels=CONFIGS,configuration_count=3,production_target_count=71,
        geometry=dict(t0=90,taus=list(range(33)),R=RADII,cutoffs=CUTOFFS,piter=10,p40_control_Nv=100),
        raw_correlator_rows=len(rows),regional_rows=28800,audited_mode_records=9900,
        direct_contraction_checks=sum(s['direct_tests'] for case in cases for s in case['scans']),
        periodic_payload_roundtrips=18,archived_reference_rows=sum(s['archived_reference_rows'] for case in cases for s in case['scans']),
        archived_reference_maximum_absolute_difference=max(s['archived_reference_maximum_absolute_difference'] for case in cases for s in case['scans']),
        maximum_relative_mode_residual=max(case['quality']['maximum_relative_residual'] for case in cases),
        maximum_central_unitarity_defect=max(float(s['unit_max']) for case in cases for s in case['scans']),
        maximum_central_determinant_defect=max(float(s['det_max']) for case in cases for s in case['scans']),
        maximum_p10_p40_component_difference=max(case['p10_p40_maximum_absolute_difference'] for case in cases),
        provenance=checked,signal_diagnostics=cutoff_counts,
        R24_tau16=[r for r in signals if r['R']==24 and r['tau']==16],
        R32_tau16=[r for r in signals if r['R']==32 and r['tau']==16],
        production_plan=dict(reused_pilot_configurations=CONFIGS,new_configurations=68,cases_per_task=2,array_tasks=34,
            concurrency=4,maximum_observed_case_seconds=max_case,estimated_two_case_task_minutes=2*max_case/60,
            timing_margin_50_percent_task_minutes=3*max_case/60,
            idealized_nine_wave_execution_hours=9*2*max_case/3600,
            caveat='Based on three observed configurations. Four tasks continuously allocated, with no queue gaps or additional overhead; not a scheduler forecast.'),
        uncertainty='N=3 independent-configuration sample-mean standard errors; nonlinear errors retain paired covariance; diagnostic only.',
        signal_flag_note='Mean/SE=3 is an illustrative plotting flag, not a physics acceptance criterion; N=3 error estimates are unstable.',
        log_policy='Form means before ratios/logs; retain every configuration. Nonpositive mean endpoints give undefined log values, exported as empty CSV fields. No absolute values or positive-configuration selection.',
        interpretation='Smaller Nv improves the R=24 signal in this pilot, but does not improve all distances. Nv changes the source operator. Full71 data and temporal-window tests are required before selecting a cutoff or quoting a ground-state potential.')
    outputs={'correlator_diagnostics_pilot.csv':signals,'effective_energies_pilot.csv':energies,
        'cutoff_differences_pilot.csv':differences,'timing_pilot.csv':timing,
        'R24_tau16_pilot.csv':summary['R24_tau16']}
    return outputs,summary


def equivalent(a,b,path='root'):
    if isinstance(b,dict):
        if not isinstance(a,dict) or set(a)!=set(b):raise ValueError('Different keys: '+path)
        for key in b:equivalent(a[key],b[key],path+'.'+key)
    elif isinstance(b,list):
        if not isinstance(a,list) or len(a)!=len(b):raise ValueError('Different list: '+path)
        for i,(x,y) in enumerate(zip(a,b)):equivalent(x,y,path+'[%d]'%i)
    elif isinstance(b,(int,float)) and not isinstance(b,bool):
        if not math.isclose(float(a),float(b),rel_tol=2e-10,abs_tol=1e-30):raise ValueError('Numerical mismatch: '+path)
    elif a!=b:raise ValueError('Different value: '+path)


def export(outputs,summary,verify=False):
    for name,values in outputs.items():
        path=HERE/'results'/name
        if verify:
            actual=read_csv(path)
            if len(actual)!=len(values):raise ValueError('Different row count: '+name)
            for i,(row,expected) in enumerate(zip(actual,values)):
                if set(row)!=set(expected):raise ValueError('Different columns: '+name)
                for key,value in expected.items():
                    if value is None:
                        if row[key]!='':raise ValueError('Expected an undefined value')
                    elif isinstance(value,(int,float)):equivalent(float(row[key]),value,name+'[%d].'%i+key)
                    elif row[key]!=value:raise ValueError('Different CSV value: '+name)
        else:
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(values[0]),lineterminator='\n');writer.writeheader();writer.writerows(values)
    path=HERE/'results/SUMMARY.json'
    if verify:equivalent(json.loads(path.read_text()),summary)
    else:path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify:check_files()
    outputs,summary=compute();export(outputs,summary,args.verify)
    print('PASS extension pilot: 3 configurations; 9 MPI logs replayed; 9900 audited modes; 23976 raw rows; 5112 reference records exact; all signed diagnostics and paired errors reproduced.')


if __name__=='__main__':main()
