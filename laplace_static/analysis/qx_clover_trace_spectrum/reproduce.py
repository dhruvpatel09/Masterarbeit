#!/usr/bin/env python3
"""Replay and compare the matched four-configuration projection diagnostic."""
import argparse
import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import sys
import tempfile
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
JOB=HERE/'inputs/job5944679'
RESULTS=HERE/'results'
CONFIGS=[100,170,670,740]
RAW_PITERS=[4,6,8,10,20,40,80]
TL_PITERS=[4,6,8,40,80]
NEW_PITERS=[6,8,80]
KEYS=('tau','axis','R','long','tr1','tr2')
IDENT=('configuration','tau','piter','axis','R','long','tr1','tr2')
D=Decimal


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def geometry():
    offsets={(r,l,x,z) for r in (1,2,4)
             for l,x,z in ({(0,1,2)}|{(l,x,0) for l in range(-1,r+2) for x in range(-2,3)})}
    return {(t,a,r,l,x,z) for t in (2,4,6) for a in (-1,0,1,2) for r,l,x,z in offsets}


def verify_inputs():
    sources=json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name,expected in sources['files'].items():
        if sha(HERE/name)!=expected:raise ValueError('Preserved input changed: '+name)
    for name,expected in sources['repository_references'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Preserved reference changed: '+name)
    for local,relative in [('ensemble_INPUTS.json','qx_clover_trace_ensemble/INPUTS.json'),
                           ('pilot_INPUTS.json','qx_clover_trace/INPUTS.json')]:
        if (HERE/'inputs'/local).read_bytes()!=(ROOT/'laplace_static/validation'/relative).read_bytes():
            raise ValueError('Prepared source pins differ from the preserved workflow')
    plan=json.loads((JOB/'PREPARED.json').read_text())
    expected=[dict(index=i,configuration=c,tau=t,piters=NEW_PITERS)
              for i,(c,t) in enumerate((c,t) for c in CONFIGS for t in (2,4,6))]
    if plan['cases']!=expected or plan['slots']!=[[i] for i in range(12)] or plan['stage']!='spectrum':
        raise ValueError('Wrong diagnostic case plan')
    if sha(JOB/'PREPARED.json')!=(JOB/'prepared.sha256').read_text().strip():
        raise ValueError('Preparation checksum differs')
    for key,name in [('source_manifest_sha256','ensemble_INPUTS.json'),
                     ('pilot_manifest_sha256','pilot_INPUTS.json'),('build_manifest_sha256','build_manifest.sha256')]:
        if plan[key]!=sha(HERE/'inputs'/name):raise ValueError('Prepared provenance differs: '+key)
    summary=json.loads((JOB/'SUMMARY.json').read_text())
    if (summary['success'],summary['stage'],summary['target_configurations'],summary['missing_cases'],
            summary['raw_records'],summary['window_records'])!=(True,'spectrum',CONFIGS,[],11952,36):
        raise ValueError('Incomplete diagnostic collection')
    return sources,plan


def log_bytes(path,expected):
    body=gzip.decompress(Path(path).read_bytes())
    if hashlib.sha256(body).hexdigest()!=expected:raise ValueError('Decompressed stdout changed')
    return body


def replay(sources,plan):
    parser=ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230'
    sys.path.insert(0,str(parser))
    import check_task as check
    spec=importlib.util.spec_from_file_location('clover_spectrum_trace_parser',
        str(ROOT/'laplace_static/validation/qx_clover_trace/campaign.py'))
    pilot=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    raw,windows,checks,launches_all=[],[],[],[]
    with tempfile.TemporaryDirectory() as temp:
        temp=Path(temp)
        for entry in plan['cases']:
            index,cfg,tau=entry['index'],entry['configuration'],entry['tau']
            case=JOB/('case_%03d'%index)
            latest=json.loads((case/'LATEST.json').read_text())
            if latest['attempt']!='attempt_0001':raise ValueError('Unexpected preserved attempt')
            attempt=case/latest['attempt']
            if sha(attempt/'SUMMARY.json')!=latest['summary_sha256']:raise ValueError('Completion checksum differs')
            if json.loads((attempt/'SUMMARY.json').read_text())!=dict(success=True,index=index,
                    configuration=cfg,tau=tau,piters=NEW_PITERS,raw_records=996,window_records=3):
                raise ValueError('Incorrect case completion metadata')
            marker=json.loads((attempt/'CASE.json').read_text())
            if marker['prepared_sha256']!=sha(JOB/'PREPARED.json'):raise ValueError('Case preparation differs')
            launches=json.loads((attempt/'LAUNCHES.json').read_text())
            if [(r['configuration'],r['tau'],r['piter'],r['exit_code']) for r in launches]!=[(cfg,tau,p,0) for p in NEW_PITERS]:
                raise ValueError('Incomplete MPI launch inventory')
            if not all(math.isfinite(r['seconds']) and r['seconds']>0 for r in launches):
                raise ValueError('Invalid launch timing')
            identity=[]
            for p,launch in zip(NEW_PITERS,launches):
                if json.loads((attempt/('launch_p%d.json'%p)).read_text())!=launch:
                    raise ValueError('Per-launch record differs')
                if (attempt/('p%d.err'%p)).read_bytes():raise ValueError('Nonempty stderr')
                compressed=attempt/('p%d.out.gz'%p)
                body=log_bytes(compressed,sources['original_stdout_sha256'][compressed.relative_to(HERE).as_posix()])
                path=temp/('p%d.out'%p)
                path.write_bytes(body)
                with contextlib.redirect_stdout(io.StringIO()):
                    _,scan,_,window=check.read_measurement(path,cfg,tau,p)
                if p==80 and max(float(window[k]) for k in ('unit_max','det_max'))>1e-12:
                    raise ValueError('High-iteration SU(3) control failed')
                text=body.decode()
                trace=pilot.trace_records(text,scan,check.avg)
                identity.append(float(check.avg.rows(text,'CLOVER_TRACE_CHECK')[0]['identity_max']))
                for row in check.raw_rows(cfg,tau,p,scan):
                    key=tuple(int(row[k]) for k in ('axis','R','long','tr1','tr2'))
                    row.update({k:trace[key][k] for k in pilot.TRACE_FIELDS})
                    raw.append({k:str(v) for k,v in row.items()})
                windows.append({k:str(v) for k,v in window.items()})
            launches_all+=launches
            checks.append(dict(configuration=cfg,tau=tau,success=True,
                raw_records=996,window_records=3,direct_clover_identity_max=max(identity)))
    raw.sort(key=lambda r:tuple(int(r[k]) for k in IDENT))
    windows.sort(key=lambda r:tuple(int(r[k]) for k in ('configuration','tau','piter')))
    if raw!=read_csv(JOB/'profiles_raw_and_trace.csv') or windows!=read_csv(JOB/'su3_windows.csv'):
        raise ValueError('Reparsed primary table differs from campaign collection')
    report=dict(success=True,job_id=5944679,configurations=CONFIGS,cases=12,mpi_launches=36,
        raw_records=len(raw),window_records=len(windows),
        maximum_direct_clover_identity_residual=max(r['direct_clover_identity_max'] for r in checks),
        checks=checks,launches=launches_all,
        scope='Complete offline replay of all 36 preserved stdout logs, case completion hashes, launch records, profile geometry, rank reductions, halo coordinates, direct samples, trace identities and collected primary tables. No gauge or saved-eigenvector files loaded.')
    return raw,windows,report


def selected_rows(path):
    with Path(path).open(newline='') as f:
        return [r for r in csv.DictReader(f) if int(r['configuration']) in CONFIGS]


def group_rows(rows,piters,fields):
    result=defaultdict(dict)
    expected_geometry=geometry()
    for row in rows:
        cfg,p=int(row['configuration']),int(row['piter'])
        key=tuple(int(row[k]) for k in KEYS)
        if cfg not in CONFIGS or p not in piters or key not in expected_geometry:
            raise ValueError('Unexpected configuration, piter or geometry')
        tau,axis=key[:2]
        if (int(row['Nv']),int(row['t0']),int(row['t1']),int(row['tmid']),int(row['Nsrc']))!=(
                100,98-tau//2,98+tau//2,98,64**3*(3 if axis==-1 else 1)):
            raise ValueError('Incorrect measurement metadata')
        if (cfg,p) in result[key]:raise ValueError('Duplicate configuration/piter/profile record')
        values=tuple(D(row[f]) for f in fields)
        if not all(v.is_finite() for v in values) or values[0]==0:raise ValueError('Invalid primary')
        result[key][cfg,p]=values
    if set(result)!=expected_geometry:raise ValueError('Incomplete profile geometry')
    expected={(c,p) for c in CONFIGS for p in piters}
    for key,values in result.items():
        if set(values)!=expected:raise ValueError('Incomplete matched cohort at '+str(key))
        for c in CONFIGS:
            if len({values[c,p][2] for p in piters})!=1:raise ValueError('Original clover changed across piters')
            if len(fields)==5 and len({values[c,p][3] for p in piters})!=1:
                raise ValueError('Original trace density changed across piters')
    return result


def evaluate(raw_sums,trace_sums,n):
    estimates={}
    for p,(L,LS,S) in raw_sums.items():
        if L==0:raise ValueError('Undefined loop ratio')
        estimates['C_RAW_p%d'%p]=LS/L-S/D(n)
    for p,(L,LS,S,T,LT) in trace_sums.items():
        if L==0:raise ValueError('Undefined loop ratio')
        estimates['C_TL_p%d'%p]=(LS-LT)/L-(S-T)/D(n)
        estimates['C_TRACE_p%d'%p]=LT/L-T/D(n)
    for definition,piters in [('RAW',sorted(raw_sums)),('TL',sorted(trace_sums))]:
        ref=estimates['C_%s_p40'%definition]
        for p in piters:
            if p==40:continue
            delta=estimates['C_%s_p%d'%(definition,p)]-ref
            estimates['delta_C_%s_p%d_minus_p40'%(definition,p)]=delta
            estimates['relative_delta_percent_%s_p%d'%(definition,p)]=100*delta/ref if ref else None
    return estimates


def paired_statistics(raw,trace,configs):
    with localcontext() as context:
        context.prec=50
        n=len(configs)
        if n<2:raise ValueError('Need at least two configurations for errors')
        raw_p=sorted({p for c,p in raw})
        trace_p=sorted({p for c,p in trace})
        if set(raw)!={(c,p) for c in configs for p in raw_p} or set(trace)!={(c,p) for c in configs for p in trace_p}:
            raise ValueError('Unmatched jackknife sample')
        sums=lambda data,piters,width:{p:tuple(sum((data[c,p][j] for c in configs),D(0)) for j in range(width)) for p in piters}
        rr,tt=sums(raw,raw_p,3),sums(trace,trace_p,5)
        point=evaluate(rr,tt,n)
        leave=[]
        for c in configs:
            r={p:tuple(rr[p][j]-raw[c,p][j] for j in range(3)) for p in raw_p}
            t={p:tuple(tt[p][j]-trace[c,p][j] for j in range(5)) for p in trace_p}
            leave.append(evaluate(r,t,n-1))
        result={k:'' if v is None else v for k,v in point.items()}
        for field in point:
            values=[r[field] for r in leave]
            if point[field] is None or any(v is None for v in values):
                result[field+'_jackknife_se']=''
            else:
                center=sum(values,D(0))/D(n)
                result[field+'_jackknife_se']=(D(n-1)/D(n)*sum(((v-center)**2 for v in values),D(0))).sqrt()
        return result


def compare(records):
    combined=[r for r in records if r['axis']==-1]
    out=[]
    with localcontext() as context:
        context.prec=50
        for definition,piters in [('RAW',RAW_PITERS),('TL',TL_PITERS)]:
            for p in piters:
                delta=lambda r:D(0) if p==40 else r['delta_C_%s_p%d_minus_p40'%(definition,p)]
                eligible=[r for r in combined if r['C_%s_p40'%definition]!=0]
                absolute=max(combined,key=lambda r:abs(delta(r)))
                relative=max(eligible,key=lambda r:abs(delta(r)/r['C_%s_p40'%definition]))
                row=dict(definition=definition,piter=p,configurations=4,combined_positions=249,
                    max_abs_delta_C_vs_p40=str(abs(delta(absolute))),
                    max_relative_delta_percent_vs_p40=str(100*abs(delta(relative)/relative['C_%s_p40'%definition])))
                for prefix,source in [('max_absolute',absolute),('max_relative',relative)]:
                    row.update({prefix+'_'+k:source[k] for k in ('tau','R','long','tr1','tr2')})
                out.append(row)
    return out


def per_configuration(raw_groups,trace_groups):
    rows=[]
    with localcontext() as context:
        context.prec=50
        for cfg in CONFIGS:
            for tau in (2,4,6):
                estimates=[]
                for key in sorted(raw_groups):
                    if key[0]!=tau or key[1]!=-1:continue
                    rr={p:raw_groups[key][cfg,p] for p in RAW_PITERS}
                    tt={p:trace_groups[key][cfg,p] for p in TL_PITERS}
                    estimates.append(dict(zip(KEYS,key),**evaluate(rr,tt,1)))
                for definition,piters in [('RAW',RAW_PITERS),('TL',TL_PITERS)]:
                    for p in piters:
                        if p==40:continue
                        field='delta_C_%s_p%d_minus_p40'%(definition,p)
                        most=max(estimates,key=lambda r:abs(r[field]))
                        nonzero=[r for r in estimates if r['C_%s_p40'%definition]!=0]
                        relative=max(nonzero,key=lambda r:abs(r[field]/r['C_%s_p40'%definition]))
                        rows.append(dict(configuration=cfg,tau=tau,definition=definition,piter=p,
                            max_absolute_delta_C=str(abs(most[field])),
                            max_relative_delta_percent=str(100*abs(relative[field]/relative['C_%s_p40'%definition])),
                            max_absolute_R=most['R'],max_absolute_long=most['long'],max_absolute_tr1=most['tr1'],max_absolute_tr2=most['tr2']))
    return rows


def export(output):
    sources,plan=verify_inputs()
    new,windows,validation=replay(sources,plan)
    baseline=selected_rows(ROOT/sources['raw_reference'])
    dual=selected_rows(ROOT/sources['dual_reference'])
    old={(tuple(int(r[k]) for k in IDENT)):r for r in baseline}
    for row in dual:
        ref=old[tuple(int(row[k]) for k in IDENT)]
        if any(D(row[k])!=D(ref[k]) for k in ('L_RE','L_IM','S_RAW','LS_RE','LS_IM')):
            raise ValueError('Raw and dual-clover controls differ')
    raw=group_rows(baseline+new,RAW_PITERS,('L_RE','LS_RE','S_RAW'))
    traced=group_rows(dual+new,TL_PITERS,('L_RE','LS_RE','S_RAW','T_TRACE','LT_RE'))
    records=[]
    for key in sorted(raw):
        row=dict(zip(KEYS,key));row['configurations']=4
        row.update(paired_statistics(raw[key],traced[key],CONFIGS))
        records.append(row)
    windows=selected_rows(ROOT/sources['window_reference'])+windows
    windows.sort(key=lambda r:tuple(int(r[k]) for k in ('configuration','tau','piter')))
    if len(windows)!=84 or {(int(r['configuration']),int(r['tau']),int(r['piter'])) for r in windows}!={(c,t,p) for c in CONFIGS for t in (2,4,6) for p in RAW_PITERS}:
        raise ValueError('Incomplete matched link diagnostics')
    comparisons=compare(records)
    individual=per_configuration(raw,traced)
    window_summary=[]
    for p in RAW_PITERS:
        selected=[r for r in windows if int(r['piter'])==p]
        unit=max(selected,key=lambda r:D(r['unit_max']))
        det=max(selected,key=lambda r:D(r['det_max']))
        window_summary.append(dict(piter=p,configurations=4,windows=12,
            unit_max=unit['unit_max'],unit_max_configuration=unit['configuration'],unit_max_tau=unit['tau'],
            det_max=det['det_max'],det_max_configuration=det['configuration'],det_max_tau=det['tau']))
    runtimes=[]
    for p in NEW_PITERS:
        times=[r['seconds'] for r in validation['launches'] if r['piter']==p]
        runtimes.append(dict(piter=p,launches=len(times),mean_seconds=math.fsum(times)/len(times),
            minimum_seconds=min(times),maximum_seconds=max(times)))
    pair_times=[sum(r['seconds'] for r in validation['launches'] if (r['configuration'],r['tau'])==(c,t)) for c in CONFIGS for t in (2,4,6)]
    summary=dict(success=True,job_id=5944679,configurations=CONFIGS,configuration_count=4,
        cases=12,new_piters=NEW_PITERS,raw_piters=RAW_PITERS,traceless_piters=TL_PITERS,
        new_primary_rows=11952,matched_raw_rows=27888,matched_trace_rows=19920,
        ensemble_geometry_rows=996,combined_geometry_rows=249,
        comparison=comparisons,central_window_defects=window_summary,
        maximum_direct_clover_identity_residual=validation['maximum_direct_clover_identity_residual'],
        runtime=dict(piters=runtimes,mean_MPI_minutes_per_array_task=math.fsum(pair_times)/len(pair_times)/60,
                     minimum_MPI_minutes_per_array_task=min(pair_times)/60,maximum_MPI_minutes_per_array_task=max(pair_times)/60,
                     scope='Measured MPI launch time, including launch I/O; excludes scheduler waiting and wrapper validation'),
        estimator='Equal-configuration ratios of means; precision50. Raw and traceless use the original six-plane clover without an additional factor 1/2.',
        uncertainties='Paired delete-one-configuration jackknife; four configurations only; preliminary without autocorrelation correction. Relative errors recompute the nonlinear ratio in every omission.',
        scope='Selected central-window convergence diagnostic. All piter curves use the same four configurations. The existing full70 study remains the ensemble result. Traceless piter10/20 primaries are unavailable; no full-lattice field-stability or physical-normalization claim.')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for name,rows in [('profiles_ensemble_4configs.csv',records),('piter_comparison_4configs.csv',comparisons),
                     ('su3_central_windows_4configs.csv',windows),('per_configuration_comparison.csv',individual),
                     ('runtime_by_piter.csv',runtimes)]:write_csv(output/name,rows)
    (output/'VALIDATION.json').write_text(json.dumps(validation,indent=2,allow_nan=False)+'\n')
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return records,summary,individual


def verify():
    with tempfile.TemporaryDirectory() as directory:
        export(Path(directory))
        for name in ('profiles_ensemble_4configs.csv','piter_comparison_4configs.csv','su3_central_windows_4configs.csv',
                     'per_configuration_comparison.csv','runtime_by_piter.csv','VALIDATION.json','SUMMARY.json'):
            if (Path(directory)/name).read_bytes()!=(RESULTS/name).read_bytes():raise ValueError('Export differs: '+name)
    print('PASS spectrum: 36 stdout logs replayed; 12 cases; seven raw and five traceless piters; matched four-configuration ratios, paired errors and diagnostics reproduced.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--plots',action='store_true')
    args=parser.parse_args()
    if args.verify:verify()
    else:
        records,summary,individual=export(RESULTS)
        if args.plots:
            import render_plots
            render_plots.render(records,summary,individual,RESULTS)
        print('PASS spectrum numerical exports')


if __name__=='__main__':main()
