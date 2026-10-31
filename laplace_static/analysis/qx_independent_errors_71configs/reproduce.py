#!/usr/bin/env python3
"""Independent-configuration error propagation for the complete qX sample."""
import argparse
import csv
import hashlib
import importlib.util
import json
import tempfile
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
UPSTREAM=HERE.parent/'qx_complete_ensemble_71configs'
RESULTS=HERE/'results'
KEYS=('tau','axis','R','long','tr1','tr2')
SPATIAL=KEYS[1:]
CONFIGS=list(range(100,801,10))
PITERS=(4,40)
TAUS=(2,4,6)
KINDS=('RAW','TRACE','TL','fraction_percent')
EXPORTS=('profiles_independent_errors_71configs.csv','temporal_independent_errors_71configs.csv',
         'temporal_independent_covariance_71configs.csv','R4_midpoint_71configs.csv','SUMMARY.json')


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,str(path));m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m);return m


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,default=str,allow_nan=False)+'\n')


def load():
    metadata=json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name,digest in metadata['repository_references'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Preserved reference changed: '+name)
    upstream=module('qx_independent_source',UPSTREAM/'reproduce.py')
    _,_,source,_,_=upstream.dependencies()
    groups=upstream.grouped([(upstream.TRACE70/'inputs/profiles_raw_and_trace_70configs.csv',PITERS),
        (upstream.EXTENSION/'profiles_raw_and_trace.csv',PITERS)],CONFIGS,PITERS,source.expected_geometry(),trace=True)
    archived={tuple(int(r[k]) for k in KEYS):r for r in read_csv(UPSTREAM/'results/profiles_raw_and_traceless_71configs.csv')}
    temporal={tuple(int(r[k]) for k in SPATIAL)+(r['definition'],int(r['piter'])):r for r in
        read_csv(UPSTREAM/'results/temporal_comparison_71configs.csv')}
    return metadata,groups,archived,temporal


def channel(values,kind):
    """Ratio-of-means point, gradient and synchronized single-sample fluctuations."""
    n=len(values)
    if n<2 or any(len(v)!=5 for v in values):raise ValueError('Need matched five-primary samples')
    with localcontext() as ctx:
        ctx.prec=50
        totals=[sum((v[j] for v in values),D(0)) for j in range(5)]
        mean=[v/n for v in totals]
        L,LS,S,T,LT=mean
        if not L:raise ValueError('Zero loop mean')
        raw=totals[1]/totals[0]-totals[2]/n
        trace=totals[4]/totals[0]-totals[3]/n
        gr=(-LS/(L*L),1/L,D(-1),D(0),D(0))
        gt=(-LT/(L*L),D(0),D(0),D(-1),1/L)
        if kind=='RAW':point,gradient=raw,gr
        elif kind=='TRACE':point,gradient=trace,gt
        elif kind=='TL':point,gradient=(totals[1]-totals[4])/totals[0]-(totals[2]-totals[3])/n,tuple(x-y for x,y in zip(gr,gt))
        elif kind=='fraction_percent':
            if not raw:raise ValueError('Zero connected raw profile in trace fraction')
            point=100*trace/raw
            gradient=tuple(100*(t*raw-trace*r)/(raw*raw) for r,t in zip(gr,gt))
        else:raise ValueError('Unknown observable')
        delta=[sum((gradient[j]*(v[j]-mean[j]) for j in range(5)),D(0)) for v in values]
        return point,delta


def covariance(x,y):
    if len(x)!=len(y) or len(x)<2:raise ValueError('Covariance requires matched samples')
    with localcontext() as ctx:
        ctx.prec=50;n=len(x);mx=sum(x,D(0))/n;my=sum(y,D(0))/n
        return sum(((a-mx)*(b-my) for a,b in zip(x,y)),D(0))/(n*(n-1))


def se(delta):
    with localcontext() as ctx:
        ctx.prec=50
        return covariance(delta,delta).sqrt()


def numerical_exports(output=RESULTS):
    metadata,groups,archived,old_temporal=load()
    profiles=[];fluctuations={};points={}
    with localcontext() as ctx:
        ctx.prec=50
        for key,group in sorted(groups.items()):
            row=dict(zip(KEYS,key));row['configurations']=71
            for p in PITERS:
                samples=[group[c,p] for c in CONFIGS]
                for kind in KINDS:
                    point,delta=channel(samples,kind)
                    field=('trace_fraction_percent' if kind=='fraction_percent' else 'C_'+kind)+'_p'+str(p)
                    if point!=D(archived[key][field]):raise ValueError('Archived point estimate changed: '+field)
                    naive=se(delta);jk=D(archived[key][field+'_jackknife_se'])
                    row[field]=point;row[field+'_naive_se']=naive;row[field+'_jackknife_se']=jk
                    row[field+'_naive_over_jackknife']=naive/jk if jk else ''
                    if kind in ('RAW','TL'):
                        fluctuations[key,kind,p]=delta;points[key,kind,p]=point
            for kind in ('RAW','TL'):
                field='delta_C_'+kind+'_p4_minus_p40'
                delta=[a-b for a,b in zip(fluctuations[key,kind,4],fluctuations[key,kind,40])]
                naive=se(delta);jk=D(archived[key][field+'_jackknife_se'])
                row[field]=D(archived[key][field]);row[field+'_naive_se']=naive
                row[field+'_jackknife_se']=jk;row[field+'_naive_over_jackknife']=naive/jk if jk else ''
            profiles.append(row)
        temporal=[];covrows=[]
        for geometry in sorted({k[1:] for k in groups}):
            for kind in ('RAW','TL'):
                for p in PITERS:
                    row=dict(zip(SPATIAL,geometry));row.update(definition=kind,piter=p,configurations=71)
                    ref=old_temporal[geometry+(kind,p)]
                    for t in TAUS:
                        field='C_tau'+str(t);key=(t,)+geometry
                        naive=se(fluctuations[key,kind,p]);jk=D(ref[field+'_jackknife_se'])
                        row[field]=points[key,kind,p];row[field+'_naive_se']=naive
                        row[field+'_jackknife_se']=jk;row[field+'_naive_over_jackknife']=naive/jk if jk else ''
                    for high,low in ((4,2),(6,4),(6,2)):
                        field='delta_tau%d_minus_tau%d'%(high,low)
                        point=points[(high,)+geometry,kind,p]-points[(low,)+geometry,kind,p]
                        delta=[a-b for a,b in zip(fluctuations[(high,)+geometry,kind,p],fluctuations[(low,)+geometry,kind,p])]
                        naive=se(delta);jk=D(ref[field+'_jackknife_se'])
                        if point!=D(ref[field]):raise ValueError('Temporal point changed')
                        row[field]=point;row[field+'_naive_se']=naive
                        row[field+'_jackknife_se']=jk;row[field+'_naive_over_jackknife']=naive/jk if jk else ''
                    temporal.append(row)
                    for t in TAUS:
                        for u in TAUS:
                            cov=covariance(fluctuations[(t,)+geometry,kind,p],fluctuations[(u,)+geometry,kind,p])
                            errors=row['C_tau%d_naive_se'%t]*row['C_tau%d_naive_se'%u]
                            result=dict(zip(SPATIAL,geometry));result.update(definition=kind,piter=p,tau_i=t,tau_j=u,
                                covariance=cov,correlation=cov/errors if errors else '')
                            covrows.append(result)
        midpoint=[r for r in temporal if tuple(r[k] for k in SPATIAL)==(-1,4,2,0,0)]
        ranges=[]
        for kind in ('RAW','TL'):
            field='C_'+kind+'_p40_naive_over_jackknife'
            ratios=[r[field] for r in profiles if r['axis']==-1]
            ranges.append(dict(definition=kind,piter=40,combined_positions=len(ratios),minimum=min(ratios),maximum=max(ratios)))
        summary=dict(success=True,configurations=CONFIGS,configuration_count=71,cases=213,
            profile_rows=len(profiles),temporal_rows=len(temporal),covariance_entries=len(covrows),
            ensemble_assessment=metadata['ensemble_assessment'],
            statistical_unit='One spatially averaged configuration. N=71, not number of sources, probes, axes or lattice sites.',
            baseline='Independent-configuration linear propagation using full covariance; sample covariance divided by N.',
            estimator='Ratio of raw ensemble means, not the mean of per-configuration ratios; no bias correction.',
            formula='delta_i=gradient_f dot (primary_i-primary_mean); covariance(mean f,mean g)=sum((delta_f_i-mean_delta_f)*(delta_g_i-mean_delta_g))/(N*(N-1)).',
            resampling_control='Existing paired delete-one-configuration jackknife retained as a nonlinear-estimator cross-check.',
            pyerrors_setting='gamma_method(S=0) is the matching independent-sample setting; automatic differentiation is checked separately.',
            assumption='The supplied saved configurations are thermalized and approximately statistically independent, as established in the ensemble discussion.',
            naive_over_jackknife_ranges=ranges,R4_midpoint=midpoint,
            scope='RAW and traceless p4/p40 matched full71 inputs. p10/p20 projection conclusions remain in the complete-ensemble study. Historical binning outputs are retained unchanged as diagnostics, not imposed as error-inflation factors.')
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    for name,table in zip(EXPORTS[:4],(profiles,temporal,covrows,midpoint)):write_csv(out/name,table)
    write_json(out/EXPORTS[4],summary)
    return profiles,temporal,summary


def verify():
    with tempfile.TemporaryDirectory() as tmp:
        numerical_exports(tmp)
        for name in EXPORTS:
            if (Path(tmp)/name).read_bytes()!=(RESULTS/name).read_bytes():raise ValueError('Numerical export differs: '+name)
    import math
    baseline={tuple(int(r[k]) for k in KEYS):r for r in read_csv(RESULTS/EXPORTS[0])}
    checks=read_csv(RESULTS/'pyerrors_S0_checks_71configs.csv')
    if len(checks)!=996:raise ValueError('Incomplete archived PyErrors profile check')
    for row in checks:
        ref=baseline[tuple(int(row[k]) for k in KEYS)]
        for kind in ('RAW','TL'):
            expected=float(ref['C_'+kind+'_p40_naive_se'])
            if not math.isclose(float(row[kind+'_pyerrors_S0_se']),expected,rel_tol=2e-8,abs_tol=2e-15):
                raise ValueError('Archived PyErrors error differs from independent propagation')
    py=json.loads((RESULTS/'PYERRORS_CHECK.json').read_text())
    if not py['success'] or py['configuration_labels']!=CONFIGS or py['profile_error_checks']!=1992 or py['temporal_contrast_checks']!=1992:
        raise ValueError('Incomplete archived PyErrors audit')
    print('PASS independent full71 errors: all archived points exact; 996 profile rows; 1328 temporal rows; 11952 covariance entries; paired nonlinear contrasts preserved.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');p.add_argument('--plots',action='store_true')
    args=p.parse_args()
    if args.verify:verify()
    else:
        profiles,temporal,summary=numerical_exports()
        if args.plots:module('qx_independent_plots',HERE/'render_plots.py').render(profiles,temporal,RESULTS)
        print('PASS independent-configuration numerical exports')


if __name__=='__main__':main()
