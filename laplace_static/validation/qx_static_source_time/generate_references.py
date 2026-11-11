#!/usr/bin/env python3
"""Freeze sources outside the pilot and export the already measured reference."""
import argparse,csv,hashlib,importlib.util,json,math,sys
import numpy as np
from plan import HERE,ROOT,BASE_COMMIT,CONFIGS,PILOT,RADII,CUTOFFS,TAUS,BASES,ANCHORS

OLD=ROOT/'laplace_static/analysis/qx_static_matrix_71configs'


def published_source_algorithm():
    previous=sys.path[:];old=sys.modules.pop('statistics_tools',None);sys.path.insert(0,str(OLD))
    try:
        spec=importlib.util.spec_from_file_location('qx_source_time_fixed_training',str(OLD/'projections.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
    finally:
        sys.path[:]=previous;sys.modules.pop('statistics_tools',None)
        if old is not None:sys.modules['statistics_tools']=old


def compute():
    pins=json.loads((HERE/'inputs/SOURCE_PINS.json').read_text())
    for name,r in pins['dependencies'].items():
        p=ROOT/name
        if p.stat().st_size!=r['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=r['sha256']:raise ValueError('Published source changed: '+name)
    cache=np.load(str(OLD/'inputs/primary_matrices.npz'),allow_pickle=False)
    if (list(cache['configurations'])!=list(CONFIGS) or list(cache['radii'])!=list(RADII)
            or list(cache['cutoffs'])!=list(CUTOFFS) or cache['C'].shape!=(71,4,33,4,4)):raise ValueError('Primary cache geometry')
    C=cache['C'];native=published_source_algorithm();training=np.array([c not in PILOT for c in CONFIGS]);records=[]
    for ri,r in enumerate(RADII):
        for basis,cutoffs in BASES.items():
            indices=[CUTOFFS.index(n) for n in cutoffs];mean=native.submatrix(C[training,ri],indices).mean(0)
            for anchor in ANCHORS:
                v,diagnostics=native.train(mean,anchor)
                records.append(dict(R=r,basis=basis,anchor=anchor,Nv=list(cutoffs),
                    coefficients_re=list(map(float,v.real)),coefficients_im=list(map(float,v.imag)),
                    reference_metric_norm=float(np.vdot(v,mean[0].dot(v)).real),diagnostics=diagnostics))
    vectors=dict(base_commit=BASE_COMMIT,training_t0=90,training_configuration_labels=[c for c in CONFIGS if c not in PILOT],
        excluded_pilot_configuration_labels=list(PILOT),source_candidates=records,
        interpretation='Fixed pilot controls conditional on the separate 68-configuration training set. These vectors are not retuned at shifted times and are not final physical GEVP levels.')
    baseline=[]
    for cfg in PILOT:
        ci=CONFIGS.index(cfg)
        for ri,r in enumerate(RADII):
            for tau in TAUS:
                for source,nsource in enumerate(CUTOFFS):
                    for sink,nsink in enumerate(CUTOFFS):
                        z=C[ci,ri,tau,sink,source]
                        baseline.append(dict(configuration=cfg,t0=90,t1=90+tau,tau=tau,R=r,axis=-1,piter=10,
                            Nv_src=nsource,Nv_sink=nsink,Nsrc=3*64**3,L_RE=float(z.real),L_IM=float(z.imag)))
    return vectors,baseline


def compare(a,b,label):
    if isinstance(b,dict):
        if set(a)!=set(b):raise ValueError('Reference keys '+label)
        for key in b:compare(a[key],b[key],label+'.'+key)
    elif isinstance(b,list):
        if len(a)!=len(b):raise ValueError('Reference length '+label)
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,label+'['+str(i)+']')
    elif isinstance(b,float):
        if not math.isfinite(float(a)) or not math.isfinite(b) or not math.isclose(float(a),b,rel_tol=5e-7,abs_tol=5e-11):raise ValueError('Reference number '+label)
    elif str(a)!=str(b):raise ValueError('Reference label '+label)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    vectors,baseline=compute();vp=HERE/'inputs/FIXED_SOURCES.json';bp=HERE/'inputs/reference_t90_axis_average.csv'
    if args.verify:
        compare(json.loads(vp.read_text()),vectors,'fixed_sources')
        with bp.open(newline='') as stream:old=list(csv.DictReader(stream))
        if len(old)!=len(baseline):raise ValueError('Reference row inventory')
        for i,(a,b) in enumerate(zip(old,baseline)):
            if set(a)!=set(b):raise ValueError('Reference row fields '+str(i))
            for key,value in b.items():
                if key in ('L_RE','L_IM'):
                    if float(a[key])!=value:raise ValueError('Raw reference differs from primary cache at row {} {}'.format(i,key))
                else:compare(a[key],value,'baseline row {} {}'.format(i,key))
        print('PASS frozen sources: 68 training configurations, three excluded pilot configurations, 36 source candidates and 6336 reference rows.')
    else:
        vp.write_text(json.dumps(vectors,indent=2,allow_nan=False)+'\n')
        with bp.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(baseline[0]),lineterminator='\n');writer.writeheader();writer.writerows(baseline)
        print('Generated fixed pilot sources and existing t0=90 reference; no new measurements.')


if __name__=='__main__':main()
