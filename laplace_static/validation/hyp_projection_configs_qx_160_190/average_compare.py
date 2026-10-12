#!/usr/bin/env python3
"""Validate averaging, point regressions, controls and projection stability."""
import math
from pathlib import Path
import sys
import point_compare as point

FIELDS=('L_RE','L_IM','S_RAW','LS_RE','LS_IM')
N=64**3

def rows(text,prefix):
    return [dict(token.split('=',1) for token in row.split()[1:])
            for row in text.splitlines() if row.startswith(prefix+' ')]

def near(a,b,label,scale=None):
    if not math.isfinite(a) or not math.isfinite(b):
        raise RuntimeError(label+': nonfinite')
    if abs(a-b)>1e-20+5e-13*(max(abs(a),abs(b)) if scale is None else scale):
        raise RuntimeError(f'{label}: {a:.17e} vs {b:.17e}')

def read_average(path):
    text=path.read_text()
    if text.count('SOURCE_AVERAGING status=PASS')!=1:
        raise RuntimeError('Missing averaging PASS')
    if text.count('HALO_COORDINATES axes=3 status=PASS')!=1:
        raise RuntimeError('Missing coordinate halo PASS')
    checks=rows(text,'DIRECT_SAMPLES')
    if len(checks)!=3 or {int(r['axis']) for r in checks}!={0,1,2}:
        raise RuntimeError('Missing direct samples')
    for r in checks:
        if int(r['count'])!=3 or r['status']!='PASS':
            raise RuntimeError('Failed direct samples')
    partials=rows(text,'RANK_SUM')
    sums=rows(text,'AXIS_SUM')
    avgs=rows(text,'AXIS_AVG')
    combined=rows(text,'COMBINED_AVG')
    if (len(partials),len(sums),len(avgs),len(combined))!=(192,3,3,1):
        raise RuntimeError('Wrong average record counts')
    expected={(rank,axis) for rank in range(64) for axis in range(3)}
    if {(int(r['rank']),int(r['axis'])) for r in partials}!=expected:
        raise RuntimeError('Missing/duplicate rank-axis records')
    if {int(r['axis']) for r in sums}!={0,1,2} or {int(r['axis']) for r in avgs}!={0,1,2}:
        raise RuntimeError('Missing axes')
    vacuum=next(r for r in rows(text,'PROBE') if r['kind']=='spatial_average')
    result={}
    for axis in range(3):
        ps=[r for r in partials if int(r['axis'])==axis]
        total=next(r for r in sums if int(r['axis'])==axis)
        avg=next(r for r in avgs if int(r['axis'])==axis)
        for r in ps:
            if (int(r['Nsrc']),int(r['Ninterior']),int(r['Nhalo']))!=(4096,3840,256) or int(r['Nwrap']) not in (0,256):
                raise RuntimeError('Wrong local counts')
        for key,want in [('Nsrc',N),('Ninterior',245760),('Nhalo',16384),('Nwrap',4096)]:
            if int(total[key])!=want or sum(int(r[key]) for r in ps)!=want:
                raise RuntimeError(f'Wrong counts {axis} {key}')
        if int(avg['Nsrc'])!=N: raise RuntimeError('Wrong average denominator')
        for k,name in enumerate(FIELDS):
            values=[float(r[f'v{k}']) for r in ps]
            near(math.fsum(values),float(total[f'v{k}']),f'rank reduction {axis} {name}',math.fsum(map(abs,values)))
            near(float(total[f'v{k}'])/N,float(avg[name]),f'normalization {axis} {name}')
        near(float(avg['S_RAW']),float(vacuum['S_RAW']),f'probe translation average {axis}')
        result[('AXIS_AVG ',axis)]='AXIS_AVG '+' '.join(f'{k}={v}' for k,v in avg.items())
    allrow=combined[0]
    if int(allrow['Nsrc'])!=3*N: raise RuntimeError('Wrong combined denominator')
    for name in FIELDS:
        vals=[float(r[name]) for r in avgs]
        near(math.fsum(vals)/3,float(allrow[name]),f'combined {name}',math.fsum(map(abs,vals))/3)
    result[('COMBINED_AVG ',0)]='COMBINED_AVG '+' '.join(f'{k}={v}' for k,v in allrow.items())
    print(f'PASS {path.stem}: counts, halo coordinates, direct samples, rank reductions, probe average')
    return result

def main():
    root=Path(sys.argv[1])
    reference=Path(sys.argv[2])
    names=('baseline','off','zero','hyp40','hyp80')
    points={name:point.read(root/(name+'.out')) for name in names}
    avgs={name:read_average(root/(name+'.out')) for name in names}
    for name in names:
        point.compare(point.read(reference/(name+'.out')),points[name],name+' vs validated job 5903809 points')
    for collection,label in [(points,'points'),(avgs,'averages')]:
        point.compare(collection['baseline'],collection['off'],label+' baseline vs off')
        point.compare(collection['baseline'],collection['zero'],label+' baseline vs zero')
        point.compare(collection['hyp40'],collection['hyp80'],label+' piter40 vs piter80')
    for name in names:
        a={k:v for k,v in points['baseline'].items() if k[0]=='PROBE '}
        b={k:v for k,v in points[name].items() if k[0]=='PROBE '}
        point.compare(a,b,'original probes vs '+name)
        for key,row in avgs[name].items():
            print(name.upper(),row)
    print('PASS qX spatial-source and three-axis averaging with temporal HYP and original clover')
if __name__=='__main__': main()
