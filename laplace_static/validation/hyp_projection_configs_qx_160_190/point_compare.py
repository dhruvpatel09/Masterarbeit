#!/usr/bin/env python3
import math
from pathlib import Path
import re
import sys

num=re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
def read(path,t0=96,t1=100):
    s=path.read_text()
    if re.search(r'(?<![A-Za-z_])[-+]?(nan|inf)(?![A-Za-z_])',s,re.I):
        raise RuntimeError(f'Nonfinite data: {path}')
    if 'status=FAIL' in s:
        raise RuntimeError(f'Factorization failed: {path}')
    expected={'baseline':(3,40),'off':(0,40),'zero':(1,40),'hyp40':(2,40),'hyp80':(2,80)}
    if path.stem in expected:
        mode,piter=expected[path.stem]
        tags=[x for x in s.splitlines() if x.startswith('HYP_INTEGRATION ')]
        if len(tags)!=1 or f'mode={mode} ' not in tags[0] or f'piter={piter} ' not in tags[0]:
            raise RuntimeError(f'Wrong HYP mode/iteration count: {path}')
    records={}
    for prefix,expected in [('TAU ',10000),('SELECTED ',16),('SUMMARY ',1),('FACTOR_R0 ',1),('PAIR_R1 ',2),('PROBE ',3),('OBSERVABLE ',2)]:
        rows=sorted(x for x in s.splitlines() if x.startswith(prefix))
        if len(rows)!=expected:
            raise RuntimeError(f'{path}: {prefix} expected {expected}, got {len(rows)}')
        for i,row in enumerate(rows):
            if prefix in ('FACTOR_R0 ','PAIR_R1 ') and 'status=PASS' not in row:
                raise RuntimeError('Missing factorization PASS')
            # Compare physical results; separately enforce factorization status.
            row=re.sub(r'\s+(?:absdiff|reldiff)=\S+','',row)
            records[(prefix,i)]=row
    pairs='\n'.join(records[k] for k in records if k[0] in ('PAIR_R1 ','OBSERVABLE '))
    if 'kind=cross_boundary' not in pairs or 'y3_halo=1' not in pairs or 'y4_halo=1' not in pairs:
        raise RuntimeError('Missing MPI halo pair coverage')
    paths=[x for x in s.splitlines() if x.startswith('PATH_CHECK ')]
    if len(paths)!=1 or f't0={t0} tmid={(t0+t1)//2} t1={t1}' not in paths[0] or 'status=PASS' not in paths[0]:
        raise RuntimeError(f'Missing/wrong path check: {path}')
    def scalar(row,key):
        return float(re.search(r'\b'+key+r'=('+num.pattern+')',row).group(1))
    for key,row in records.items():
        if key[0]!='OBSERVABLE ': continue
        kind=re.search(r'kind=(\S+)',row).group(1)
        probe=next(v for k,v in records.items() if k[0]=='PROBE ' and 'kind='+kind+' ' in v)
        if scalar(row,'S_RAW')!=scalar(probe,'S_RAW'):
            raise RuntimeError('Probe/observable mismatch')
        for part in ('RE','IM'):
            expected=scalar(row,'L_'+part)*scalar(row,'S_RAW')
            actual=scalar(row,'LS_'+part)
            if abs(actual-expected)>1e-25+2e-15*abs(expected):
                raise RuntimeError('LS multiplication mismatch')
    return records

def compare(a,b,label):
    worst=0.
    if a.keys()!=b.keys(): raise RuntimeError('Record keys differ')
    for key in a:
        x,y=a[key],b[key]
        if num.sub('#',x)!=num.sub('#',y): raise RuntimeError(f'{label}: record shape {key}')
        for sx,sy in zip(num.findall(x),num.findall(y)):
            vx,vy=float(sx),float(sy)
            d=abs(vx-vy); worst=max(worst,d)
            if not math.isfinite(d) or d>1e-18+1e-11*max(abs(vx),abs(vy)):
                raise RuntimeError(f'{label}: {key}: {sx} vs {sy}')
    print(f'PASS {label}: max_absolute_difference={worst:.16e}')

def main():
    root=Path(sys.argv[1])
    data={n:read(root/(n+'.out')) for n in ('baseline','off','zero','hyp40','hyp80')}
    compare(data['baseline'],data['off'],'baseline vs disabled')
    compare(data['baseline'],data['zero'],'baseline vs zero-alpha')
    compare(data['hyp40'],data['hyp80'],'HYP piter40 vs piter80')
    for name,records in data.items():
        a={k:v for k,v in data['baseline'].items() if k[0]=='PROBE '}
        b={k:v for k,v in records.items() if k[0]=='PROBE '}
        compare(a,b,'original probe vs '+name)
    for name,records in data.items():
        print(name.upper(),records[('SUMMARY ',0)])
        for k in records:
            if k[0] in ('PAIR_R1 ','OBSERVABLE '): print(name.upper(),records[k])
    print('PASS qX multi-link temporal-HYP and original-field clover integration')
if __name__=='__main__': main()
