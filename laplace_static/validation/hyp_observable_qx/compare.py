#!/usr/bin/env python3
import math
from pathlib import Path
import re
import sys

num=re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
def read(path):
    s=path.read_text()
    if re.search(r'(?<![A-Za-z_])[-+]?(nan|inf)(?![A-Za-z_])',s,re.I):
        raise RuntimeError(f'Nonfinite data: {path}')
    if 'status=FAIL' in s:
        raise RuntimeError(f'Factorization failed: {path}')
    records={}
    for prefix,expected in [('TAU ',10000),('SELECTED ',16),('SUMMARY ',1),('FACTOR_R0 ',1),('PAIR_R1 ',2)]:
        rows=sorted(x for x in s.splitlines() if x.startswith(prefix))
        if len(rows)!=expected:
            raise RuntimeError(f'{path}: {prefix} expected {expected}, got {len(rows)}')
        for i,row in enumerate(rows):
            if prefix in ('FACTOR_R0 ','PAIR_R1 ') and 'status=PASS' not in row:
                raise RuntimeError('Missing factorization PASS')
            # Compare physical results; separately enforce factorization status.
            row=re.sub(r'\s+(?:absdiff|reldiff)=\S+','',row)
            records[(prefix,i)]=row
    pairs='\n'.join(records[k] for k in records if k[0]=='PAIR_R1 ')
    if 'kind=cross_boundary' not in pairs or 'y3_halo=1' not in pairs or 'y4_halo=1' not in pairs:
        raise RuntimeError('Missing MPI halo pair coverage')
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
        print(name.upper(),records[('SUMMARY ',0)])
        for k in records:
            if k[0]=='PAIR_R1 ': print(name.upper(),records[k])
    print('PASS qX temporal-HYP one-link observable integration')
if __name__=='__main__': main()
