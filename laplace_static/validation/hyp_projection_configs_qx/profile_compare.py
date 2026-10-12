#!/usr/bin/env python3
"""Validate endpoint shifts and raw profile averages; export per-config CSV."""
import csv
import math
from pathlib import Path
import sys
import average_compare as avg
import point_compare as point

FIELDS=avg.FIELDS

def offsets(r):
    return {(0,1,2)}|{(l,t,0) for l in range(-1,r+2) for t in range(-2,3)}

def read_scan(path):
    text=path.read_text()
    if text.count('PROFILE_SCAN status=PASS')!=1:
        raise RuntimeError('Missing profile PASS')
    checks=avg.rows(text,'SHIFT_CHECK')
    if len(checks)!=9 or {(int(x['axis']),int(x['R'])) for x in checks}!={(a,r) for a in range(3) for r in (1,2,4)}:
        raise RuntimeError('Missing/duplicate shift checks')
    for row in checks:
        r=int(row['R'])
        if (int(row['Nsrc']),int(row['Ncross']),int(row['Nwrap']),int(row['direct_samples']))!=(64**3,r*16384,r*4096,3) or row['status']!='PASS':
            raise RuntimeError('Shift counts/direct checks failed')
        if not math.isfinite(float(row['max_direct_difference'])):
            raise RuntimeError('Nonfinite direct check')
    records={}
    for prefix in ('SCAN_AXIS','SCAN_COMBINED'):
        for row in avg.rows(text,prefix):
            axis=int(row['axis']) if prefix=='SCAN_AXIS' else -1
            r=int(row['R']); location=tuple(int(row[k]) for k in ('long','tr1','tr2'))
            key=(axis,r,*location)
            if key in records: raise RuntimeError('Duplicate profile record')
            if int(row['Nsrc'])!=(64**3 if axis!=-1 else 3*64**3):
                raise RuntimeError('Wrong profile denominator')
            if not all(math.isfinite(float(row[k])) for k in FIELDS):
                raise RuntimeError('Nonfinite profile')
            records[key]=row
    wanted={(a,r,*p) for a in (-1,0,1,2) for r in (1,2,4) for p in offsets(r)}
    if records.keys()!=wanted: raise RuntimeError('Missing/unexpected profile offsets')
    vacuum=float(next(x for x in avg.rows(text,'PROBE') if x['kind']=='spatial_average')['S_RAW'])
    for (a,r,l,t,u),row in records.items():
        avg.near(float(row['S_RAW']),vacuum,'translated probe average')
        ref=records[(a,r,0,1,2)]
        for field in ('L_RE','L_IM'):
            avg.near(float(row[field]),float(ref[field]),'L independent of probe offset')
        if a==-1:
            for field in FIELDS:
                vals=[float(records[(b,r,l,t,u)][field]) for b in range(3)]
                avg.near(float(row[field]),math.fsum(vals)/3,'combined profile '+field,math.fsum(map(abs,vals))/3)
    old={int(x['axis']):x for x in avg.rows(text,'AXIS_AVG')}
    old[-1]=avg.rows(text,'COMBINED_AVG')[0]
    for a in (-1,0,1,2):
        for field in FIELDS:
            avg.near(float(records[(a,1,0,1,2)][field]),float(old[a][field]),'shifted R1 vs halo R1 '+field)
    print(f'PASS {path.stem}: shifts R=1,2,4; 332 profile records; counts; probe translations; R1 regression')
    return records

def main():
    root=Path(sys.argv[1]); reference=Path(sys.argv[2])
    all_records={}
    for name in ('baseline','hyp40'):
        newpath=root/(name+'.out'); refpath=reference/(name+'.out')
        point.compare(point.read(refpath),point.read(newpath),name+' point regression vs 5903813')
        point.compare(avg.read_average(refpath),avg.read_average(newpath),name+' averages regression vs 5903813')
        all_records[name]=read_scan(newpath)
    for key in all_records['baseline']:
        avg.near(float(all_records['baseline'][key]['S_RAW']),float(all_records['hyp40'][key]['S_RAW']),'original probe baseline vs HYP')
    target=root/'profile_raw_qX_p2n100_t96_100.csv'
    columns=['configuration','mode','piter','Nv','t0','t1','tmid','axis','R','long','tr1','tr2','long_from_midpoint','Nsrc',*FIELDS]
    with target.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns); writer.writeheader()
        for name in ('baseline','hyp40'):
            for key,row in sorted(all_records[name].items()):
                axis,r,l,t,u=key
                data=dict(configuration='qX_p2n100',mode=name,piter=0 if name=='baseline' else 40,Nv=100,
                          t0=96,t1=100,tmid=98,axis=axis,R=r,long=l,tr1=t,tr2=u,
                          long_from_midpoint=l-r/2,Nsrc=row['Nsrc'])
                data.update({k:row[k] for k in FIELDS});writer.writerow(data)
            for r in (1,2,4):
                row=all_records[name][(-1,r,0,1,2)]
                print(name.upper(),'R='+str(r),' '.join(k+'='+row[k] for k in FIELDS))
    print(f'CSV: {target} (664 rows; axis=-1 means combined; raw single-configuration averages)')
    print('PASS qX separation and probe-plane scan with temporal HYP and original clover')
if __name__=='__main__': main()
