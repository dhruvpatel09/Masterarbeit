#!/usr/bin/env python3
"""Validate scan structure and original probes; report piter4/40 differences."""
import csv
import math
import sys
from pathlib import Path
import average_compare as avg
import point_compare as point
import profile_compare as profile

FIELDS=avg.FIELDS
TAUS=(2,4,6)


def checked(path,tau,piter):
    text=path.read_text()
    tags=avg.rows(text,'HYP_INTEGRATION')
    if len(tags)!=1 or any(int(tags[0][k])!=v for k,v in
                           [('mode',2),('steps',1),('piter',piter)]):
        raise RuntimeError('Wrong HYP settings: '+str(path))
    if text.count('Nv=100 cfg=100 temporal_bc=open')!=1:
        raise RuntimeError('Wrong configuration or Nv: '+str(path))
    points=point.read(path,98-tau//2,98+tau//2)
    averages=avg.read_average(path)
    records=profile.read_scan(path)
    return points,averages,records


def compare_control(reference,control):
    point.compare(reference[0],control[0],'fresh piter40 point regression')
    point.compare(reference[1],control[1],'fresh piter40 average regression')
    if reference[2].keys()!=control[2].keys():
        raise RuntimeError('Fresh control profile keys differ')
    for key in reference[2]:
        for field in FIELDS:
            avg.near(float(reference[2][key][field]),float(control[2][key][field]),
                     'fresh piter40 profile regression '+field)
    print('PASS fresh piter40 profile regression')


def compare_probes(reference,measurement):
    ref_probes={k:v for k,v in reference[0].items() if k[0]=='PROBE '}
    new_probes={k:v for k,v in measurement[0].items() if k[0]=='PROBE '}
    point.compare(ref_probes,new_probes,'original clover probes piter4 vs 40')
    if reference[2].keys()!=measurement[2].keys():
        raise RuntimeError('Profile keys differ')
    for key in reference[2]:
        avg.near(float(reference[2][key]['S_RAW']),float(measurement[2][key]['S_RAW']),
                 'original clover average piter4 vs 40')
    print('PASS unchanged original clover average across all profile records')


def fractional_change(new,reference):
    """Signed change divided by |reference|; undefined if reference is zero."""
    return (new-reference)/abs(reference) if reference!=0 else None


def export(root,ref,measured):
    raw=root/'profile_raw_qX_p2n100_mid98_piter4_40.csv'
    columns=['configuration','mode','piter','Nv','tau','t0','t1','tmid','axis',
             'R','long','tr1','tr2','long_from_midpoint','Nsrc',*FIELDS]
    with raw.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        for tau in TAUS:
            for label,piter,collection in [('hyp4',4,measured),('hyp40_reference',40,ref)]:
                for (axis,r,l,t,u),row in sorted(collection[tau][2].items()):
                    out=dict(configuration='qX_p2n100',mode=label,piter=piter,Nv=100,
                             tau=tau,t0=98-tau//2,t1=98+tau//2,tmid=98,axis=axis,
                             R=r,long=l,tr1=t,tr2=u,long_from_midpoint=l-r/2,
                             Nsrc=row['Nsrc'])
                    out.update({k:row[k] for k in FIELDS});writer.writerow(out)
    difference=root/'profile_piter4_minus40_qX_p2n100_mid98.csv'
    cols=['configuration','tau','axis','R','long','tr1','tr2','Nsrc']
    cols += [f'{field}_{suffix}' for field in FIELDS for suffix in ('p4','p40','difference','fractional_change')]
    with difference.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=cols);writer.writeheader()
        for tau in TAUS:
            for (axis,r,l,t,u),row in sorted(measured[tau][2].items()):
                other=ref[tau][2][(axis,r,l,t,u)]
                out=dict(configuration='qX_p2n100',tau=tau,axis=axis,R=r,
                         long=l,tr1=t,tr2=u,Nsrc=row['Nsrc'])
                for field in FIELDS:
                    a,b=float(row[field]),float(other[field])
                    fraction=fractional_change(a,b)
                    delta=a-b
                    if not math.isfinite(delta) or (fraction is not None and not math.isfinite(fraction)):
                        raise RuntimeError('Nonfinite comparison value')
                    out.update({field+'_p4':a,field+'_p40':b,
                                field+'_difference':delta,
                                field+'_fractional_change':fraction})
                writer.writerow(out)
            for r in (1,2,4):
                key=(-1,r,0,1,2)
                a,b=measured[tau][2][key],ref[tau][2][key]
                print('PITER_OBSERVABLE',f'tau={tau} R={r}',
                      ' '.join(f'{f}_P4={a[f]} {f}_P40={b[f]}' for f in FIELDS))
    print('CSV:',raw,'(1992 raw single-configuration records)')
    print('CSV:',difference,'(996 signed-difference records; fractional_change uses abs(p40))')


def main():
    root,reference=map(Path,sys.argv[1:])
    ref={tau:checked(reference/f'tau{tau}.out',tau,40) for tau in TAUS}
    measured={tau:checked(root/f'p4_tau{tau}.out',tau,4) for tau in TAUS}
    control=checked(root/'p40_tau4_control.out',4,40)
    compare_control(ref[4],control)
    for tau in TAUS:
        compare_probes(ref[tau],measured[tau])
    export(root,ref,measured)
    print('PASS scan structure, original probes and fresh piter40 regression; '
          'piter4-versus40 observable differences reported without an acceptance threshold')


if __name__=='__main__':main()
