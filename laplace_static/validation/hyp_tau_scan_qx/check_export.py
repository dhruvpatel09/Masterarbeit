#!/usr/bin/env python3
import csv
import sys
from pathlib import Path
import point_compare as point
import average_compare as avg
import profile_compare as profile

def checked(path,tau):
    text=path.read_text()
    tags=avg.rows(text,'HYP_INTEGRATION')
    if len(tags)!=1 or any(int(tags[0][k])!=v for k,v in [('mode',2),('steps',1),('piter',40)]):
        raise RuntimeError('Wrong HYP settings')
    if text.count('Nv=100 cfg=100 temporal_bc=open')!=1:
        raise RuntimeError('Wrong configuration/Nv')
    points=point.read(path,98-tau//2,98+tau//2)
    averages=avg.read_average(path)
    records=profile.read_scan(path)
    return points,averages,records

def main():
    root,reference=map(Path,sys.argv[1:])
    ref=checked(reference,4)
    results={}
    for tau in (2,4,6):
        value=checked(root/f'tau{tau}.out',tau)
        results[tau]=value[2]
        if tau==4:
            point.compare(ref[0],value[0],'tau4 point regression vs 5903817')
            point.compare(ref[1],value[1],'tau4 averages regression vs 5903817')
            for key,row in value[2].items():
                for field in avg.FIELDS:
                    avg.near(float(row[field]),float(ref[2][key][field]),'tau4 profile regression')
        for key,row in value[2].items():
            avg.near(float(row['S_RAW']),float(ref[2][key]['S_RAW']),'fixed-midpoint original probe')
    target=root/'profile_raw_qX_p2n100_mid98_tau2_4_6.csv'
    cols=['configuration','mode','piter','Nv','tau','t0','t1','tmid','axis','R','long','tr1','tr2','long_from_midpoint','Nsrc',*avg.FIELDS]
    with target.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for tau,records in results.items():
            for (a,r,l,t,u),row in sorted(records.items()):
                out=dict(configuration='qX_p2n100',mode='hyp40',piter=40,Nv=100,tau=tau,t0=98-tau//2,t1=98+tau//2,tmid=98,axis=a,R=r,long=l,tr1=t,tr2=u,long_from_midpoint=l-r/2,Nsrc=row['Nsrc'])
                out.update({k:row[k] for k in avg.FIELDS});w.writerow(out)
    print('CSV:',target,'(996 raw single-configuration records)')
    print('PASS fixed-midpoint temporal scan tau=2,4,6')
if __name__=='__main__':main()
