"""Validate structure and export paired raw profiles, without a piter4 bias gate."""
import csv
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
import common as c
import point_compare as point
import average_compare as avg
import profile_compare as profile

FIELDS = avg.FIELDS
RAW_COLUMNS = ['configuration','piter','Nv','tau','t0','t1','tmid','axis',
               'R','long','tr1','tr2','long_from_midpoint','Nsrc',*FIELDS]
PAIRED_COLUMNS = ['configuration','tau','axis','R','long','tr1','tr2',
                  *[f'{k}_p{p}' for p in (4,40) for k in FIELDS],
                  'C_p4','C_p40','C_difference','C_relative_difference']
SU3_COLUMNS = ['configuration','tau','piter','t','t0','t1','count',
               'unit_max','unit_rms','det_max','det_rms',
               'unit_above_1em2','det_above_1em2']

def write_csv(path, columns, data):
    with Path(path).open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns)
        writer.writeheader()
        writer.writerows(data)

def only(rows, label):
    if len(rows) != 1:
        raise RuntimeError('Expected one '+label)
    return rows[0]

def numeric(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise RuntimeError('Nonfinite '+key)
    return value

def su3(text,cfg,tau,piter):
    t0,t1 = 98-tau//2,98+tau//2
    times = avg.rows(text,'SU3_TIME')
    window = only(avg.rows(text,'SU3_WINDOW'),'SU3_WINDOW')
    if len(times)!=tau or {int(r['t']) for r in times}!=set(range(t0,t1)):
        raise RuntimeError('Wrong SU3 time coverage')
    if (int(window['t0']),int(window['t1']))!=(t0,t1):
        raise RuntimeError('Wrong SU3 window')
    for r in times+[window]:
        count = int(r['count'])
        if (int(r['cfg']),int(r['piter']))!=(cfg,piter) or count!=(64**3*(tau if r is window else 1)):
            raise RuntimeError('Wrong SU3 metadata/count')
        for key in ('unit','det'):
            maximum,rms = numeric(r,key+'_max'),numeric(r,key+'_rms')
            if not 0 <= rms <= maximum*(1+1e-12)+1e-20:
                raise RuntimeError('Invalid SU3 max/RMS')
            tail = int(r[key+'_above_1em2'])
            if not 0<=tail<=count:
                raise RuntimeError('Invalid SU3 tail count')
            if (maximum <= 1e-2 and tail != 0) or (maximum > 1e-2 and tail == 0):
                raise RuntimeError('Inconsistent SU3 tail count')
    for key in ('unit','det'):
        avg.near(numeric(window,key+'_max'),max(numeric(r,key+'_max') for r in times),'window maximum')
        expected = math.sqrt(math.fsum(numeric(r,key+'_rms')**2 for r in times)/tau)
        avg.near(numeric(window,key+'_rms'),expected,'window RMS')
        if int(window[key+'_above_1em2'])!=sum(int(r[key+'_above_1em2']) for r in times):
            raise RuntimeError('Window tail sum differs')
    if piter==40 and max(numeric(window,'unit_max'),numeric(window,'det_max'))>1e-12:
        raise RuntimeError('piter40 control fails its 1e-12 SU3 tolerance')
    def export(r):
        return dict(configuration=cfg,tau=tau,piter=piter,t=r.get('t',''),t0=t0,t1=t1,
                    **{k:r[k] for k in SU3_COLUMNS[6:]})
    return [export(r) for r in sorted(times,key=lambda r:int(r['t']))],export(window)

def read_measurement(path,cfg,tau,piter):
    text = Path(path).read_text()
    if text.count(f'Nv=100 cfg={cfg} temporal_bc=open\n')!=1 or text.count('P=1 x 4 x 4 x 4\n')!=1:
        raise RuntimeError('Wrong configuration or observable process grid')
    support = only(avg.rows(text,'MPI_SUPPORT'),'MPI_SUPPORT')
    if int(support['cfg'])!=cfg or int(support['funneled'])<1 or int(support['threads'])!=1:
        raise RuntimeError('Wrong MPI/OpenMP initialization')
    hyp = only(avg.rows(text,'HYP_INTEGRATION'),'HYP_INTEGRATION')
    if (int(hyp['mode']),int(hyp['steps']),int(hyp['piter']))!=(2,1,piter):
        raise RuntimeError('Wrong temporal HYP prescription')
    if numeric(hyp,'preserved_change')>1e-13 or numeric(hyp,'temporal_change')<=1e-13:
        raise RuntimeError('Temporal HYP field isolation failed')
    points = point.read(Path(path),98-tau//2,98+tau//2)
    avg.read_average(Path(path))
    records = profile.read_scan(Path(path))
    times,window = su3(text,cfg,tau,piter)
    return points,records,times,window

def raw_rows(cfg,tau,piter,records):
    result = []
    for (axis,r,l,t,u),row in sorted(records.items()):
        result.append(dict(configuration=cfg,piter=piter,Nv=100,tau=tau,
                           t0=98-tau//2,t1=98+tau//2,tmid=98,axis=axis,R=r,
                           long=l,tr1=t,tr2=u,long_from_midpoint=l-r/2,
                           Nsrc=row['Nsrc'],**{k:row[k] for k in FIELDS}))
    return result

def connected(row):
    with localcontext() as context:
        context.prec=50
        l=Decimal(row['L_RE'])
        if not l:
            return None
        return Decimal(row['LS_RE'])/l-Decimal(row['S_RAW'])

def compare_pair(outdir,cfg,tau,prefix=''):
    data = {p:read_measurement(outdir/f'{prefix}p{p}.out',cfg,tau,p) for p in (4,40)}
    probes = [{k:v for k,v in data[p][0].items() if k[0]=='PROBE '} for p in (4,40)]
    point.compare(probes[0],probes[1],'unchanged original-field clover')
    scans = {p:data[p][1] for p in (4,40)}
    if scans[4].keys()!=scans[40].keys():
        raise RuntimeError('Paired profile keys differ')
    paired=[]
    max_relative = {'L_RE':0.0,'LS_RE':0.0,'C_all_axes':0.0,'C_combined':0.0}
    undefined=0
    for key,row4 in sorted(scans[4].items()):
        row40=scans[40][key]
        avg.near(float(row4['S_RAW']),float(row40['S_RAW']),'original clover across prescriptions')
        axis,r,l,t,u=key
        row=dict(configuration=cfg,tau=tau,axis=axis,R=r,long=l,tr1=t,tr2=u)
        for p,source in ((4,row4),(40,row40)):
            row.update({f'{k}_p{p}':source[k] for k in FIELDS})
        a,b=connected(row4),connected(row40)
        with localcontext() as context:
            context.prec=50
            difference=None if a is None or b is None else a-b
            relative=None if difference is None or b==0 else abs(difference)/abs(b)
        row.update(C_p4='' if a is None else str(a),C_p40='' if b is None else str(b),
                   C_difference='' if difference is None else str(difference),
                   C_relative_difference='' if relative is None else str(relative))
        if relative is None: undefined+=1
        else:
            max_relative['C_all_axes']=max(max_relative['C_all_axes'],float(relative))
            if axis==-1:max_relative['C_combined']=max(max_relative['C_combined'],float(relative))
        for name in ('L_RE','LS_RE'):
            a0,b0=float(row4[name]),float(row40[name])
            if b0:max_relative[name]=max(max_relative[name],abs(a0-b0)/abs(b0))
        paired.append(row)
    raw=[x for p in (4,40) for x in raw_rows(cfg,tau,p,scans[p])]
    times=[x for p in (4,40) for x in data[p][2]]
    windows=[data[p][3] for p in (4,40)]
    write_csv(outdir/(prefix+'profiles_raw.csv'),RAW_COLUMNS,raw)
    write_csv(outdir/(prefix+'profiles_paired.csv'),PAIRED_COLUMNS,paired)
    write_csv(outdir/(prefix+'su3_times.csv'),SU3_COLUMNS,times)
    write_csv(outdir/(prefix+'su3_windows.csv'),SU3_COLUMNS,windows)
    return dict(configuration=cfg,tau=tau,raw_records=len(raw),paired_records=len(paired),
                max_relative_difference=max_relative,undefined_connected_ratios=undefined,
                p4_window=windows[0],p40_window=windows[1]),scans

def regression_check(scans):
    with (c.CASE/'n100_tau4_reference.csv').open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    expected={}
    for row in rows:
        if (row['configuration'],int(row['Nv']),int(row['tau']),int(row['t0']),int(row['t1']))!=('qX_p2n100',100,4,96,100):
            raise RuntimeError('Wrong n100 reference metadata')
        p=int(row['piter'])
        key=tuple(int(row[k]) for k in ('axis','R','long','tr1','tr2'))
        if (p,key) in expected:raise RuntimeError('Duplicate n100 reference key')
        if int(row['Nsrc'])!=(3*64**3 if key[0]==-1 else 64**3):
            raise RuntimeError('Wrong n100 reference source count')
        expected[(p,key)]=row
    wanted={(p,key) for p in (4,40) for key in scans[p]}
    if expected.keys()!=wanted:raise RuntimeError('Incomplete n100 reference')
    worst=0.0
    for p,key in sorted(wanted):
        for name in FIELDS:
            a,b=float(scans[p][key][name]),float(expected[(p,key)][name])
            avg.near(a,b,'n100 tau4 regression '+name)
            worst=max(worst,abs(a-b))
    print('PASS n100 tau4 regression against job5915129; max difference',worst)
    return dict(success=True,reference_job=5915129,records=664,max_absolute_difference=worst)

def check(outdir,cfg,tau,regression=False):
    summary,scans=compare_pair(outdir,cfg,tau)
    summary.update(success=True,scope='Structural/profile checks and piter40 SU3 control; '
                   'piter4 defects and observable sensitivity are reported without an acceptance gate.')
    if regression:
        ignored,ref_scans=compare_pair(outdir,100,4,'n100_')
        summary['n100_regression']=regression_check(ref_scans)
    return summary
