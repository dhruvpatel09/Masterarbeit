"""Strict replay of all source/sink pairs, controls and spatial reductions."""
import csv
import importlib.util
import math
import sys
from pathlib import Path
from matrix_plan import CUTOFFS,RADII,TAUS,AXES,expected_keys,direct_keys,legacy_keys,block_keys

HERE=Path(__file__).resolve().parent
parent=HERE.parent/'qx_static_extension/check_extension.py'
previous=sys.path[:];sys.path.insert(0,str(parent.parent))
try:
    spec=importlib.util.spec_from_file_location('qx_matrix_saved_mode_checks',str(parent))
    legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
finally:sys.path[:]=previous
audit=legacy.audit
near=legacy.near


def rows(text,prefix):
    result=[]
    for line in text.splitlines():
        if not line.startswith(prefix+' '):continue
        row={}
        for token in line.split()[1:]:
            if '=' not in token:raise ValueError('Malformed '+prefix+' token')
            key,value=token.split('=',1)
            if key in row:raise ValueError('Duplicate '+prefix+' field')
            row[key]=value
        result.append(row)
    return result


def require(row,fields):
    if not set(fields.split()).issubset(row):raise ValueError('Missing diagnostic field')


def load_references(path,configs):
    result={};configs=set(configs)
    with Path(path).open(newline='') as f:
        for row in csv.DictReader(f):
            if row['scheme']!='fixed' or row['piter']!='10' or int(row['cfg']) not in configs:continue
            n,t,r,a=(int(row[k]) for k in ('Nv','tau','R','axis'))
            if n not in CUTOFFS or t not in TAUS or r not in RADII or a not in AXES:continue
            key=(int(row['cfg']),n,t,r,a)
            if key in result or row['t0']!='90' or int(row['Nsrc'])!=64**3*(3 if a==-1 else 1):raise ValueError('Invalid diagonal reference metadata')
            if not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite reference')
            result[key]=row
    expected={(c,n,t,r,a) for c in configs for n in CUTOFFS for t in TAUS for r in RADII for a in AXES}
    if set(result)!=expected:raise ValueError('Incomplete published diagonal reference')
    return result


def equal_time_matrix(matrix):
    scale=max(abs(matrix[i][i]) for i in range(4));tol=1e-20+1e-10*scale;defect=0.
    for i in range(4):
        for j in range(4):
            difference=abs(matrix[i][j]-matrix[j][i].conjugate());defect=max(defect,difference)
            if difference>tol:raise ValueError('Equal-time matrix is not Hermitian within roundoff')
    lower=[[0j]*4 for _ in range(4)];pivots=[]
    for j in range(4):
        value=matrix[j][j].real-math.fsum(abs(lower[j][k])**2 for k in range(j))
        if value < -tol:raise ValueError('Equal-time matrix is not positive semidefinite')
        pivots.append(value);lower[j][j]=math.sqrt(max(value,0.)) if value>tol else 0.
        for i in range(j+1,4):
            residual=matrix[i][j]-sum(lower[i][k]*lower[j][k].conjugate() for k in range(j))
            if lower[j][j]:lower[i][j]=residual/lower[j][j]
            elif abs(residual)>tol:raise ValueError('Equal-time zero-pivot consistency failed')
    return dict(relative_hermitian_defect=defect/scale if scale else 0.,minimum_relative_pivot=min(pivots)/scale if scale else 0.)


def direct_checks(text,prefix,expected,cfg,piter):
    found=rows(text,prefix);inventory=set()
    for row in found:
        require(row,'cfg piter rank Nv_src Nv_sink tau R axis site fast_re fast_im direct_re direct_im abs_error tolerance status')
        if int(row['cfg'])!=cfg or int(row['piter'])!=piter or row['rank']!='0':raise ValueError('Incorrect direct-test metadata')
        key=tuple(int(row[k]) for k in ('Nv_src','Nv_sink','tau','R','axis','site'))
        if key in inventory:raise ValueError('Duplicate direct check')
        inventory.add(key)
        values=[float(row[k]) for k in ('fast_re','fast_im','direct_re','direct_im','abs_error','tolerance')]
        if not all(math.isfinite(x) for x in values):raise ValueError('Nonfinite direct check')
        a,b,c,d,error,tol=values;actual=math.hypot(a-c,b-d);bound=1e-20+5e-12*max(math.hypot(a,b),math.hypot(c,d))
        near(error,actual,rtol=1e-12,atol=1e-30);near(tol,bound,rtol=1e-12,atol=1e-30)
        if row['status']!='PASS' or not 0<=error<=bound or tol<=0:raise ValueError('Failed direct/legacy contraction check')
    if len(found)!=len(expected) or inventory!=expected:raise ValueError('Missing direct/legacy check coverage')
    return len(found)


def projection(text,cfg,piter):
    times=rows(text,'SU3_TIME');window=rows(text,'SU3_WINDOW')
    if len(times)!=32 or {int(r['t']) for r in times}!=set(range(90,122)) or len(window)!=1:raise ValueError('Incomplete central transporter diagnostics')
    for row in times+window:
        require(row,'cfg piter count unit_max unit_rms det_max det_rms unit_above_1em2 det_above_1em2')
        count=64**3*(32 if row is window[0] else 1)
        if int(row['cfg'])!=cfg or int(row['piter'])!=piter or int(row['count'])!=count:raise ValueError('Invalid central-link count')
        for field in ('unit','det'):
            maximum,rms=float(row[field+'_max']),float(row[field+'_rms'])
            if not math.isfinite(maximum) or not math.isfinite(rms) or not 0<=rms<=maximum<=1e-12:raise ValueError('Central-link numerical tolerance failed')
            if row[field+'_above_1em2']!='0':raise ValueError('Inconsistent central-link counters')
    if (window[0]['t0'],window[0]['t1'])!=('90','122'):raise ValueError('Wrong transporter time window')
    for field in ('unit','det'):
        near(float(window[0][field+'_max']),max(float(r[field+'_max']) for r in times),atol=1e-30)
        near(float(window[0][field+'_rms']),math.sqrt(math.fsum(float(r[field+'_rms'])**2 for r in times)/32),atol=1e-30)
    return window[0]


def scan(path,cfg,piter,references):
    text=Path(path).read_text();begin=rows(text,'MATRIX_BEGIN');end=rows(text,'MATRIX_END');expected=expected_keys()
    if len(begin)!=1 or len(end)!=1 or end[0]!=dict(cfg=str(cfg),piter=str(piter),records=str(len(expected)),status='PASS'):raise ValueError('Incomplete matrix scan')
    header=dict(cfg=str(cfg),piter=str(piter),saved_Nv='100',source_cutoffs='40,60,80,100',sink_cutoffs='40,60,80,100',
        ranks='64',grid='1x4x4x4',t0='90',tau_min='0',tau_max='32',radii='4,12,24,32',temporal_bc='open',hyp_steps='1')
    if any(begin[0].get(k)!=v for k,v in header.items()):raise ValueError('Wrong matrix geometry or cutoff orientation')
    require(begin[0],'temporal_change preserved_change')
    change,preserved=float(begin[0]['temporal_change']),float(begin[0]['preserved_change'])
    if not math.isfinite(change) or change<=1e-13 or not math.isfinite(preserved) or not 0<=preserved<=1e-13:raise ValueError('Temporal HYP isolation failed')
    values=rows(text,'MATRIX_RT');lookup={};reference_max=0.;reference_count=0
    fields={'cfg','piter','scheme','t0','tau','R','axis','Nv_src','Nv_sink','Nsrc','L_RE','L_IM'}
    for row in values:
        if set(row)!=fields or row['scheme']!='fixed' or row['t0']!='90' or int(row['cfg'])!=cfg or int(row['piter'])!=piter:raise ValueError('Invalid matrix record metadata')
        key=tuple(int(row[k]) for k in ('Nv_src','Nv_sink','tau','R','axis'))
        if key in lookup:raise ValueError('Duplicate matrix record')
        if int(row['Nsrc'])!=64**3*(3 if key[-1]==-1 else 1):raise ValueError('Wrong source count')
        if not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite matrix value')
        lookup[key]=row
    if len(values)!=len(expected) or set(lookup)!=expected:raise ValueError('Incomplete 4x4 source/sink matrix coverage')
    for key,row in lookup.items():
        source,sink,t,r,axis=key
        if axis==-1:
            for field in ('L_RE','L_IM'):near(float(row[field]),math.fsum(float(lookup[source,sink,t,r,a][field]) for a in (0,1,2))/3)
        if source==sink:
            reference=references[cfg,source,t,r,axis];reference_count+=1
            for field in ('L_RE','L_IM'):
                x,y=float(row[field]),float(reference[field]);near(x,y,rtol=5e-12 if piter==10 else 1e-10);reference_max=max(reference_max,abs(x-y))
        if r==32 and abs(float(row['L_IM']))>1e-18+1e-10*abs(float(row['L_RE'])):raise ValueError('Half-box imaginary cancellation failed')
    gram=[]
    for r in RADII:
        for axis in AXES:
            matrix=[[complex(float(lookup[s,k,0,r,axis]['L_RE']),float(lookup[s,k,0,r,axis]['L_IM'])) for k in CUTOFFS] for s in CUTOFFS]
            gram.append(equal_time_matrix(matrix))
    direct=direct_checks(text,'MATRIX_DIRECT',direct_keys(),cfg,piter)
    legacy_count=direct_checks(text,'MATRIX_LEGACY',legacy_keys(),cfg,piter)
    roundtrip=rows(text,'PERIODIC_ROUNDTRIP')
    if len(roundtrip)!=3 or {int(r['axis']) for r in roundtrip}!={0,1,2} or any(r['cells']!='64' or r['status']!='PASS' for r in roundtrip):raise ValueError('Incomplete periodic payload roundtrips')
    central=projection(text,cfg,piter);blocks=rows(text,'MATRIX_BLOCK');grouped={};rank_positions={}
    for row in blocks:
        require(row,'cfg piter t0 Nsrc Nv_src Nv_sink tau R axis rank px py pz L_RE L_IM')
        if int(row['cfg'])!=cfg or int(row['piter'])!=piter or row['t0']!='90' or row['Nsrc']!='4096':raise ValueError('Wrong matrix block metadata')
        key=tuple(int(row[k]) for k in ('Nv_src','Nv_sink','tau','R','axis'));rank=int(row['rank']);pos=tuple(int(row[k]) for k in ('px','py','pz'))
        if rank in rank_positions and rank_positions[rank]!=pos:raise ValueError('Rank-coordinate mapping changed')
        rank_positions[rank]=pos
        if not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite block reduction')
        grouped.setdefault(key,[]).append(row)
    if set(grouped)!=block_keys() or len(blocks)!=len(block_keys())*64:raise ValueError('Incomplete matrix rank-block inventory')
    for key,selected in grouped.items():
        if len(selected)!=64 or {int(r['rank']) for r in selected}!=set(range(64)) or {tuple(int(r[k]) for k in ('px','py','pz')) for r in selected}!={(x,y,z) for x in range(4) for y in range(4) for z in range(4)}:raise ValueError('Duplicate or missing rank coordinates')
        for field in ('L_RE','L_IM'):near(math.fsum(float(r[field]) for r in selected)/64,float(lookup[key][field]))
    return values,blocks,dict(configuration=cfg,piter=piter,rows=len(values),diagonal_reference_rows=reference_count,
        diagonal_reference_maximum_absolute_difference=reference_max,direct_tests=direct,legacy_tests=legacy_count,
        periodic_roundtrips=3,rank_block_rows=len(blocks),equal_time_gram_checks=len(gram),
        maximum_equal_time_relative_hermitian_defect=max(r['relative_hermitian_defect'] for r in gram),
        minimum_equal_time_relative_pivot=min(r['minimum_relative_pivot'] for r in gram),unit_max=float(central['unit_max']),det_max=float(central['det_max']))


def projection_pair(low,high):
    def index(values):return {tuple(r[k] for k in ('cfg','tau','R','axis','Nv_src','Nv_sink')):r for r in values}
    a,b=index(low),index(high)
    if not a or set(a)!=set(b):raise ValueError('Projection matrix controls do not match')
    maximum=0.
    for key in a:
        for field in ('L_RE','L_IM'):
            x,y=float(a[key][field]),float(b[key][field]);near(x,y,rtol=1e-10);maximum=max(maximum,abs(x-y))
    return maximum
