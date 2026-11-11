"""Replay physical endpoints, complete matrices, native audits and reductions."""
import csv,math
from pathlib import Path
from parents import legacy,design
from pilot_plan import PILOT,CUTOFFS,RADII,TAUS,AXES,expected_keys,direct_keys,legacy_keys,block_keys

rows=legacy.rows
near=legacy.near


def fields(row,names):
    if set(row)!=set(names.split()):raise ValueError('Unexpected diagnostic fields')


def load_references(path):
    with Path(path).open(newline='') as stream:values=list(csv.DictReader(stream))
    reference={}
    for row in values:
        key=tuple(int(row[k]) for k in ('configuration','Nv_src','Nv_sink','tau','R'))
        cfg,s,k,t,r=key
        if key in reference or row['t0']!='90' or int(row['t1'])!=90+t or row['piter']!='10' or row['axis']!='-1' or int(row['Nsrc'])!=3*64**3:
            raise ValueError('Wrong baseline reference metadata')
        if not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite reference')
        reference[key]=row
    expected={(c,s,k,t,r) for c in PILOT for s in CUTOFFS for k in CUTOFFS for t in TAUS for r in RADII}
    if len(values)!=6336 or set(reference)!=expected:raise ValueError('Incomplete published matrix reference')
    return reference


def audit(path,cfg):
    text=Path(path).read_text();times=design.audit_times(cfg)
    start=rows(text,'QUALITY_START');end=rows(text,'QUALITY_END');modes=rows(text,'MODE');slices=rows(text,'SLICE')
    if len(start)!=1 or len(end)!=1 or end[0]!=dict(cfg=str(cfg),slices=str(len(times)),status='PASS'):raise ValueError('Incomplete selected-time audit')
    s=start[0];fields(s,'cfg grid Nv APE FUNNELED smear_seconds')
    if s['cfg']!=str(cfg) or s['grid']!='1x4x4x4' or s['Nv']!='100' or s['APE']!='20,0.5' or int(s['FUNNELED'])<1 or not math.isfinite(float(s['smear_seconds'])) or float(s['smear_seconds'])<0:
        raise ValueError('Wrong saved-mode audit setup')
    expected={(t,i) for t in times for i in range(100)};lookup={}
    for row in modes:
        fields(row,'cfg t i lambda norm residual_l2 residual_relative rayleigh_re rayleigh_im')
        key=(int(row['t']),int(row['i']))
        if key in lookup or row['cfg']!=str(cfg):raise ValueError('Duplicate/wrong mode audit identity')
        values=[float(row[k]) for k in ('lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im')]
        lam,norm,residual,relative,rayre,rayim=values
        if not all(math.isfinite(v) for v in values) or not 0<lam<=12.00000001 or norm<=0 or residual<0 or not 0<=relative<=1e-9 or abs(norm-1)>1e-9 or abs(rayim)>1e-10:
            raise ValueError('Saved mode fails unchanged residual/norm criteria')
        near(relative,residual/(lam*norm),rtol=5e-12,atol=1e-25)
        if abs(rayre-lam)>residual/norm+5e-14*max(lam,1.):raise ValueError('Rayleigh value inconsistent with mode residual')
        lookup[key]=row
    if len(modes)!=len(expected) or set(lookup)!=expected:raise ValueError('Wrong saved-mode time/mode inventory')
    if len(slices)!=len(times) or {int(r['t']) for r in slices}!=set(times):raise ValueError('Wrong saved-mode slice inventory')
    for row in slices:
        fields(row,'cfg t Nv max_relative_residual max_norm_defect gram_diag_defect gram_offdiag max_rayleigh_im lambda_min lambda_max status')
        if row['cfg']!=str(cfg) or row['Nv']!='100' or row['status']!='PASS':raise ValueError('Failed saved-mode slice')
        t=int(row['t']);selected=[lookup[t,i] for i in range(100)];lambdas=[float(r['lambda']) for r in selected]
        if any(lambdas[i+1]<lambdas[i]-1e-12 for i in range(99)):raise ValueError('Saved modes not ordered by eigenvalue')
        for key,limit in [('max_relative_residual',1e-9),('max_norm_defect',1e-9),('gram_diag_defect',1e-9),('gram_offdiag',1e-9),('max_rayleigh_im',1e-10)]:
            v=float(row[key])
            if not math.isfinite(v) or not 0<=v<=limit:raise ValueError('Invalid saved-mode slice quality')
        near(float(row['max_relative_residual']),max(float(r['residual_relative']) for r in selected),atol=1e-25)
        near(float(row['max_norm_defect']),max(abs(float(r['norm'])-1) for r in selected),atol=5e-16)
        near(float(row['max_rayleigh_im']),max(abs(float(r['rayleigh_im'])) for r in selected),atol=1e-25)
        near(float(row['lambda_min']),lambdas[0]);near(float(row['lambda_max']),lambdas[-1])
    return dict(configuration=cfg,slices=len(times),modes=len(modes),maximum_relative_residual=max(float(r['residual_relative']) for r in modes),saved_eigenvalue_order_checked=True)


def projection(text,cfg,piter,t0):
    times=rows(text,'SU3_TIME');window=rows(text,'SU3_WINDOW')
    if len(times)!=32 or {int(r['t']) for r in times}!=set(range(t0,t0+32)) or len(window)!=1:raise ValueError('Wrong shifted SU3 time inventory')
    if (window[0]['t0'],window[0]['t1'])!=(str(t0),str(t0+32)):raise ValueError('Wrong shifted transporter window')
    for row in times+window:
        required='cfg piter count unit_max unit_rms det_max det_rms unit_above_1em2 det_above_1em2 '+('t0 t1' if row is window[0] else 't')
        fields(row,required)
        if row['cfg']!=str(cfg) or row['piter']!=str(piter) or int(row['count'])!=64**3*(32 if row is window[0] else 1):raise ValueError('Wrong shifted SU3 count')
        for prefix in ('unit','det'):
            maximum,rms=float(row[prefix+'_max']),float(row[prefix+'_rms'])
            if not math.isfinite(maximum) or not math.isfinite(rms) or not 0<=rms<=maximum<=1e-12 or row[prefix+'_above_1em2']!='0':raise ValueError('Transporter projection tolerance failed')
    for prefix in ('unit','det'):
        near(float(window[0][prefix+'_max']),max(float(r[prefix+'_max']) for r in times),atol=1e-30)
        near(float(window[0][prefix+'_rms']),math.sqrt(math.fsum(float(r[prefix+'_rms'])**2 for r in times)/32),atol=1e-30)
    return window[0]


def scan(path,cfg,spec,reference):
    piter,t0=spec['piter'],spec['t0']
    if spec not in design.scans(cfg):raise ValueError('Unplanned native scan')
    text=Path(path).read_text();begin=rows(text,'MATRIX_BEGIN');end=rows(text,'MATRIX_END')
    expected=expected_keys();header=dict(protocol='source_time_v1',cfg=str(cfg),piter=str(piter),saved_Nv='100',
        source_cutoffs='40,60,80,100',sink_cutoffs='40,60,80,100',ranks='64',grid='1x4x4x4',t0=str(t0),t1_max=str(t0+32),
        tau_min='0',tau_max='32',radii='4,12,24,32',temporal_bc='open',hyp_steps='1')
    ending=dict(protocol='source_time_v1',cfg=str(cfg),piter=str(piter),t0=str(t0),t1_max=str(t0+32),records=str(len(expected)),status='PASS')
    if len(begin)!=1 or len(end)!=1 or end[0]!=ending:raise ValueError('Incomplete source-time matrix scan')
    fields(begin[0],' '.join(header)+' temporal_change preserved_change')
    if any(begin[0][k]!=v for k,v in header.items()):raise ValueError('Wrong source-time scan geometry')
    changed,preserved=float(begin[0]['temporal_change']),float(begin[0]['preserved_change'])
    if not math.isfinite(changed) or changed<=1e-13 or not math.isfinite(preserved) or not 0<=preserved<=1e-13:raise ValueError('Temporal HYP isolation failed')
    loaded=rows(text,'MODE_LOAD');seen=set()
    for row in loaded:
        fields(row,'cfg piter t0 loaded_t role status')
        key=(int(row['loaded_t']),row['role'])
        if key in seen or row['cfg']!=str(cfg) or row['piter']!=str(piter) or row['t0']!=str(t0) or row['status']!='PASS':raise ValueError('Wrong HDF5 endpoint load identity')
        seen.add(key)
    if len(loaded)!=33 or seen!={(t0,'source')}|{(t0+t,'sink') for t in range(1,33)}:raise ValueError('Missing/duplicate endpoint HDF5 loads')
    values=rows(text,'MATRIX_RT');lookup={}
    for row in values:
        fields(row,'cfg piter scheme t0 t1 tau R axis Nv_src Nv_sink Nsrc L_RE L_IM')
        key=tuple(int(row[k]) for k in ('Nv_src','Nv_sink','tau','R','axis'))
        if key in lookup or row['cfg']!=str(cfg) or row['piter']!=str(piter) or row['scheme']!='fixed' or row['t0']!=str(t0) or int(row['t1'])!=t0+key[2]:raise ValueError('Wrong physical matrix endpoint/identity')
        if int(row['Nsrc'])!=64**3*(3 if key[-1]==-1 else 1) or not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Wrong source count/nonfinite matrix')
        lookup[key]=row
    if len(values)!=len(expected) or set(lookup)!=expected:raise ValueError('Incomplete ordered source/sink matrix')
    matches=0;maximum=0.
    for (s,k,t,r,a),row in lookup.items():
        if a==-1:
            for field in ('L_RE','L_IM'):near(float(row[field]),math.fsum(float(lookup[s,k,t,r,q][field]) for q in (0,1,2))/3)
            if spec['purpose']=='fresh_reference_regression':
                ref=reference[cfg,s,k,t,r];matches+=1
                for field in ('L_RE','L_IM'):
                    x,y=float(row[field]),float(ref[field]);near(x,y);maximum=max(maximum,abs(x-y))
        if r==32 and abs(float(row['L_IM']))>1e-18+1e-10*abs(float(row['L_RE'])):raise ValueError('Half-box imaginary cancellation failed')
    if matches!=(2112 if spec['purpose']=='fresh_reference_regression' else 0):raise ValueError('Fresh reference match inventory')
    gram=[]
    for r in RADII:
        for a in AXES:
            matrix=[[complex(float(lookup[s,k,0,r,a]['L_RE']),float(lookup[s,k,0,r,a]['L_IM'])) for k in CUTOFFS] for s in CUTOFFS]
            gram.append(legacy.equal_time_matrix(matrix))
    for prefix,expected_direct in (('MATRIX_DIRECT',direct_keys()),('MATRIX_LEGACY',legacy_keys())):
        for row in rows(text,prefix):
            fields(row,'cfg piter t0 t1 tau R axis rank site Nv_src Nv_sink fast_re fast_im direct_re direct_im abs_error tolerance status')
            if row['t0']!=str(t0) or int(row['t1'])!=t0+int(row['tau']):raise ValueError('Wrong independent-test physical endpoints')
        legacy.direct_checks(text,prefix,expected_direct,cfg,piter)
    roundtrip=rows(text,'PERIODIC_ROUNDTRIP')
    if len(roundtrip)!=3 or {int(r['axis']) for r in roundtrip}!={0,1,2} or any(r['cells']!='64' or r['status']!='PASS' for r in roundtrip):raise ValueError('Incomplete periodic payload roundtrips')
    central=projection(text,cfg,piter,t0);blocks=rows(text,'MATRIX_BLOCK');grouped={};positions={}
    for row in blocks:
        fields(row,'cfg piter t0 t1 tau R axis Nv_src Nv_sink rank px py pz Nsrc L_RE L_IM')
        key=tuple(int(row[k]) for k in ('Nv_src','Nv_sink','tau','R','axis'));rank=int(row['rank']);position=tuple(int(row[k]) for k in ('px','py','pz'))
        if row['cfg']!=str(cfg) or row['piter']!=str(piter) or row['t0']!=str(t0) or int(row['t1'])!=t0+key[2] or row['Nsrc']!='4096':raise ValueError('Wrong rank-block endpoint/identity')
        if rank in positions and positions[rank]!=position:raise ValueError('Rank-coordinate mapping changed')
        positions[rank]=position
        if not all(math.isfinite(float(row[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite rank block')
        grouped.setdefault(key,[]).append(row)
    if len(blocks)!=len(block_keys())*64 or set(grouped)!=block_keys():raise ValueError('Incomplete rank-block inventory')
    for key,group in grouped.items():
        if len(group)!=64 or {int(r['rank']) for r in group}!=set(range(64)) or {tuple(int(r[k]) for k in ('px','py','pz')) for r in group}!={(x,y,z) for x in range(4) for y in range(4) for z in range(4)}:raise ValueError('Duplicate/missing rank coordinates')
        for field in ('L_RE','L_IM'):near(math.fsum(float(r[field]) for r in group)/64,float(lookup[key][field]))
    return values,blocks,dict(configuration=cfg,label=spec['label'],piter=piter,t0=t0,rows=len(values),rank_block_rows=len(blocks),
        fresh_reference_matches=matches,reference_maximum_absolute_difference=maximum,direct_tests=len(direct_keys()),
        legacy_tests=len(legacy_keys()),equal_time_gram_checks=len(gram),periodic_roundtrips=3,endpoint_loads=len(loaded),
        maximum_equal_time_relative_hermitian_defect=max(r['relative_hermitian_defect'] for r in gram),
        minimum_equal_time_relative_pivot=min(r['minimum_relative_pivot'] for r in gram),
        unit_max=float(central['unit_max']),det_max=float(central['det_max']))


def projection_pair(low,high):
    def index(values):return {tuple(row[k] for k in ('cfg','scheme','t0','t1','tau','R','axis','Nv_src','Nv_sink')):row for row in values}
    a,b=index(low),index(high)
    if len(a)!=8448 or len(b)!=8448 or set(a)!=set(b):raise ValueError('Unmatched physical geometry in projection control')
    maximum=0.
    for key in a:
        for field in ('L_RE','L_IM'):
            x,y=float(a[key][field]),float(b[key][field]);near(x,y,rtol=1e-10);maximum=max(maximum,abs(x-y))
    return maximum
