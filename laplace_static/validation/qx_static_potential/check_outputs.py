#!/usr/bin/env python3
"""Validate the extended static scan, original-geometry controls and saved modes."""
import csv
import math
from pathlib import Path


def rows(text,prefix):
    return [dict(token.split('=',1) for token in line.split()[1:])
            for line in text.splitlines() if line.startswith(prefix+' ')]


def near(a,b,rtol=5e-12,atol=1e-20):
    if not math.isfinite(a) or not math.isfinite(b) or abs(a-b)>atol+rtol*max(abs(a),abs(b)):
        raise ValueError('Numerical comparison failed: %.17e versus %.17e'%(a,b))


def audit(path,cfg):
    text=Path(path).read_text();modes=rows(text,'MODE');slices=rows(text,'SLICE');end=rows(text,'QUALITY_END')
    if len(end)!=1 or end[0]!=dict(cfg=str(cfg),slices='17',status='PASS'):raise ValueError('Incomplete saved-mode audit')
    if len(modes)!=1700 or {(int(r['t']),int(r['i'])) for r in modes}!={(t,i) for t in range(90,107) for i in range(100)}:
        raise ValueError('Wrong saved-mode audit coverage')
    for r in modes:
        values=[float(r[k]) for k in ('lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im')]
        if int(r['cfg'])!=cfg or not all(math.isfinite(x) for x in values) or not 0<values[0]<=12.00000001:
            raise ValueError('Invalid mode audit value')
        if values[2]<0 or values[3]<0 or values[3]>1e-9 or abs(values[1]-1)>1e-9 or abs(values[5])>1e-10:
            raise ValueError('Saved mode fails unchanged residual/norm criteria')
    if len(slices)!=17 or {int(r['t']) for r in slices}!=set(range(90,107)):raise ValueError('Missing slice diagnostics')
    for r in slices:
        if int(r['cfg'])!=cfg or r['status']!='PASS' or int(r['Nv'])!=100:raise ValueError('Failed slice')
        for field,limit in [('max_relative_residual',1e-9),('max_norm_defect',1e-9),('gram_diag_defect',1e-9),('gram_offdiag',1e-9),('max_rayleigh_im',1e-10)]:
            value=float(r[field])
            if not math.isfinite(value) or not 0<=value<=limit:raise ValueError('Invalid slice quality')
    return dict(configuration=cfg,slices=17,modes=1700,maximum_relative_residual=max(float(r['residual_relative']) for r in modes))


def scan(path,cfg,piter,references):
    text=Path(path).read_text();begin=rows(text,'STATIC_BEGIN');end=rows(text,'STATIC_END')
    if len(begin)!=1 or len(end)!=1 or end[0]!=dict(cfg=str(cfg),piter=str(piter),records='2212',status='PASS'):
        raise ValueError('Incomplete static scan')
    b=begin[0]
    expected=dict(cfg=str(cfg),piter=str(piter),Nv='100',ranks='64',grid='1x4x4x4',t0='90',tau_min='0',tau_max='16',R_min='1',R_max='32',temporal_bc='open',hyp_steps='1')
    change=float(b['preserved_change'])
    if any(b.get(k)!=v for k,v in expected.items()) or not math.isfinite(change) or not 0<=change<=1e-13:raise ValueError('Changed static geometry or HYP setup')
    temporal_change=float(b['temporal_change'])
    if not math.isfinite(temporal_change) or temporal_change<=1e-13:raise ValueError('Temporal HYP did not change the field')
    values=rows(text,'RT');lookup={}
    for r in values:
        key=(r['scheme'],int(r['tau']),int(r['R']),int(r['axis']))
        if key in lookup or int(r['cfg'])!=cfg or int(r['piter'])!=piter or int(r['Nv'])!=100:raise ValueError('Duplicate/invalid scan metadata')
        if int(r['Nsrc'])!=64**3*(3 if key[3]==-1 else 1):raise ValueError('Wrong spatial source count')
        if int(r['t0'])!=(90 if key[0]=='fixed' else 98-key[1]//2):raise ValueError('Wrong source time')
        if not all(math.isfinite(float(r[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite correlator')
        lookup[key]=r
    expected_keys={('fixed',t,r,a) for t in range(17) for r in range(1,33) for a in (-1,0,1,2)}
    expected_keys|={('regression',t,r,a) for t in (2,4,6) for r in (1,2,4) for a in (-1,0,1,2)}
    if set(lookup)!=expected_keys or len(values)!=2212:raise ValueError('Incomplete static correlator geometry')
    for scheme,t,r,a in expected_keys:
        row=lookup[scheme,t,r,a]
        if a==-1:
            for field in ('L_RE','L_IM'):near(float(row[field]),math.fsum(float(lookup[scheme,t,r,q][field]) for q in (0,1,2))/3)
        if scheme=='regression':
            ref=references[cfg,piter,t,r,a]
            for field in ('L_RE','L_IM'):near(float(row[field]),float(ref[field]))
    direct=rows(text,'DIRECT_TEST');roundtrip=rows(text,'PERIODIC_ROUNDTRIP')
    direct_keys={(r['scheme'],int(r['tau']),int(r['R']),int(r['axis']),int(r['site'])) for r in direct}
    expected_direct={('fixed',t,r,a,s) for t in (0,4,16) for r in (1,4,16,17,31,32) for a in range(3) for s in (0,4095)}
    expected_direct|={('regression',4,r,a,s) for r in (1,4) for a in range(3) for s in (0,4095)}
    if len(direct)!=120 or direct_keys!=expected_direct or any(r['status']!='PASS' or not math.isfinite(float(r['tolerance'])) or not 0<=float(r['abs_error'])<=float(r['tolerance']) for r in direct):
        raise ValueError('Direct/factorized contraction check missing or failed')
    if len(roundtrip)!=3 or {int(r['axis']) for r in roundtrip}!={0,1,2} or any(r['status']!='PASS' or r['cells']!='64' for r in roundtrip):
        raise ValueError('Periodic payload roundtrip failed')
    times=rows(text,'SU3_TIME');window=rows(text,'SU3_WINDOW')
    if len(times)!=16 or {int(r['t']) for r in times}!=set(range(90,106)) or len(window)!=1:raise ValueError('Incomplete transporter SU(3) diagnostic')
    for r in times+window:
        count=64**3*(16 if r is window[0] else 1)
        if int(r['cfg'])!=cfg or int(r['piter'])!=piter or int(r['count'])!=count:raise ValueError('Incorrect SU(3) count')
        for prefix in ('unit','det'):
            maximum,rms=float(r[prefix+'_max']),float(r[prefix+'_rms'])
            if not math.isfinite(maximum) or not math.isfinite(rms) or not 0<=rms<=maximum or maximum>1e-12:
                raise ValueError('Projection does not meet the central-link numerical tolerance')
        if int(r['unit_above_1em2'])!=0 or int(r['det_above_1em2'])!=0:raise ValueError('Inconsistent projection counters')
    if (int(window[0]['t0']),int(window[0]['t1']))!=(90,106):raise ValueError('Wrong transporter window')
    for prefix in ('unit','det'):
        near(float(window[0][prefix+'_max']),max(float(r[prefix+'_max']) for r in times),atol=1e-30)
        near(float(window[0][prefix+'_rms']),math.sqrt(math.fsum(float(r[prefix+'_rms'])**2 for r in times)/16),atol=1e-30)
    blocks=rows(text,'RT_BLOCK');grouped={}
    for r in blocks:
        if int(r['cfg'])!=cfg or int(r['piter'])!=piter or r['t0']!='90' or int(r['Nsrc'])!=4096:raise ValueError('Invalid spatial block metadata')
        key=tuple(int(r[k]) for k in ('tau','R','axis'));grouped.setdefault(key,[]).append(r)
    expected_blocks={(t,r,a) for t in (0,4,8,16) for r in (1,4,8,16,32) for a in (0,1,2)}
    if set(grouped)!=expected_blocks or len(blocks)!=3840:raise ValueError('Incomplete spatial-block diagnostics')
    for key,selected in grouped.items():
        if len(selected)!=64 or {int(r['rank']) for r in selected}!=set(range(64)) or {tuple(int(r[k]) for k in ('px','py','pz')) for r in selected}!={(x,y,z) for x in range(4) for y in range(4) for z in range(4)}:
            raise ValueError('Duplicate or missing spatial block')
        for field in ('L_RE','L_IM'):near(math.fsum(float(r[field]) for r in selected)/64,float(lookup[('fixed',)+key][field]))
    return values,blocks,dict(configuration=cfg,piter=piter,rows=len(values),regression_rows=36,direct_tests=120,periodic_roundtrips=3,spatial_block_rows=len(blocks),unit_max=window[0]['unit_max'],det_max=window[0]['det_max'])


def projection_pair(low,high):
    a={(r['scheme'],r['tau'],r['R'],r['axis']):r for r in low}
    b={(r['scheme'],r['tau'],r['R'],r['axis']):r for r in high}
    if set(a)!=set(b):raise ValueError('Projection controls do not share the geometry')
    maximum=0.
    for k in a:
        for field in ('L_RE','L_IM'):
            x,y=float(a[k][field]),float(b[k][field]);near(x,y,rtol=1e-10)
            maximum=max(maximum,abs(x-y))
    return maximum
