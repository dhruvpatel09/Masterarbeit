#!/usr/bin/env python3
"""Replay temporal/cutoff measurements against archived matched correlators."""
import math
from pathlib import Path

from measurement_plan import (AUDIT_TIMES, BLOCK_RADII, BLOCK_TAUS, RADII,
                              cutoffs, direct_keys, expected_keys)


def rows(text,prefix):
    return [dict(token.split('=',1) for token in line.split()[1:])
            for line in text.splitlines() if line.startswith(prefix+' ')]


def near(a,b,rtol=5e-12,atol=1e-20):
    if not math.isfinite(a) or not math.isfinite(b) or abs(a-b)>atol+rtol*max(abs(a),abs(b)):
        raise ValueError('Numerical comparison failed: %.17e versus %.17e'%(a,b))


def audit(path,cfg):
    text=Path(path).read_text();modes=rows(text,'MODE');slices=rows(text,'SLICE');end=rows(text,'QUALITY_END')
    if len(end)!=1 or end[0]!=dict(cfg=str(cfg),slices='33',status='PASS'):
        raise ValueError('Incomplete 33-slice saved-mode audit')
    expected={(t,i) for t in AUDIT_TIMES for i in range(100)}
    if len(modes)!=3300 or {(int(r['t']),int(r['i'])) for r in modes}!=expected:
        raise ValueError('Wrong saved-mode audit coverage')
    eigenvalues={}
    for r in modes:
        values=[float(r[k]) for k in ('lambda','norm','residual_l2','residual_relative','rayleigh_re','rayleigh_im')]
        if int(r['cfg'])!=cfg or not all(math.isfinite(x) for x in values) or not 0<values[0]<=12.00000001:
            raise ValueError('Invalid mode audit value')
        if values[2]<0 or values[3]<0 or values[3]>1e-9 or abs(values[1]-1)>1e-9 or abs(values[5])>1e-10:
            raise ValueError('Saved mode fails unchanged residual/norm criteria')
        eigenvalues[int(r['t']),int(r['i'])]=values[0]
    # The cutoff must select the lowest saved modes in their stored order.
    for t in AUDIT_TIMES:
        if any(eigenvalues[t,i+1]<eigenvalues[t,i]-1e-12 for i in range(99)):
            raise ValueError('Saved modes are not ordered by eigenvalue')
    if len(slices)!=33 or {int(r['t']) for r in slices}!=set(AUDIT_TIMES):
        raise ValueError('Missing slice diagnostics')
    for r in slices:
        if int(r['cfg'])!=cfg or r['status']!='PASS' or int(r['Nv'])!=100:raise ValueError('Failed slice')
        for field,limit in [('max_relative_residual',1e-9),('max_norm_defect',1e-9),('gram_diag_defect',1e-9),('gram_offdiag',1e-9),('max_rayleigh_im',1e-10)]:
            value=float(r[field])
            if not math.isfinite(value) or not 0<=value<=limit:raise ValueError('Invalid slice quality')
    return dict(configuration=cfg,slices=33,modes=3300,
        maximum_relative_residual=max(float(r['residual_relative']) for r in modes),
        saved_eigenvalue_order_checked=True)


def scan(path,cfg,piter,references):
    text=Path(path).read_text();begin=rows(text,'EXTENSION_BEGIN');end=rows(text,'EXTENSION_END')
    keys=expected_keys(piter)
    if len(begin)!=1 or len(end)!=1 or end[0]!=dict(cfg=str(cfg),piter=str(piter),records=str(len(keys)),status='PASS'):
        raise ValueError('Incomplete temporal/cutoff scan')
    b=begin[0]
    expected=dict(cfg=str(cfg),piter=str(piter),saved_Nv='100',cutoffs=','.join(map(str,cutoffs(piter))),
        ranks='64',grid='1x4x4x4',t0='90',tau_min='0',tau_max='32',
        radii=','.join(map(str,RADII)),temporal_bc='open',hyp_steps='1')
    change=float(b['preserved_change']);temporal_change=float(b['temporal_change'])
    if any(b.get(k)!=v for k,v in expected.items()) or not math.isfinite(change) or not 0<=change<=1e-13:
        raise ValueError('Changed scan geometry or HYP setup')
    if not math.isfinite(temporal_change) or temporal_change<=1e-13:raise ValueError('Temporal HYP did not change the field')
    values=rows(text,'RT');lookup={}
    for r in values:
        key=(r['scheme'],int(r['Nv']),int(r['tau']),int(r['R']),int(r['axis']))
        if key in lookup or int(r['cfg'])!=cfg or int(r['piter'])!=piter:raise ValueError('Duplicate/invalid scan metadata')
        if int(r['Nsrc'])!=64**3*(3 if key[-1]==-1 else 1):raise ValueError('Wrong spatial source count')
        if int(r['t0'])!=(90 if key[0]=='fixed' else 98-key[2]//2):raise ValueError('Wrong source time')
        if not all(math.isfinite(float(r[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite correlator')
        lookup[key]=r
    if set(lookup)!=keys or len(values)!=len(keys):raise ValueError('Incomplete temporal/cutoff geometry')
    baseline_count=0;baseline_max=0.
    for scheme,n,t,r,a in keys:
        row=lookup[scheme,n,t,r,a]
        if a==-1:
            for field in ('L_RE','L_IM'):
                near(float(row[field]),math.fsum(float(lookup[scheme,n,t,r,q][field]) for q in (0,1,2))/3)
        if n==100 and (scheme=='regression' or t<=16):
            ref=references[cfg,piter,scheme,n,t,r,a];baseline_count+=1
            for field in ('L_RE','L_IM'):
                x,y=float(row[field]),float(ref[field]);near(x,y);baseline_max=max(baseline_max,abs(x-y))
        if r==32:
            if abs(float(row['L_IM']))>1e-18+1e-10*abs(float(row['L_RE'])):
                raise ValueError('Half-box imaginary cancellation failed')
    direct=rows(text,'DIRECT_TEST');roundtrip=rows(text,'PERIODIC_ROUNDTRIP')
    direct_inventory={(r['scheme'],int(r['Nv']),int(r['tau']),int(r['R']),int(r['axis']),int(r['site'])) for r in direct}
    if len(direct)!=len(direct_keys(piter)) or direct_inventory!=direct_keys(piter):
        raise ValueError('Missing or duplicate direct-contraction checks')
    for r in direct:
        error,tolerance=float(r['abs_error']),float(r['tolerance'])
        if r['status']!='PASS' or not math.isfinite(error) or not math.isfinite(tolerance) or not 0<=error<=tolerance:
            raise ValueError('Direct/factorized contraction failed')
    if len(roundtrip)!=3 or {int(r['axis']) for r in roundtrip}!={0,1,2} or any(r['status']!='PASS' or r['cells']!='64' for r in roundtrip):
        raise ValueError('Periodic payload roundtrip failed')
    times=rows(text,'SU3_TIME');window=rows(text,'SU3_WINDOW')
    if len(times)!=32 or {int(r['t']) for r in times}!=set(range(90,122)) or len(window)!=1:
        raise ValueError('Incomplete extended transporter SU(3) diagnostic')
    for r in times+window:
        count=64**3*(32 if r is window[0] else 1)
        if int(r['cfg'])!=cfg or int(r['piter'])!=piter or int(r['count'])!=count:raise ValueError('Incorrect SU(3) count')
        for prefix in ('unit','det'):
            maximum,rms=float(r[prefix+'_max']),float(r[prefix+'_rms'])
            if not math.isfinite(maximum) or not math.isfinite(rms) or not 0<=rms<=maximum or maximum>1e-12:
                raise ValueError('Projection does not meet the central-link numerical tolerance')
        if int(r['unit_above_1em2'])!=0 or int(r['det_above_1em2'])!=0:raise ValueError('Inconsistent projection counters')
    if (int(window[0]['t0']),int(window[0]['t1']))!=(90,122):raise ValueError('Wrong transporter window')
    for prefix in ('unit','det'):
        near(float(window[0][prefix+'_max']),max(float(r[prefix+'_max']) for r in times),atol=1e-30)
        near(float(window[0][prefix+'_rms']),math.sqrt(math.fsum(float(r[prefix+'_rms'])**2 for r in times)/32),atol=1e-30)
    blocks=rows(text,'RT_BLOCK');grouped={}
    for r in blocks:
        if int(r['cfg'])!=cfg or int(r['piter'])!=piter or r['Nv']!='100' or r['t0']!='90' or int(r['Nsrc'])!=4096:
            raise ValueError('Invalid spatial block metadata')
        if not all(math.isfinite(float(r[k])) for k in ('L_RE','L_IM')):raise ValueError('Nonfinite spatial block')
        key=tuple(int(r[k]) for k in ('tau','R','axis'));grouped.setdefault(key,[]).append(r)
    expected_blocks={(t,r,a) for t in BLOCK_TAUS for r in BLOCK_RADII for a in (0,1,2)}
    if set(grouped)!=expected_blocks or len(blocks)!=4800:raise ValueError('Incomplete spatial-block diagnostics')
    for key,selected in grouped.items():
        if len(selected)!=64 or {int(r['rank']) for r in selected}!=set(range(64)) or {tuple(int(r[k]) for k in ('px','py','pz')) for r in selected}!={(x,y,z) for x in range(4) for y in range(4) for z in range(4)}:
            raise ValueError('Duplicate or missing spatial block')
        for field in ('L_RE','L_IM'):
            near(math.fsum(float(r[field]) for r in selected)/64,float(lookup[('fixed',100)+key][field]))
    return values,blocks,dict(configuration=cfg,piter=piter,rows=len(values),regression_rows=36,
        archived_reference_rows=baseline_count,archived_reference_maximum_absolute_difference=baseline_max,
        cutoffs=cutoffs(piter),direct_tests=len(direct),periodic_roundtrips=3,
        spatial_block_rows=len(blocks),unit_max=window[0]['unit_max'],det_max=window[0]['det_max'])


def projection_pair(low,high):
    def subset(values):
        return {(r['cfg'],r['scheme'],r['Nv'],r['tau'],r['R'],r['axis']):r
                for r in values if int(r['Nv'])==100}
    a,b=subset(low),subset(high)
    if not a or set(a)!=set(b):raise ValueError('Projection controls do not share the Nv=100 geometry')
    maximum=0.
    for key in a:
        for field in ('L_RE','L_IM'):
            x,y=float(a[key][field]),float(b[key][field]);near(x,y,rtol=1e-10)
            maximum=max(maximum,abs(x-y))
    return maximum
