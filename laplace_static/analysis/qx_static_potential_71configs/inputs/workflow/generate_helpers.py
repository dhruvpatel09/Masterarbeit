#!/usr/bin/env python3
"""Reuse the validated qX correlator kernels without modifying their parent source."""
import re
import shutil
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PARENT=ROOT/'laplace_static/validation/qx_clover_trace/test_profiles_trace_qx.c'


def function(text,name):
    match=re.search(r'^static [^;{}]*\b'+re.escape(name)+r'\s*\(',text,re.M)
    if not match:raise ValueError('Missing parent function '+name)
    start=match.start();brace=text.find('{',match.end());depth=0;state='code';i=brace
    while i<len(text):
        c=text[i];n=text[i+1:i+2]
        if state=='code':
            if c=='/' and n=='*':state='comment';i+=2;continue
            if c=='/' and n=='/':state='line';i+=2;continue
            if c in ('"',"'"):state=c
            elif c=='{':depth+=1
            elif c=='}':
                depth-=1
                if not depth:return text[start:i+1]+'\n'
        elif state=='comment':
            if c=='*' and n=='/':state='code';i+=2;continue
        elif state=='line':
            if c=='\n':state='code'
        elif c=='\\':i+=2;continue
        elif c==state:state='code'
        i+=1
    raise ValueError('Unclosed parent function '+name)


def generate(output):
    out=Path(output);out.mkdir(parents=True,exist_ok=True);s=PARENT.read_text()
    packet=s[s.index('typedef struct {\n    qcd_complex_16 v[2][NV][3];'):s.index('static qcd_complex_16 packet_pair(')]
    pure=['mat3_zero','mat3_mul','mat3_adjoint','mat3_trace']
    (out/'kernels.h').write_text('/* Exact kernels extracted from the pinned qX profile source. */\n'+
        '\n'.join(function(s,n) for n in pure)+packet+'\n'+
        '\n'.join(function(s,n) for n in ('packet_pair','packet_tau','packet_direct')))
    (out/'io_helpers.h').write_text('\n'.join(function(s,n) for n in ('read_env_int','env_or_default','read_evec_time','temporal_product')))
    (out/'mpi_helpers.h').write_text('\n'.join(function(s,n) for n in ('rank_neighbours','shift_endpoint','verify_shift','kahan_add')))
    source=ROOT/'laplace_static/validation/qx_n640_recovery_pilot/inputs/test_saved_modes.c'
    audit=source.read_text()
    old='static const int selected_times[]={0,32,64,95,96,97,98,99,100,101,128,160,191};'
    new='static const int selected_times[]={'+','.join(map(str,range(90,107)))+'};'
    if audit.count(old)!=1 or audit.count('1.1e-9')!=2:raise ValueError('Saved-mode audit template changed')
    audit=audit.replace(old,new).replace('1.1e-9','1e-9')
    (out/'audit_selected.c').write_text(audit)
    shutil.copyfile(ROOT/'laplace_static/validation/qx_n640_recovery_pilot/metrics.h',out/'metrics.h')


if __name__=='__main__':
    import sys
    generate(sys.argv[1])
