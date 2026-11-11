#!/usr/bin/env python3
"""Reuse exact contraction/I/O helpers and generate strict 66/97-slice audits."""
import sys
from pathlib import Path
from parents import generator,design
from pilot_plan import BLOCK_RADII


def generate(output):
    out=Path(output);generator.generate(out)
    geometry=(out/'matrix_plan.h').read_text()
    old='#define MATRIX_T0 90\n'
    if geometry.count(old)!=1:raise ValueError('Parent source-time definition changed')
    geometry=geometry.replace(old,'/* Actual source time is a required runtime input. */\n')
    old='static const int matrix_block_radii[]={4,32};'
    if geometry.count(old)!=1:raise ValueError('Parent rank-block radii changed')
    geometry=geometry.replace(old,'static const int matrix_block_radii[]={'+','.join(map(str,BLOCK_RADII))+'};')
    geometry+='#define MATRIX_NBLOCK_RADII '+str(len(BLOCK_RADII))+'\n'
    (out/'matrix_plan.h').write_text(geometry)
    audit=(out/'audit_selected.c').read_text()
    old='static const int selected_times[]={'+','.join(map(str,range(90,123)))+'};'
    if audit.count(old)!=1 or '1.1e-9' in audit:raise ValueError('Parent saved-mode audit changed')
    for cfg,size in ((170,66),(650,97)):
        times=design.audit_times(cfg)
        if len(times)!=size:raise ValueError('Unexpected saved-mode audit size')
        new='static const int selected_times[]={'+','.join(map(str,times))+'};'
        (out/('audit%d.c'%size)).write_text(audit.replace(old,new))
    header='/* Exact supported pilot scans; no open-time wrapping. */\n'
    header+='static int source_time_scan_allowed(int cfg,int piter,int t0){\n'
    header+='    if(t0<1||t0+32>=191)return 0;\n'
    for cfg in design.PILOT:
        tests=['(piter==%d&&t0==%d)'%(s['piter'],s['t0']) for s in design.scans(cfg)]
        header+='    if(cfg==%d)return %s;\n'%(cfg,'||'.join(tests))
    header+='    return 0;\n}\n'
    (out/'source_time_plan.h').write_text(header)


if __name__=='__main__':generate(sys.argv[1])
