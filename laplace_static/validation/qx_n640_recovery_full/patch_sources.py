#!/usr/bin/env python3
"""Select disjoint time blocks in the successful saved-mode pilot sources."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Changed source anchor: '+old[:80])
    return text.replace(old,new,1)


def header(reused):
    return '''#ifndef QX_RECOVERY_SELECTION_H
#define QX_RECOVERY_SELECTION_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int qx_slot(void)
{
    static int slot=-3;
    if(slot==-3) {
        const char *s=getenv("QX_RECOVERY_SLOT");char *end;
        if(s && strcmp(s,"all")==0)slot=-1;
        else if(s && *s) {
            long value=strtol(s,&end,10);
            if(*end || value<0 || value>15) {fprintf(stderr,"Invalid QX_RECOVERY_SLOT\\n");exit(2);}
            slot=(int)value;
        } else {fprintf(stderr,"Missing QX_RECOVERY_SLOT\\n");exit(2);}
    }
    return slot;
}
static int qx_selected(int t)
{
    static const int reused[]={%s};
    int slot=qx_slot();
    if(t<0 || t>=192)return 0;
    if(slot==-1)return 1;
    if(t/12!=slot)return 0;
    for(unsigned i=0;i<sizeof(reused)/sizeof(reused[0]);++i)
        if(t==reused[i])return 0;
    return 1;
}
static int qx_selected_count(void)
{
    int count=0;for(int t=0;t<192;++t)count+=qx_selected(t);return count;
}
#endif
''' % ','.join(map(str,reused))


def transform_writer(source):
    text=once(source,'(unsigned long)(sizeof(selected_times)/sizeof(selected_times[0]))',
              '(unsigned long)qx_selected_count()')
    return once(text,'qcd_errorLoc(params.P[0]!=1','qcd_errorLoc(qx_slot()<0 || params.P[0]!=1')


def transform_checker(source):
    old='''static int selected(int t)
{
    for(size_t k=0;k<sizeof(selected_times)/sizeof(selected_times[0]);k++)
        if(t==selected_times[k])return 1;
    return 0;
}'''
    text=once(source,old,'static int selected(int t) {return qx_selected(t);}')
    return once(text,'(int)(sizeof(selected_times)/sizeof(selected_times[0]))','qx_selected_count()')


def generate(destination):
    spec=json.loads((HERE/'INPUTS.json').read_text());destination=Path(destination)
    destination.mkdir(exist_ok=True)
    for name in ('mental_recovery.c','check_saved.c'):
        if hashlib.sha256((HERE/'inputs'/name).read_bytes()).hexdigest()!=spec['preserved_inputs']['inputs/'+name]:
            raise ValueError('Preserved successful source changed')
    (destination/'selection.h').write_text(header(spec['reused_times']))
    (destination/'mental_recovery.c').write_text(transform_writer((HERE/'inputs/mental_recovery.c').read_text()))
    (destination/'check_saved.c').write_text(transform_checker((HERE/'inputs/check_saved.c').read_text()))


if __name__=='__main__':
    import sys
    if len(sys.argv)!=2:raise SystemExit('usage: patch_sources.py BUILD_DIRECTORY')
    generate(sys.argv[1])
