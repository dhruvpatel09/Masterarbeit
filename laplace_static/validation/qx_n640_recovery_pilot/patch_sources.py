#!/usr/bin/env python3
"""Derive a selected-time writer and independent saved-mode checker."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Source anchor changed: ' + old[:80])
    return text.replace(old, new, 1)


def header(times):
    if times != sorted(set(times)) or not times or min(times) < 0 or max(times) >= 192:
        raise ValueError('Invalid selected times')
    return '''#ifndef QX_RECOVERY_SELECTION_H
#define QX_RECOVERY_SELECTION_H
static const int selected_times[]={%s};
static int qx_selected(int t)
{
    for(unsigned k=0;k<sizeof(selected_times)/sizeof(selected_times[0]);++k)
        if(t==selected_times[k])return 1;
    return 0;
}
#endif
''' % ','.join(map(str, times))


def writer(source):
    text = replace_once(source, '#include <DistRedelemsHdf5Writer.h>',
                        '#include "selection.h"\n#include <DistRedelemsHdf5Writer.h>')
    anchor = '      for( lt=0; lt<geo.lL[0]; lt++ )\n      {\n'
    guard = '''      qcd_errorLoc(params.P[0]!=1 || params.P[1]!=4 || params.P[2]!=4 || params.P[3]!=4 ||
                   numprocs!=64 || params.Nv!=100 ||
                   params.compute_ev!=1 || params.write_ev!=1 || nc!=640 ||
                   eval_only || lpar.res!=1e-9 || lpar.maxRestarts!=30 ||
                   lpar.deg!=10 || lpar.nLanczosFactor!=7 || lpar.upper!=0.15 ||
                   lpar.redTol!=1e-12 || lpar.redFreqFactor!=10 || lpar.tolReortho!=1e-10 ||
                   params.laplace_nAPE!=20 || params.laplace_alphaAPE!=0.5,
                   1,"n640 recovery pilot","Unexpected recovery setup");
      if(myid==0) {printf("RECOVERY_START cfg=640 ranks=64 P=1,4,4,4 selected=%lu res=%.17e\\n",
          (unsigned long)(sizeof(selected_times)/sizeof(selected_times[0])),lpar.res);fflush(stdout);}
'''
    text = replace_once(text, anchor, guard + anchor +
                        '         if(!qx_selected(lt))continue;\n')
    anchor = '            lret=qcd_lanczosCheby(lambda, v, &lstat, &u3d, params.Nv, &lpar);\n'
    text = replace_once(text, anchor, anchor + '''            if(myid==0) {printf("RECOVERY_SOLVE cfg=640 t=%d code=%d restarts=%d nRed=%d maxRes=%.17e lambda_Nv=%.17e\\n",
                lt,lret,lstat.restarts,lstat.nRed,lstat.maxRes,lambda[params.Nv-1]);fflush(stdout);}
''')
    anchor = '''               "DistEigvalsHdf5Writer_WriteDistEigvals failed.");
            }
         }else'''
    text = replace_once(text, anchor, '''               "DistEigvalsHdf5Writer_WriteDistEigvals failed.");
            }
            if(myid==0) {printf("RECOVERY_WRITTEN cfg=640 t=%d Nv=100\\n",lt);fflush(stdout);}
            continue; /* This pilot writes eigenpairs only. */
         }else''')
    return text


def checker(source):
    text = replace_once(source, '#include "solver.h"', '#include "selection.h"')
    text = replace_once(text, 'static const int selected_times[]={0,32,64,95,96,97,98,99,100,101,128,160,191};\n', '')
    text = replace_once(text, '''        seen[t]=1;
        for(int i=0;i<NV;i++) {''', '''        seen[t]=1;
        if(!qx_selected(t))continue; /* Unwritten rows are not valid data. */
        for(int i=0;i<NV;i++) {''')
    text = replace_once(text, '    solver_mpiInit(argc,argv);', '''    int requested;
    if(MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&requested)!=MPI_SUCCESS)return 1;''')
    text = replace_once(text, '''    require(argc!=5 || ranks!=64 || minprovided<MPI_THREAD_FUNNELED ||
            solver_getNThread()!=1,"need 64 ranks, FUNNELED, one thread; CFG GAUGE EV_DIR EVAL",MPI_COMM_WORLD);''',
        '''    require(argc!=5 || ranks!=64 || minprovided<MPI_THREAD_FUNNELED,
            "need 64 ranks, FUNNELED; CFG GAUGE EV_DIR EVAL",MPI_COMM_WORLD);''')
    text = replace_once(text, 'require(*end || nc<100 || nc>800 || nc%10,"invalid configuration",MPI_COMM_WORLD);',
                        'require(*end || nc!=640,"invalid configuration",MPI_COMM_WORLD);')
    text = text.replace('1.1e-9', '1e-9')
    return text


def generate(destination):
    spec = json.loads((HERE/'INPUTS.json').read_text())
    destination = Path(destination)
    destination.mkdir(exist_ok=True)
    for name in ('mental_funneled.c', 'test_saved_modes.c'):
        data = (HERE/'inputs'/name).read_bytes()
        expected = spec['preserved_inputs']['inputs/'+name]
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('Changed preserved source: '+name)
    (destination/'selection.h').write_text(header(spec['selected_times']))
    (destination/'mental_recovery.c').write_text(writer((HERE/'inputs/mental_funneled.c').read_text()))
    (destination/'check_saved.c').write_text(checker((HERE/'inputs/test_saved_modes.c').read_text()))


if __name__ == '__main__':
    import sys
    if len(sys.argv)!=2:raise SystemExit('usage: patch_sources.py BUILD_DIRECTORY')
    generate(sys.argv[1])
