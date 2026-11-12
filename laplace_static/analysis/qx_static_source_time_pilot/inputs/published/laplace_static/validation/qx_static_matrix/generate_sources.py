#!/usr/bin/env python3
"""Reuse the pinned 100-mode I/O, MPI payload and strict saved-mode audit."""
import importlib.util
import sys
from pathlib import Path
from matrix_plan import CUTOFFS,RADII,DIRECT_TAUS,BLOCK_TAUS,BLOCK_RADII,expected_keys

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'qx_static_extension/generate_sources.py'


def load_parent():
    previous=sys.path[:];sys.path.insert(0,str(PARENT.parent))
    try:
        spec=importlib.util.spec_from_file_location('qx_matrix_parent_generator',str(PARENT))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module
    finally:sys.path[:]=previous


def generate(output):
    out=Path(output);load_parent().generate(out)
    text='/* Generated pilot geometry; source and sink cutoffs are independent. */\n'
    text+='#define MATRIX_NCUT 4\n#define MATRIX_T0 90\n#define MATRIX_TMAX 32\n#define MATRIX_RECORDS %d\n'%len(expected_keys())
    for name,values in [('matrix_cutoffs',CUTOFFS),('matrix_radii',RADII),('matrix_direct_taus',DIRECT_TAUS),('matrix_block_taus',BLOCK_TAUS),('matrix_block_radii',BLOCK_RADII)]:
        text+='static const int %s[]={%s};\n'%(name,','.join(map(str,values)))
    text+='static int matrix_contains(const int *v,int n,int x){for(int i=0;i<n;i++)if(v[i]==x)return 1;return 0;}\n'
    (out/'matrix_plan.h').write_text(text)


if __name__=='__main__':generate(sys.argv[1])
