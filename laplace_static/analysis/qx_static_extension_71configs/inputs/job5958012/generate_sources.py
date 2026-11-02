#!/usr/bin/env python3
"""Extract pinned parent kernels; vary only the number of summed saved modes."""
import importlib.util
from pathlib import Path

from measurement_plan import AUDIT_TIMES, CUTOFFS, RADII

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'qx_static_potential/generate_helpers.py'
spec=importlib.util.spec_from_file_location('qx_static_parent_helpers',str(PARENT))
parent=importlib.util.module_from_spec(spec)
spec.loader.exec_module(parent)


def generate(output):
    out=Path(output)
    parent.generate(out)
    kernel=(out/'kernels.h').read_text()
    # HDF5 buffers and MPI payloads always retain all 100 saved modes. Only
    # the two contraction sums change; the O(Nv^2) check stays independent.
    bounds=[('n<NV','n<active_modes'),('i<NV','i<active_modes'),
            ('j<NV','j<active_modes')]
    for old,new in bounds:
        if kernel.count(old)!=1:raise ValueError('Parent contraction changed: '+old)
        kernel=kernel.replace(old,new)
    (out/'kernels.h').write_text('/* First active_modes saved eigenvectors; I/O size remains NV. */\n'
        'static int active_modes=NV;\n'+kernel)
    audit=(out/'audit_selected.c').read_text()
    old='static const int selected_times[]={'+','.join(map(str,range(90,107)))+'};'
    new='static const int selected_times[]={'+','.join(map(str,AUDIT_TIMES))+'};'
    if audit.count(old)!=1:raise ValueError('Parent audit geometry changed')
    (out/'audit_selected.c').write_text(audit.replace(old,new))
    (out/'measurement_plan.h').write_text(
        '/* Generated from measurement_plan.py. */\n'
        'static const int mode_cutoffs[]={'+','.join(map(str,CUTOFFS))+'};\n'
        'static int selected_radius(int r) {\n'
        '    static const int radii[]={'+','.join(map(str,RADII))+'};\n'
        '    for(unsigned int i=0;i<sizeof(radii)/sizeof(radii[0]);i++)if(r==radii[i])return 1;\n'
        '    return 0;\n}\n')


if __name__=='__main__':
    import sys
    generate(sys.argv[1])
