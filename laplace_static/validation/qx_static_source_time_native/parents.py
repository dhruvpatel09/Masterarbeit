"""Load pinned parent helpers without leaking their generic module names."""
import importlib.util,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DESIGN=HERE.parent/'qx_static_source_time'
MATRIX=HERE.parent/'qx_static_matrix'


def load(path,name):
    path=Path(path);previous=sys.path[:]
    names=('plan','matrix_plan','measurement_plan')
    saved={key:sys.modules.pop(key,None) for key in names}
    sys.path[:0]=[str(path.parent),str(path.parent.parent)]
    try:
        spec=importlib.util.spec_from_file_location(name,str(path))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
    finally:
        sys.path[:]=previous
        for key,value in saved.items():
            sys.modules.pop(key,None)
            if value is not None:sys.modules[key]=value


design=load(DESIGN/'plan.py','qx_published_source_time_design')
inventory=load(DESIGN/'preflight.py','qx_published_source_time_inventory')
legacy=load(MATRIX/'check_outputs.py','qx_published_matrix_output_helpers')
generator=load(MATRIX/'generate_sources.py','qx_published_matrix_generator')
