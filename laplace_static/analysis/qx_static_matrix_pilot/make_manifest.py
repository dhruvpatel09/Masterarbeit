#!/usr/bin/env python3
"""Record the validated package and byte-exact published validation dependencies."""
import json
from pathlib import Path
from replay import digest
from run_analysis import atomic_text

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def make():
    paths=[p for p in HERE.rglob('*') if p.is_file() and p.name!='FILE_MANIFEST.json' and '__pycache__' not in p.parts and not p.name.endswith('.pending')]
    parent='laplace_static/validation/qx_static_matrix/INPUTS.json';dependencies=json.loads((ROOT/parent).read_text())['files']
    dependencies=dict(dependencies);dependencies[parent]=digest(ROOT/parent)
    for name,expected in dependencies.items():
        if digest(ROOT/name)!=expected:raise ValueError('Published dependency does not match its source pin: '+name)
    record=lambda p:dict(sha256=digest(p),bytes=p.stat().st_size)
    data=dict(base_commit='646175fd325c07ec5aafcb8cbb0322be4c3ba9dc',files={str(p.relative_to(HERE)):record(p) for p in sorted(paths)},
        published_dependencies={name:record(ROOT/name) for name in sorted(dependencies)})
    atomic_text(HERE/'FILE_MANIFEST.json',json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('Recorded %d package files and %d published dependencies.'%(len(paths),len(dependencies)))


if __name__=='__main__':make()
