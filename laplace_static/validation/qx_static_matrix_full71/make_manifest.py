#!/usr/bin/env python3
"""Refresh receipts after intentionally editing the full71 workflow."""
import json
from campaign import HERE, ROOT, digest


if __name__=='__main__':
    spec=json.loads((HERE/'INPUTS.json').read_text())
    dependencies={}
    for name,expected in sorted(spec['files'].items()):
        path=ROOT/name
        if digest(path)!=expected:raise ValueError('Native source/reference changed: '+name)
        dependencies[name]=dict(bytes=path.stat().st_size,sha256=expected)
    files={}
    for p in sorted(HERE.iterdir()):
        if p.is_file() and p.name!='FILE_MANIFEST.json' and not p.name.endswith(('.pending','.out','.err')):
            files[p.name]=dict(bytes=p.stat().st_size,sha256=digest(p))
    (HERE/'FILE_MANIFEST.json').write_text(json.dumps(dict(files=files,dependencies=dependencies),indent=2)+'\n')
