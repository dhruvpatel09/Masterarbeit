#!/usr/bin/env python3
"""Refresh integrity receipts after intentionally regenerating this package."""
import json
from data import HERE, check_inputs, digest


if __name__=='__main__':
    source=check_inputs()
    files={}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name=='FILE_MANIFEST.json' or p.name.endswith('.pending'):
            continue
        files[p.relative_to(HERE).as_posix()] = dict(bytes=p.stat().st_size,sha256=digest(p))
    (HERE/'FILE_MANIFEST.json').write_text(json.dumps(dict(files=files,dependencies=source['files']),indent=2)+'\n')
