#!/usr/bin/env python3
"""Refresh package receipts after a reviewed analysis change."""
import json
from import_review import HERE,digest

source=json.loads((HERE/'inputs/SOURCES.json').read_text())
files={}
for p in sorted(HERE.rglob('*')):
    if not p.is_file():continue
    name=p.relative_to(HERE).as_posix()
    if name=='FILE_MANIFEST.json' or any(part in ('work','__pycache__','review_chunks') for part in p.relative_to(HERE).parts):continue
    files[name]=dict(bytes=p.stat().st_size,sha256=digest(p))
manifest=dict(files=files,source_chunks={r['path']:{k:r[k] for k in ('bytes','sha256')} for r in source['chunks']},
              dependencies=source['dependencies'])
(HERE/'FILE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Manifest:',len(files),'files')
