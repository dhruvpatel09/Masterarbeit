#!/usr/bin/env python3
"""Receipt for this design, its reference data and unchanged published inputs."""
import hashlib,json
from plan import HERE
files={}
for p in sorted(HERE.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(HERE)
    if rel.as_posix()=='FILE_MANIFEST.json' or any(s in ('work','__pycache__') for s in rel.parts):continue
    files[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
source=json.loads((HERE/'inputs/SOURCE_PINS.json').read_text())
(HERE/'FILE_MANIFEST.json').write_text(json.dumps(dict(files=files,dependencies=source['dependencies']),indent=2)+'\n')
print('Manifest:',len(files),'design/reference files')
