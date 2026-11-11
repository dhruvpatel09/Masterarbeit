#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
files={}
for path in sorted(HERE.rglob('*')):
    if not path.is_file():continue
    rel=path.relative_to(HERE)
    if rel.as_posix()=='FILE_MANIFEST.json' or any(p in ('build','campaigns','__pycache__') for p in rel.parts):continue
    if path.name.startswith('source_time_') and path.suffix in ('.out','.err'):continue
    files[rel.as_posix()]=dict(bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
source=json.loads((HERE/'SOURCE_PINS.json').read_text())
(HERE/'FILE_MANIFEST.json').write_text(json.dumps(dict(files=files,dependencies=source['dependencies']),indent=2)+'\n')
print('Manifest:',len(files),'native source files')
