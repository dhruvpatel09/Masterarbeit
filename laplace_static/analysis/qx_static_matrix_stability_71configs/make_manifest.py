#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
files={}
for p in sorted(HERE.rglob('*')):
    if not p.is_file():continue
    name=p.relative_to(HERE).as_posix()
    if name=='FILE_MANIFEST.json' or any(s in ('work','__pycache__') for s in p.relative_to(HERE).parts):continue
    files[name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
source=json.loads((HERE/'inputs/SOURCES.json').read_text())
(HERE/'FILE_MANIFEST.json').write_text(json.dumps(dict(files=files,dependencies=source['dependencies']),indent=2)+'\n')
print('Manifest:',len(files),'content files')
