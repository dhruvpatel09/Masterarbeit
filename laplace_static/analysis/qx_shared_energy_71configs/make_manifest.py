#!/usr/bin/env python3
"""Record an intentionally reviewed analysis package."""
import json
from verify_files import HERE,digest
from analysis import write_json


def main():
    files={}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name=='FILE_MANIFEST.json':continue
        if p.name.endswith('.pending'):raise ValueError('Unfinished output '+str(p))
        files[str(p.relative_to(HERE))]={'bytes':p.stat().st_size,'sha256':digest(p)}
    source=json.loads((HERE/'inputs/SOURCES.json').read_text())
    write_json(HERE/'FILE_MANIFEST.json',{'schema':1,'files':files,'dependencies':source['files']})
    print('Recorded %d package files.'%len(files))


if __name__=='__main__':main()
