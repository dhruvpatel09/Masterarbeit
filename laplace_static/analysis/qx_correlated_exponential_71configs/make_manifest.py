#!/usr/bin/env python3
"""Record an intentionally reviewed set of analysis files."""
import json
from verify_files import HERE,digest
from run_analysis import write_json


def main():
    files={}
    for path in sorted(HERE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.name=='FILE_MANIFEST.json':continue
        if path.name.endswith('.pending'):raise ValueError('Unfinished output: '+str(path))
        files[str(path.relative_to(HERE))]={'bytes':path.stat().st_size,'sha256':digest(path)}
    source=json.loads((HERE/'inputs/SOURCES.json').read_text())
    write_json(HERE/'FILE_MANIFEST.json',{'schema':1,'files':files,'dependencies':source['files']})
    print('Recorded %d package files and two published dependencies.'%len(files))


if __name__=='__main__':main()
