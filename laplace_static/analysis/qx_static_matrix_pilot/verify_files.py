#!/usr/bin/env python3
"""Verify the exact package and published dependency bytes using only stdlib."""
import json
from pathlib import Path
from replay import digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def verify():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for prefix,entries in [(HERE,manifest['files']),(ROOT,manifest['published_dependencies'])]:
        for name,record in entries.items():
            path=prefix/name
            if path.stat().st_size!=record['bytes'] or digest(path)!=record['sha256']:raise ValueError('Changed file: '+str(path))
    print('PASS matrix-pilot package: %d files; %d published dependencies; SHA-256 and byte counts verified.'%(len(manifest['files']),len(manifest['published_dependencies'])))


if __name__=='__main__':verify()
