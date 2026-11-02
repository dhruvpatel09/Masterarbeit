#!/usr/bin/env python3
"""Verify the package and published dependencies with the standard library."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def main():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for name,r in manifest['files'].items():
        path=HERE/name
        if path.stat().st_size!=r['bytes'] or digest(path)!=r['sha256']:raise ValueError('Changed package file: '+name)
    actual={str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='FILE_MANIFEST.json'}
    if actual!=set(manifest['files']):raise ValueError('Unexpected or missing package files: '+str(actual.symmetric_difference(manifest['files'])))
    source=json.loads((HERE/'inputs/SOURCES.json').read_text())
    if source['files']!=manifest['dependencies']:raise ValueError('Changed dependency list')
    for name,expected in source['files'].items():
        if digest(ROOT/name)!=expected:raise ValueError('Changed published dependency: '+name)
    print('PASS shared-energy package: %d files; %d published dependencies; SHA-256 and byte counts verified.'%(len(manifest['files']),len(source['files'])))


if __name__=='__main__':main()
