#!/usr/bin/env python3
"""Verify package and published dependency bytes with the standard library."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def verify():
    manifest = json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for group,base in (('files',HERE),('dependencies',ROOT)):
        for name,receipt in manifest[group].items():
            path = base/name
            if path.stat().st_size!=receipt['bytes'] or digest(path)!=receipt['sha256']:
                raise ValueError('Integrity mismatch: '+name)
    actual = set(p.relative_to(HERE).as_posix() for p in HERE.rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and p.name!='FILE_MANIFEST.json' and not p.name.endswith('.pending'))
    if actual!=set(manifest['files']):
        raise ValueError('Package file inventory mismatch')
    print('PASS qX radius package: {} files; {} published dependencies; SHA-256 and byte counts verified.'.format(len(manifest['files']),len(manifest['dependencies'])))


if __name__=='__main__':verify()
