#!/usr/bin/env python3
"""Check native sources, the unchanged design and all published dependencies."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def verify():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for key,root in (('files',HERE),('dependencies',ROOT)):
        for name,record in manifest[key].items():
            path=root/name
            if path.stat().st_size!=record['bytes'] or digest(path)!=record['sha256']:raise ValueError('Pinned native source/reference changed: '+name)
    print('PASS source-time native package:',len(manifest['files']),'files;',len(manifest['dependencies']),'published dependencies.')


if __name__=='__main__':verify()
