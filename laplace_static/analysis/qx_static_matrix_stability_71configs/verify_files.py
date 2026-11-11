#!/usr/bin/env python3
"""Standard-library integrity verification, including published dependencies."""
import hashlib,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def verify():
    m=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for group,base in (('files',HERE),('dependencies',ROOT)):
        for name,r in m[group].items():
            p=base/name
            if p.stat().st_size!=r['bytes'] or digest(p)!=r['sha256']:raise ValueError('File receipt mismatch '+name)
    print('PASS matrix stability package:',len(m['files']),'files;',len(m['dependencies']),'published dependencies.')


if __name__=='__main__':verify()
