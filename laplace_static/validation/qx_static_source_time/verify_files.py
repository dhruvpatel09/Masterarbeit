#!/usr/bin/env python3
"""Standard-library integrity checks for the design and its published inputs."""
import hashlib,json
from pathlib import Path
from plan import HERE,ROOT,build_plan


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def verify():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for group,root in (('files',HERE),('dependencies',ROOT)):
        for name,record in manifest[group].items():
            p=root/name
            if p.stat().st_size!=record['bytes'] or digest(p)!=record['sha256']:raise ValueError('Pinned design/input changed: '+name)
    if json.loads((HERE/'PLAN.json').read_text())!=build_plan():raise ValueError('Exported plan differs from executable geometry')
    print('PASS source-time design:',len(manifest['files']),'files;',len(manifest['dependencies']),'published dependencies; no native submission command.')


if __name__=='__main__':verify()
