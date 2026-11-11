#!/usr/bin/env python3
"""Check package, four immutable source chunks and published dependencies."""
import json
from import_review import HERE, check_file

ROOT=HERE.parents[2]


def verify():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for group,base in (('files',HERE),('source_chunks',HERE),('dependencies',ROOT)):
        for name,receipt in manifest[group].items():
            check_file(base/name,receipt)
    print('PASS full71 matrix package: {} files; four source chunks; {} published dependencies.'.format(
        len(manifest['files']),len(manifest['dependencies'])))


if __name__=='__main__':verify()
