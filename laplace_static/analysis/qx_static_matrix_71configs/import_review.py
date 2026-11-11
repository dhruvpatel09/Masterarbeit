#!/usr/bin/env python3
"""Import the immutable native review as four regular-Git-sized chunks."""
import argparse
import hashlib
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def check_file(path, receipt):
    path = Path(path)
    if path.stat().st_size != receipt['bytes'] or digest(path) != receipt['sha256']:
        raise ValueError('Size or SHA-256 mismatch: ' + str(path))


def import_archive(archive, here=HERE):
    here = Path(here)
    source = json.loads((here/'inputs/SOURCES.json').read_text())
    check_file(archive, dict(bytes=source['archive_bytes'], sha256=source['archive_sha256']))
    # Check every existing destination before writing anything.
    for record in source['chunks']:
        target = here/record['path']
        if target.exists():
            check_file(target, record)
    with Path(archive).open('rb') as stream:
        for record in source['chunks']:
            data = stream.read(record['bytes'])
            if len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
                raise ValueError('Invalid source chunk: ' + record['path'])
            target = here/record['path']
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            pending = target.with_suffix('.pending')
            with pending.open('xb') as out:
                out.write(data)
                out.flush()
                os.fsync(out.fileno())
            pending.replace(target)
        if stream.read(1):
            raise ValueError('Unconsumed archive bytes')
    return len(source['chunks'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True, type=Path)
    args = parser.parse_args()
    print('PASS review import: {} verified chunks; source archive retained.'.format(import_archive(args.archive)))
