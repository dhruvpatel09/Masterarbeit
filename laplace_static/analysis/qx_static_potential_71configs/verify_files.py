#!/usr/bin/env python3
"""Standard-library integrity check for the archived full71 result package."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    if data['configuration_labels']!=list(range(100,801,10)):raise ValueError('Wrong cohort')
    for name,expected in data['files'].items():
        path=(HERE/name).resolve()
        if HERE not in path.parents:raise ValueError('Invalid manifest path')
        h=hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
        if h.hexdigest()!=expected:raise ValueError('File changed: '+name)
    print('PASS package integrity: %d files; complete 71-configuration result set.'%len(data['files']))


if __name__=='__main__':main()
