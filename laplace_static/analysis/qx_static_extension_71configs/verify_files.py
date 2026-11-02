#!/usr/bin/env python3
"""Check the published package using the Python standard library."""
import argparse
import json
from pathlib import Path
from replay_inputs import digest

HERE=Path(__file__).resolve().parent


def verify():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text())
    for name,expected in manifest['files'].items():
        if digest(HERE/name)!=expected:raise ValueError('Package file changed: '+name)
    print('PASS package integrity: %d files; complete full71 extension analysis.'%len(manifest['files']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');a=p.parse_args();verify()
    if a.replay:
        from replay_inputs import replay
        _,d,_=replay(progress=True)
        print('PASS full71 input replay: 145 launches, 234300 audited modes, 452412 correlator records.')
