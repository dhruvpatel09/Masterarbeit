#!/usr/bin/env python3
"""Optionally replay the full stdout archive using the preserved measurement checkers."""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path
import reproduce as a


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    args = parser.parse_args()
    sources, report = a.verify_inputs()
    if a.sha(args.archive) != sources['review_archive_sha256']:
        raise ValueError('Review archive differs from the preserved source')
    for name, expected in report['parser_source_sha256'].items():
        if a.sha(a.ROOT/name) != expected:
            raise ValueError('Published measurement checker changed: '+name)
    old = a.ROOT/'laplace_static/validation/hyp_projection_configs_qx_200_230'
    sys.path.insert(0,str(old))
    import check_task as check
    spec = importlib.util.spec_from_file_location('clover_trace_pilot_replay',
              str(a.ROOT/'laplace_static/validation/qx_clover_trace/campaign.py'))
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    expected, expected_windows = defaultdict(list), defaultdict(list)
    for row in a.read_csv(a.HERE/'inputs/profiles_raw_and_trace_70configs.csv'):
        expected[int(row['configuration']),int(row['tau'])].append(row)
    for row in a.read_csv(a.HERE/'inputs/su3_windows_70configs.csv'):
        expected_windows[int(row['configuration']),int(row['tau'])].append(row)
    base = 'laplace_static/validation/qx_clover_trace_ensemble/campaigns/full70/'
    with zipfile.ZipFile(str(args.archive)) as archive, tempfile.TemporaryDirectory() as temp:
        if len(archive.namelist()) != len(set(archive.namelist())):
            raise ValueError('Duplicate archive paths')
        for name in ('PREPARED.json','prepared.sha256','SUBMISSION_001.json','SUMMARY.json'):
            if archive.read(base+name) != (a.HERE/'inputs/job5943513'/name).read_bytes():
                raise ValueError('Campaign metadata changed: '+name)
        for name, target in [('profiles_raw_and_trace.csv','profiles_raw_and_trace_70configs.csv'),
                             ('su3_windows.csv','su3_windows_70configs.csv')]:
            original = archive.read(base+name).decode().replace('\r\n','\n')
            if original != (a.HERE/'inputs'/target).read_text():
                raise ValueError('Collected campaign primary table changed')
        plan = json.loads(archive.read(base+'PREPARED.json'))
        if plan['stage'] != 'full70' or len(plan['cases']) != 209 or len(plan['slots']) != 35:
            raise ValueError('Wrong full-cohort plan')
        cases = list(plan['cases'])+[dict(index=None,configuration=100,tau=4,piters=[4,40])]
        count = 0
        maximum_identity = 0.0
        for entry in cases:
            cfg,tau = entry['configuration'],entry['tau']
            if entry['index'] is None:
                seed = a.ROOT/'laplace_static/analysis/qx_clover_trace_pilot/inputs/job5942805'
                logs = seed/'outputs'
                launches = json.loads((seed/'LAUNCHES.json').read_text())
            else:
                folder = base+'case_%03d/' % entry['index']
                latest = json.loads(archive.read(folder+'LATEST.json'))
                if latest['attempt'] != 'attempt_0001':
                    raise ValueError('Unexpected attempt in this preserved campaign')
                prefix = folder+latest['attempt']+'/'
                summary = archive.read(prefix+'SUMMARY.json')
                if hashlib.sha256(summary).hexdigest() != latest['summary_sha256']:
                    raise ValueError('Case completion checksum differs')
                wanted = dict(success=True,index=entry['index'],configuration=cfg,tau=tau,
                              piters=[4,40],raw_records=664,window_records=2)
                if json.loads(summary) != wanted:
                    raise ValueError('Case validation metadata differs')
                marker = json.loads(archive.read(prefix+'CASE.json'))
                if marker['prepared_sha256'] != a.sha(a.HERE/'inputs/job5943513/PREPARED.json'):
                    raise ValueError('Case preparation checksum differs')
                launches = json.loads(archive.read(prefix+'LAUNCHES.json'))
                logs = Path(temp)
                for p in (4,40):
                    if archive.read(prefix+'p%d.err' % p):
                        raise ValueError('Nonempty launch stderr')
                    if json.loads(archive.read(prefix+'launch_p%d.json' % p)) != launches[(0 if p == 4 else 1)]:
                        raise ValueError('Individual launch record differs')
                    (logs/('p%d.out' % p)).write_bytes(archive.read(prefix+'p%d.out' % p))
            if [(r['configuration'],r['tau'],r['piter'],r['exit_code']) for r in launches] != [(cfg,tau,p,0) for p in (4,40)]:
                raise ValueError('Incomplete launch inventory')
            raw, windows = [], []
            for p in (4,40):
                path = logs/('p%d.out' % p)
                with contextlib.redirect_stdout(io.StringIO()):
                    _,scan,_,window = check.read_measurement(path,cfg,tau,p)
                text = path.read_text()
                traces = pilot.trace_records(text,scan,check.avg)
                maximum_identity = max(maximum_identity,float(check.avg.rows(text,'CLOVER_TRACE_CHECK')[0]['identity_max']))
                for row in check.raw_rows(cfg,tau,p,scan):
                    key = tuple(int(row[k]) for k in ('axis','R','long','tr1','tr2'))
                    row.update({k:traces[key][k] for k in pilot.TRACE_FIELDS})
                    raw.append({k:str(v) for k,v in row.items()})
                windows.append({k:str(v) for k,v in window.items()})
            if raw != expected[cfg,tau] or windows != expected_windows[cfg,tau]:
                raise ValueError('Reparsed primary or window table differs at '+str((cfg,tau)))
            count += 1
            if count % 30 == 0:
                print('PASS replay',count,'of 210 cases',flush=True)
        if maximum_identity != report['maximum_direct_clover_identity_residual']:
            raise ValueError('Direct clover identity diagnostic differs')
    print('PASS review: all 210 cases and 420 stdout measurements reproduce the preserved primary tables.')


if __name__ == '__main__':
    main()
