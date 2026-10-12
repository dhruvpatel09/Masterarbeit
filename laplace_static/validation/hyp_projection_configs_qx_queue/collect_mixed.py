#!/usr/bin/env python3
"""Recheck fifteen cases mapped to distinct Slurm parents and make a review ZIP."""
import argparse
import csv
import fcntl
import json
import zipfile
from pathlib import Path
import queue_common as q


def exports(entry):
    summary = q.validate_case(entry)
    folder = q.output_directory(entry)
    rows = {}
    for name, count in [('profiles_raw.csv',664), ('profiles_paired.csv',332),
                        ('su3_times.csv',2*entry['tau']), ('su3_windows.csv',2)]:
        with (folder / name).open(newline='') as stream:
            items = list(csv.DictReader(stream))
        if len(items) != count:
            raise RuntimeError('Wrong export row count: '+str(folder / name))
        rows[name] = items
    return summary, rows


def collect(data, partial=False):
    q.validate_ledger(data)
    mapping = {str(e['index']):e['array_job_id'] for e in data['cases']}
    name = 'mixed_'+data['first_job']
    if partial:
        name = 'partial_'+name+'_'+q.now().replace(':','').replace('+','_')
    target = q.HERE / ('collection_'+name)
    archive = q.c.ROOT / ('hyp_projection_configs_review_'+name+'.zip')
    if target.exists() or archive.exists():
        if not partial and target.is_dir() and archive.is_file():
            saved = json.loads((target / 'SUMMARY.json').read_text())
            if saved.get('success') is True and saved.get('job_mapping') == mapping:
                # Resume after a completed collection: verify records and archive.
                for entry in data['cases']:
                    exports(entry)
                if (target / 'archive.sha256').read_text().split()[0] != q.c.digest(archive):
                    raise RuntimeError('Existing review archive changed')
                with zipfile.ZipFile(archive) as bundle:
                    if bundle.testzip() is not None:
                        raise RuntimeError('Existing review archive is damaged')
                q.say('REVIEW_ARCHIVE', path=str(archive), reused=True)
                return archive
        raise RuntimeError('Collection path exists; outputs are preserved: '+str(target))
    summaries, errors = [], []
    tables = {name:[] for name in ('profiles_raw.csv','profiles_paired.csv','su3_times.csv','su3_windows.csv')}
    for entry in data['cases']:
        try:
            if not entry.get('array_job_id'):
                raise RuntimeError('Case has not been submitted')
            summary, rows = exports(entry)
            summaries.append(summary)
            for filename, items in rows.items():
                tables[filename].extend(items)
        except Exception as exc:
            errors.append(dict(index=entry['index'],configuration=entry['configuration'],
                               tau=entry['tau'],job=entry.get('array_job_id'),error=str(exc)))
    success = (len(summaries)==15 and not errors and
               summaries[0].get('n100_regression',{}).get('success') is True)
    if not success and not partial:
        raise RuntimeError('Incomplete/failed collection; results retained. Use --partial to '
                           'create a diagnostic archive. Errors: '+json.dumps(errors))
    if success and [len(tables[n]) for n in tables] != [9960,4980,120,30]:
        raise RuntimeError('Aggregate coverage differs from fifteen expected cases')
    target.mkdir()
    summary = dict(success=success, job_mapping=mapping, validated_cases=len(summaries),
                   expected_cases=15, cases=summaries, errors=errors,
                   scope='Structural/profile checks, original-field clover isolation, piter40 '
                         'SU3 control and n100 regression. piter4 sensitivity is reported '
                         'without a physics acceptance threshold.',
                   warning='Connected columns are per-configuration diagnostics. Form ensemble '
                           'averages of raw L, LS and S before the ratio/subtraction, with '
                           'correlated uncertainty analysis. Sites, axes and tau are not '
                           'independent ensemble samples.')
    q.c.write_json(target / 'SUMMARY.json', summary)
    for name, columns in [('profiles_raw.csv',q.checks.RAW_COLUMNS),
                          ('profiles_paired.csv',q.checks.PAIRED_COLUMNS),
                          ('su3_times.csv',q.checks.SU3_COLUMNS),
                          ('su3_windows.csv',q.checks.SU3_COLUMNS)]:
        q.checks.write_csv(target / name, columns, tables[name])
    report = ['# Multi-configuration projection comparison',
              f'Validated cases: {len(summaries)}/15; complete check success: {success}',
              '', 'piter4 acceptance remains a physics decision. Connected values are '
              'per-configuration diagnostics; compare absolute differences as well as '
              'relative changes when the connected signal is small.', '',
              '| cfg | tau | parent job | p4 max unitary defect (Frobenius) | '
              'p4 max determinant defect | max relative L real difference | '
              'max relative connected difference, combined axes |',
              '|---|---|---|---|---|---|---|']
    for s in summaries:
        report.append(f"| {s['configuration']} | {s['tau']} | {s['array_job_id']} | "
                      f"{float(s['p4_window']['unit_max']):.6e} | "
                      f"{float(s['p4_window']['det_max']):.6e} | "
                      f"{s['max_relative_difference']['L_RE']:.6e} | "
                      f"{s['max_relative_difference']['C_combined']:.6e} |")
    if errors:
        report += ['', 'Incomplete/failed cases:', *[json.dumps(e) for e in errors]]
    (target / 'REPORT.md').write_text('\n'.join(report)+'\n')
    q.c.write_json(target / 'JOB_MAPPING.json', mapping)
    include = [p for p in target.rglob('*') if p.is_file()]
    for entry in data['cases']:
        if not entry.get('array_job_id'):
            continue
        folder = q.output_directory(entry)
        include += [p for p in folder.rglob('*') if p.is_file()]
        include += [q.c.CASE / f"profiles_{entry['array_job_id']}_{entry['index']}.{ext}"
                    for ext in ('out','err')]
    include += [p for p in q.c.CASE.iterdir() if p.is_file()]
    include += [q.c.CASE / 'build' / name for name in ('manifest.sha256','modules.txt','compiler.txt','ldd.txt')]
    include += [p for p in q.HERE.iterdir() if p.is_file() and p.name != 'controller.lock']
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as bundle:
        for p in sorted(set(include)):
            if not p.is_file() or p.is_symlink():
                continue
            if p.suffix in ('.h5','.o','.a','.zip','.pyc') or p.stat().st_size > 8*1024**2:
                raise RuntimeError('Unexpected review file: '+str(p))
            bundle.write(p,p.relative_to(q.c.ROOT).as_posix())
    with zipfile.ZipFile(archive) as bundle:
        if bundle.testzip() is not None:
            raise RuntimeError('Review archive integrity failed')
    (target / 'archive.sha256').write_text(q.c.digest(archive)+'  '+archive.name+'\n')
    q.say('REVIEW_ARCHIVE', path=str(archive), success=success,
          raw_rows=len(tables['profiles_raw.csv']), paired_rows=len(tables['profiles_paired.csv']))
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partial',action='store_true',help='Export incomplete/failed cases explicitly as partial')
    args = parser.parse_args()
    with (q.HERE / 'controller.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Controller is active; it collects automatically after completion')
        q.verify_environment()
        collect(q.load_ledger(), partial=args.partial)
    if args.partial:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
