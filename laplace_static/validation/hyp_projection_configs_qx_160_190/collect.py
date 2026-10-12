#!/usr/bin/env python3
"""Require complete array coverage; retain partial/failed cases in review ZIP."""
import csv
import json
import sys
import zipfile
from pathlib import Path
import common as c
import check_task as check

def main():
    job=c.array_id(sys.argv[1])
    # Refuse to reinterpret results with modified package/build/preparation.
    c.prepared()
    target=c.CASE/f'collection_{job}'
    if target.exists():raise RuntimeError(f'Collection already exists: {target}')
    target.mkdir()
    summaries=[];errors=[];raw=[];paired=[];times=[];windows=[]
    for i in range(c.NCASES):
        selected=c.case_for(i)
        outdir=c.CASE/f'results_{job}_{i}'
        try:
            case=json.loads((outdir/'CASE.json').read_text())
            if case!=dict(selected,array_job_id=job):raise RuntimeError('Wrong case identity')
            for source in ('package.sha256','prepared.sha256','build/manifest.sha256'):
                if (outdir/source.replace('/','_')).read_bytes()!=(c.CASE/source).read_bytes():
                    raise RuntimeError('Measurement provenance changed')
            launches=json.loads((outdir/'LAUNCHES.json').read_text())
            expected=[(selected['configuration'],selected['tau'],p,f'p{p}') for p in (4,40)]
            if i==0:expected += [(100,4,p,f'n100_p{p}') for p in (4,40)]
            actual=[(r['configuration'],r['tau'],r['piter'],r['label']) for r in launches]
            if actual!=expected or any(r['exit_code']!=0 for r in launches):
                raise RuntimeError('MPI measurement inventory/status failed')
            # Recheck outputs rather than trusting a success marker alone.
            summary=check.check(outdir,selected['configuration'],selected['tau'],regression=(i==0))
            summary['index']=i;summary['mpi_seconds']=sum(r['seconds'] for r in launches)
            summaries.append(summary)
            for name,rows,count in [('profiles_raw.csv',raw,664),('profiles_paired.csv',paired,332),
                                    ('su3_times.csv',times,2*selected['tau']),('su3_windows.csv',windows,2)]:
                with (outdir/name).open(newline='') as stream:items=list(csv.DictReader(stream))
                if len(items)!=count:raise RuntimeError('Unexpected export count: '+name)
                rows.extend(items)
        except Exception as exc:
            errors.append(dict(index=i,configuration=selected['configuration'],tau=selected['tau'],error=str(exc)))
    success=(len(summaries)==c.NCASES and not errors and
             summaries[0].get('n100_regression',{}).get('success') is True)
    summary=dict(success=success,array_job_id=job,validated_cases=len(summaries),expected_cases=c.NCASES,
                 cases=summaries,errors=errors,
                 scope='Four configurations at fixed midpoint98, tau2/4/6, R1/2/4. '
                       'piter4 sensitivity diagnostics; no ensemble estimator or physics acceptance verdict.',
                 warning='Connected columns are per-configuration diagnostics. Spatial sites/axes/tau '
                         'are not independent ensemble samples. Form ensemble averages of raw L, LS, S '
                         'before the ratio/subtraction, with correlated uncertainty analysis.')
    c.write_json(target/'SUMMARY.json',summary)
    check.write_csv(target/'profiles_raw.csv',check.RAW_COLUMNS,raw)
    check.write_csv(target/'profiles_paired.csv',check.PAIRED_COLUMNS,paired)
    check.write_csv(target/'su3_times.csv',check.SU3_COLUMNS,times)
    check.write_csv(target/'su3_windows.csv',check.SU3_COLUMNS,windows)
    report=['# Multi-configuration projection comparison',f'Array job: {job}',
            f'Validated cases: {len(summaries)}/12; complete structural/control/regression success: {success}',
            '', 'piter4 acceptance is not decided by this checker. Inspect central-window maxima and '
            'paired observable differences together. Connected columns are single-configuration diagnostics.',
            '', '| cfg | tau | p4 max unitary defect (Frobenius) | p4 max determinant defect | '
            'max relative L real difference | max relative connected difference, combined axes |',
            '|---|---|---|---|---|---|']
    for s in summaries:
        report.append(f"| {s['configuration']} | {s['tau']} | {float(s['p4_window']['unit_max']):.6e} | "
                      f"{float(s['p4_window']['det_max']):.6e} | {s['max_relative_difference']['L_RE']:.6e} | "
                      f"{s['max_relative_difference']['C_combined']:.6e} |")
    if errors:report += ['','Incomplete/failed cases:',*[json.dumps(e) for e in errors]]
    (target/'REPORT.md').write_text('\n'.join(report)+'\n')
    archive=c.ROOT/f'hyp_projection_configs_160_190_review_{job}.zip'
    include=[]
    for folder in [outdir for outdir in sorted(c.CASE.glob(f'results_{job}_*'))]+[target]:
        include += [p for p in folder.rglob('*') if p.is_file()]
    include += [p for p in c.CASE.iterdir() if p.is_file()]
    include += [c.CASE/'build'/name for name in ('manifest.sha256','modules.txt','compiler.txt','ldd.txt')]
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as bundle:
        for p in sorted(set(include)):
            if not p.is_file() or p.is_symlink():continue
            if p.suffix in ('.h5','.o','.a','.zip') or p.stat().st_size>8*1024**2:
                raise RuntimeError(f'Unexpected review file: {p}')
            bundle.write(p,p.relative_to(c.ROOT).as_posix())
    with zipfile.ZipFile(archive) as bundle:
        if bundle.testzip() is not None:raise RuntimeError('Review ZIP integrity failure')
    print('REVIEW_ARCHIVE',archive)
    if not success:
        print('INCOMPLETE/FAILED: review retained; see SUMMARY.json and REPORT.md',file=sys.stderr)
        raise SystemExit(1)
    print('PASS complete coverage: 7968 raw rows, 3984 paired rows; inspect sensitivity before accepting piter4.')

if __name__=='__main__':main()
