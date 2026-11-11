#!/usr/bin/env python3
"""Inspect physical input availability; never build, run MPI or submit jobs."""
import argparse,hashlib,json,os
from pathlib import Path
from plan import HERE,ROOT,BASE_COMMIT,CONFIGS,PILOT,build_plan,audit_times


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def fingerprint(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size<=0:raise ValueError('Missing or empty input: '+str(path))
    s=path.stat()
    return dict(path=str(path),resolved=str(path.resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def inspect_inputs(runs,gauge_root):
    if not isinstance(runs,list) or [r['configuration'] for r in runs]!=list(CONFIGS):raise ValueError('RUNS must contain all 71 ordered configuration labels exactly once')
    report=dict(base_commit=BASE_COMMIT,stage='filesystem_inventory_only',success=False,cases=[],
        native_mode_audit_completed=False,native_extension_built=False,jobs_submitted=False,
        interpretation='File availability and admission receipts only; no HDF5 residual, orthogonality, transporter or native contraction validation has been performed.')
    for cfg in PILOT:
        record=dict(configuration=cfg,available=False,required_times=audit_times(cfg),errors=[],vectors=[])
        try:
            source=next(r for r in runs if r['configuration']==cfg)
            if cfg==640:
                eigenvalues=Path(source['eigenvalues']);vectors=Path(source['eigenvectors'])
                record['admission']=dict(kind='recovered_n640_eigenvalue_checksum',expected_sha256=source['eigenvalue_sha256'])
                if digest(eigenvalues)!=source['eigenvalue_sha256']:raise ValueError('Recovered n640 eigenvalue admission receipt changed')
            else:
                run=Path(source['run']);marker=run/'COMPLETED.json'
                if digest(marker)!=source['completion_marker_sha256']:raise ValueError('Saved-mode completion marker changed')
                values=sorted((run/'eigenvalues').glob('*.h5'))
                if len(values)!=1:raise ValueError('Expected one admitted eigenvalue table')
                eigenvalues=values[0];vectors=run/'eigenvectors'
                record['admission']=dict(kind='saved_mode_completion_marker',marker=str(marker),sha256=digest(marker))
            record['gauge']=fingerprint(Path(gauge_root)/('qX_p2n%d'%cfg))
            record['eigenvalues']=dict(fingerprint(eigenvalues),sha256=digest(eigenvalues))
            for t in record['required_times']:
                path=vectors/('qX_p2n%d_evec_t%d.h5'%(cfg,t))
                try:record['vectors'].append(dict(fingerprint(path),time=t))
                except (OSError,ValueError) as exc:record['errors'].append(str(exc))
            record['available']=not record['errors']
        except (OSError,ValueError,KeyError,StopIteration) as exc:record['errors'].append(str(exc))
        report['cases'].append(record)
    report['success']=all(c['available'] for c in report['cases'])
    report['available_configurations']=[c['configuration'] for c in report['cases'] if c['available']]
    report['inspected_mode_files']=sum(len(c['vectors']) for c in report['cases'])
    report['estimated_total_mode_bytes']=sum(v['bytes'] for c in report['cases'] for v in c['vectors'])
    return report


def save_report(path,value):
    path=Path(path).resolve()
    if ROOT not in path.parents or (HERE/'work').resolve() not in path.parents:
        raise ValueError('Keep the preflight report within this package work directory')
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_name(path.name+'.pending')
    with tmp.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    tmp.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'work/PREFLIGHT.json')
    parser.add_argument('--gauge-root',type=Path,default=Path('/hppfs/scratch/00/di98qes3/qX'))
    args=parser.parse_args()
    from verify_files import verify
    verify()
    runs=json.loads((ROOT/'laplace_static/validation/qx_static_extension/inputs/RUNS.json').read_text())
    report=inspect_inputs(runs,args.gauge_root);save_report(args.output,report)
    for case in report['cases']:
        print('n{}: {} ({}/{} mode files)'.format(case['configuration'],'AVAILABLE' if case['available'] else 'BLOCKED',len(case['vectors']),len(case['required_times'])))
        for error in case['errors']:print('  '+error)
    print('PREFLIGHT_REPORT',args.output.resolve())
    print('No build, MPI launch or Slurm submission performed. A successful inventory is not native validation.')
    if not report['success']:raise SystemExit(1)


if __name__=='__main__':main()
