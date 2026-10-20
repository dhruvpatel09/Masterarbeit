#!/usr/bin/env python3
"""Audit existing n100 modes, with n650 as a preserved numerical control."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRESERVED = ROOT/'laplace_static/analysis/qx_eigenpair_completion/inputs'
RESULTS = HERE/'results_n100'
TIMES = [0, 32, 64, 95, 96, 97, 98, 99, 100, 101, 128, 160, 191]
EXE = HERE/'build/test_saved_modes'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


legacy = module('qx_n100_legacy', ROOT/'laplace_static/validation/qx_projection_spectrum/campaign.py')
paths = module('qx_n100_paths', ROOT/'laplace_static/validation/hyp_projection_configs_qx/common.py')
parser = module('qx_n100_parser', PRESERVED/'job5925408/check_export.py')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def verify_sources():
    legacy.verify_sources()
    for name, expected in json.loads((HERE/'INPUTS.json').read_text())['files'].items():
        if digest(ROOT/name) != expected:
            raise ValueError('Audit source changed: '+name)


def verify_build():
    verify_sources()
    manifest = HERE/'build/manifest.sha256'
    for line in manifest.read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if digest(Path(name)) != expected:
            raise ValueError('Audit build changed: '+name)
    if not EXE.is_file():
        raise ValueError('Audit executable is missing')
    return digest(manifest)


def file_state(path):
    st = Path(path).stat()
    if st.st_size <= 0:
        raise ValueError('Empty input: '+str(path))
    return dict(bytes=st.st_size, mtime_ns=st.st_mtime_ns)


def eigenvalues(run, override=None):
    if override:
        selected = Path(override).resolve()
        file_state(selected)
        return selected
    candidates = sorted((Path(run)/'eigenvalues').glob('*.h5'))
    if len(candidates) != 1:
        raise ValueError('Expected one eigenvalue HDF5 file in %s; found %s. Use --eigenvalues for n100.'
                         % (Path(run)/'eigenvalues', [p.name for p in candidates]))
    file_state(candidates[0])
    return candidates[0]


def prepare(args):
    build_hash = verify_build()
    observable_hash = legacy.build_manifest()
    source = paths.N100_RUN
    n100_values = eigenvalues(source, args.eigenvalues)
    control = dict(legacy.choose_runs()[650])
    marker650 = Path(control['run'])/'COMPLETED.json'
    if digest(marker650) != control['completion_marker_sha256']:
        raise ValueError('n650 control completion marker changed')
    control_values = eigenvalues(control['run'])
    inputs = {}
    for cfg, run, values, times in [(100, source, n100_values, range(192)),
                                   (650, Path(control['run']), control_values, TIMES)]:
        selected = [legacy.GAUGES/('qX_p2n%d' % cfg), values]
        selected += [run/'eigenvectors'/('qX_p2n%d_evec_t%d.h5' % (cfg, t)) for t in times]
        for p in selected:
            inputs[str(p)] = file_state(p)
    inputs[str(marker650)] = file_state(marker650)
    RESULTS.mkdir(exist_ok=False)
    view = RESULTS/'input_run'
    view.mkdir()
    (view/'eigenvectors').symlink_to(source/'eigenvectors', target_is_directory=True)
    (view/'eigenvalues').symlink_to(n100_values.parent, target_is_directory=True)
    marker = view/'COMPLETED.json'
    write_json(marker, dict(configuration=100, timeslices=192, eigenvector_files=192,
                            eigenvector_bytes=sum(inputs[str(source/'eigenvectors'/('qX_p2n100_evec_t%d.h5' % t))]['bytes'] for t in range(192)),
                            eigenvalue_files=[n100_values.name], source_run=str(source),
                            scope='Inventory of pre-existing n100 files. This metadata view uses symlinks; '
                                  'no eigenvectors are generated or copied. Numerical quality is established '
                                  'separately by the successful audit SUMMARY.json.'))
    inputs[str(marker)] = file_state(marker)
    n100 = dict(configuration=100, run=str(view), source_run=str(source),
                gauge=str(legacy.GAUGES/'qX_p2n100'), eigenvectors=str(view/'eigenvectors'),
                eigenvalues=str(n100_values), completion_marker_sha256=digest(marker))
    control.update(gauge=str(legacy.GAUGES/'qX_p2n650'),
                   eigenvectors=str(Path(control['run'])/'eigenvectors'), eigenvalues=str(control_values))
    data = dict(configurations=[650, 100], Nv=100, times=TIMES, grid=[1, 4, 4, 4],
                runs=[control, n100], inputs=inputs, audit_build_sha256=build_hash,
                observable_build_sha256=observable_hash, source_manifest_sha256=digest(HERE/'INPUTS.json'),
                scope='Read-only saved-mode audit: n650 layout control, then n100. All 192 n100 file '
                      'names are present; 100 modes on 13 selected time slices are numerically checked.')
    write_json(RESULTS/'PREPARED.json', data)
    (RESULTS/'prepared.sha256').write_text(digest(RESULTS/'PREPARED.json')+'\n')
    print('Prepared n100 audit with n650 control: 64 ranks, 8 nodes, 13 selected slices each')


def prepared():
    build_hash = verify_build()
    if digest(RESULTS/'PREPARED.json') != (RESULTS/'prepared.sha256').read_text().strip():
        raise ValueError('Audit preparation changed')
    data = json.loads((RESULTS/'PREPARED.json').read_text())
    if (data['audit_build_sha256'] != build_hash or
            data['observable_build_sha256'] != legacy.build_manifest() or
            data['source_manifest_sha256'] != digest(HERE/'INPUTS.json')):
        raise ValueError('Prepared audit implementation changed')
    if data['configurations'] != [650, 100] or data['times'] != TIMES or data['Nv'] != 100:
        raise ValueError('Unexpected audit coverage')
    for name, expected in data['inputs'].items():
        if file_state(Path(name)) != expected:
            raise ValueError('Prepared audit input changed: '+name)
    return data


def submit(args):
    prepared()
    command = ['sbatch', '--parsable', str(HERE/'submit_micro.sh')]
    record = RESULTS/'SUBMISSION.json'
    write_json(record, dict(state='submitting', command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, check=False)
    record.write_text(json.dumps(dict(command=command, returncode=result.returncode,
                                     stdout=result.stdout, stderr=result.stderr), indent=2)+'\n')
    if result.returncode:
        raise ValueError('Audit submission rejected: '+result.stderr.strip())
    if not re.fullmatch(r'\d+(?:;\S+)?\n?', result.stdout):
        raise ValueError('Ambiguous audit submission; inspect SUBMISSION.json before retrying')
    print('SUBMITTED', result.stdout.strip(), 'n650 control and n100 saved-mode audit')


def control_regression(modes):
    with (PRESERVED/'job5925408/modes.csv').open(newline='') as stream:
        reference = [r for r in csv.DictReader(stream) if int(r['configuration']) == 650]
    if len(modes) != len(reference) or len(modes) != 1300:
        raise ValueError('n650 control inventory differs from preserved audit')
    for current, old in zip(modes, reference):
        if (current['t'], current['mode']) != (int(old['t']), int(old['mode'])):
            raise ValueError('n650 control coordinate mismatch')
        for name in parser.MODE_FIELDS:
            tolerance = 5e-12 if name == 'residual_relative' else 1e-12
            if not math.isclose(current[name], float(old[name]), rel_tol=5e-12, abs_tol=tolerance):
                raise ValueError('n650 layout control differs: '+name)


def reduce_outputs(data, launches):
    if [r['configuration'] for r in launches] != [650, 100] or any(r['mpi_exit_code'] for r in launches):
        raise ValueError('Audit launch inventory failed')
    modes, slices, runs = [], [], []
    for cfg in (650, 100):
        passed, current_modes, current_slices = parser.parse(RESULTS/('n%d.out' % cfg), cfg)
        if not passed:
            raise ValueError('Saved-mode numerical audit failed for n%d' % cfg)
        if cfg == 650:
            control_regression(current_modes)
        modes.extend(current_modes)
        slices.extend(current_slices)
        runs.append(dict(configuration=cfg, success=True, mpi_exit_code=0,
                         maximum_relative_residual=max(r['residual_relative'] for r in current_modes),
                         modes_above_requested_1e_minus9=sum(r['residual_relative'] > 1e-9 for r in current_modes)))
    summary = dict(success=True, runs=runs, mode_records=len(modes), slice_records=len(slices),
                   n650_layout_control_passed=True, residual_tolerance=1.1e-9,
                   requested_solver_res=1e-9, normalization_and_gram_tolerance=1e-9,
                   rayleigh_im_tolerance=1e-10, grid=[1, 4, 4, 4],
                   scope='Existing qcd Laplacian and preserved audit numerics on 64 ranks. '
                         'All 100 saved modes on 13 selected slices per configuration; full-time '
                         'eigenvalue metadata and spectra checked. n650 reproduces the preserved '
                         'control before n100 is admitted. This is not eigenpair generation.')
    return summary, modes, slices


def save_exports(summary, modes, slices):
    with tempfile.TemporaryDirectory() as name:
        temp = Path(name)
        for filename, rows, fields in [
                ('modes.csv', modes, ['configuration', 't', 'mode']+parser.MODE_FIELDS),
                ('slices.csv', slices, ['configuration', 't', 'status']+parser.SLICE_FIELDS)]:
            with (temp/filename).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
        write_json(temp/'SUMMARY.json', summary)
        for source in temp.iterdir():
            target = RESULTS/source.name
            if target.exists() and target.read_bytes() != source.read_bytes():
                raise ValueError('Existing audit export differs: '+source.name)
            target.write_bytes(source.read_bytes())


def run(args):
    data = prepared()
    if int(os.environ.get('SLURM_NTASKS', '0')) != 64 or int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != 8:
        raise ValueError('Run audit inside its 8-node, 64-rank allocation')
    launches = []
    try:
        for run in data['runs']:
            cfg = run['configuration']
            command = ['mpiexec', '-n', '64', str(EXE), str(cfg), run['gauge'], run['eigenvectors'], run['eigenvalues']]
            start = time.monotonic()
            with (RESULTS/('n%d.out' % cfg)).open('x') as stdout, (RESULTS/('n%d.err' % cfg)).open('x') as stderr:
                code = subprocess.run(command, stdout=stdout, stderr=stderr, check=False).returncode
            launches.append(dict(configuration=cfg, mpi_exit_code=code, elapsed_seconds=time.monotonic()-start, command=command))
            write_json(RESULTS/('launch_%d.json' % cfg), launches[-1])
            if code:
                raise ValueError('Saved-mode MPI audit failed for n%d' % cfg)
            passed, modes, _ = parser.parse(RESULTS/('n%d.out' % cfg), cfg)
            if not passed:
                raise ValueError('Numerical audit failed')
            if cfg == 650:
                control_regression(modes)
        write_json(RESULTS/'LAUNCHES.json', launches)
        save_exports(*reduce_outputs(data, launches))
    except Exception as exc:
        write_json(RESULTS/'FAILED.json', dict(success=False, error=str(exc), launches=launches))
        raise
    print('PASS n100 saved-mode audit and n650 layout control')


def collect(args):
    data = prepared()
    if (RESULTS/'FAILED.json').exists():
        raise ValueError('Audit has a failure record')
    launches = json.loads((RESULTS/'LAUNCHES.json').read_text())
    save_exports(*reduce_outputs(data, launches))
    archive = ROOT/'qx_saved_modes_n100_review.zip'
    files = [p for p in RESULTS.iterdir() if p.is_file()]
    files += [RESULTS/'input_run/COMPLETED.json', HERE/'build/manifest.sha256']
    files += [p for p in HERE.iterdir() if p.is_file() and p.suffix not in ('.out', '.err')]
    if archive.exists():
        with zipfile.ZipFile(archive) as z:
            expected = {p.relative_to(ROOT).as_posix(): p for p in files}
            if set(z.namelist()) != set(expected) or any(z.read(n) != expected[n].read_bytes() for n in z.namelist()):
                raise ValueError('Existing audit review differs; it was preserved')
    else:
        descriptor, name = tempfile.mkstemp(prefix='.qx_n100_review_', suffix='.tmp', dir=str(ROOT))
        os.close(descriptor)
        temporary = Path(name)
        try:
            with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as z:
                for p in sorted(files):
                    if p.is_symlink() or p.stat().st_size > 8*1024**2 or p.suffix == '.h5':
                        raise ValueError('Unexpected audit review file')
                    z.write(p, p.relative_to(ROOT).as_posix())
            with zipfile.ZipFile(temporary) as z:
                if z.testzip() is not None or len(z.namelist()) != len(files):
                    raise ValueError('Audit archive failed integrity verification')
            os.replace(str(temporary), str(archive))
        finally:
            if temporary.exists():
                temporary.unlink()
    print('PASS n100 audit review')
    print('REVIEW_ARCHIVE', archive)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest='action')
    for action in ('verify-sources', 'submit', 'run', 'collect'):
        sub.add_parser(action)
    p = sub.add_parser('prepare')
    p.add_argument('--eigenvalues', type=Path)
    args = ap.parse_args()
    if args.action is None:
        ap.error('Select verify-sources, prepare, submit, run or collect')
    if args.action == 'verify-sources':
        verify_sources()
        print('PASS audit source pins')
    else:
        globals()[args.action](args)


if __name__ == '__main__':
    main()
