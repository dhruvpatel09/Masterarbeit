"""Cold replay of native stdout, receipts, references and every rank reduction.

No MPI execution, physical gauge/eigenvector access or native-binary execution.
Captured checker code is admitted against independent published source pins.
"""
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import math
import shutil
import stat
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
import numpy as np
from import_review import HERE, check_file, digest

ROOT = HERE.parents[2]
CONFIGS = list(range(100, 801, 10))
RADII = [4, 12, 24, 32]
CUTOFFS = [40, 60, 80, 100]
NATIVE = 'laplace_static/validation/qx_static_matrix'
FULL = 'laplace_static/validation/qx_static_matrix_full71'
BUILD_HASH = '4138399da51637487d768ea10778557aab6a4c04b7106ae9a633fd3648db33e3'
WORKFLOW_HASH = 'e6ae1c16785d841a91747beb40517c76a6217cefd57e51ad5f41bc2705e328c7'
PREPARED_HASH = '6457c812002be0740cf178dbd49bdcf014e01f5a5fa4445abb304847c5d00801'


def read(path):
    return json.loads(Path(path).read_text())


def safe_name(name):
    p = PurePosixPath(name)
    if not name or p.is_absolute() or '..' in p.parts or '\\' in name or name != p.as_posix():
        raise ValueError('Unsafe archive path: ' + name)
    return p


@contextlib.contextmanager
def extracted_review(here=HERE):
    source = read(here/'inputs/SOURCES.json')
    (here/'work').mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='replay_', dir=str(here/'work')) as tmp:
        tmp = Path(tmp)
        archive = tmp/'review.zip'
        whole = hashlib.sha256()
        with archive.open('wb') as out:
            for receipt in source['chunks']:
                p = here/receipt['path']
                check_file(p, receipt)
                with p.open('rb') as stream:
                    for block in iter(lambda: stream.read(1048576), b''):
                        out.write(block)
                        whole.update(block)
        if archive.stat().st_size != source['archive_bytes'] or whole.hexdigest() != source['archive_sha256']:
            raise ValueError('Reassembled archive mismatch')
        root = tmp/'source'
        with zipfile.ZipFile(str(archive)) as z:
            names = z.namelist()
            if len(names) != source['archive_members'] or len(set(names)) != len(names):
                raise ValueError('Archive member inventory differs')
            if sum(i.file_size for i in z.infolist()) != source['uncompressed_bytes']:
                raise ValueError('Archive uncompressed byte count differs')
            for info in z.infolist():
                safe_name(info.filename)
                if stat.S_ISLNK(info.external_attr >> 16) or info.is_dir():
                    raise ValueError('Unexpected linked/directory member')
            inventory = json.loads(z.read(source['inventory_path']).decode('utf-8'))
            expected = {r['path']: r for r in inventory['files']}
            if len(expected) != len(inventory['files']) or set(names) != set(expected) | {source['inventory_path']}:
                raise ValueError('Archive receipt member set differs')
            if inventory['omitted']:
                raise ValueError('Native campaign evidence was omitted')
            for name in names:
                target = root/name
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(name) as inp, target.open('wb') as out:
                    h = hashlib.sha256()
                    size = 0
                    for block in iter(lambda: inp.read(1048576), b''):
                        h.update(block)
                        size += len(block)
                        out.write(block)
                if name in expected and (size != expected[name]['bytes'] or h.hexdigest() != expected[name]['sha256']):
                    raise ValueError('Archived member mismatch: ' + name)
        yield root, inventory


def checker(root, here=HERE):
    for name, expected in read(here/'inputs/CHECKER_PINS.json').items():
        if digest(root/name) != expected:
            raise ValueError('Captured checker differs from published code: ' + name)
    names = ('matrix_plan', 'check_outputs', 'measurement_plan')
    previous = sys.path[:]
    modules = {name: sys.modules.get(name) for name in names}
    for name in names:
        sys.modules.pop(name, None)
    sys.path.insert(0, str(root/NATIVE))
    try:
        spec = importlib.util.spec_from_file_location('qx_matrix_archived_native', str(root/NATIVE/'campaign.py'))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path[:] = previous
        for name in names:
            sys.modules.pop(name, None)
            if modules[name] is not None:
                sys.modules[name] = modules[name]


def checked_csv(path, rows):
    with Path(path).open(newline='') as stream:
        if list(csv.DictReader(stream)) != rows:
            raise ValueError('CSV differs from native stdout: ' + str(path))


def compare_case_summary(archived, recomputed, configuration, path=()):
    """Exact native metadata, with roundoff allowance for one derived field.

    The native checker has already rerun all Gram/positivity checks. Only its
    dimensionless minimum Cholesky pivot is recalculated through platform math
    routines; counts, audit values and every other field remain exact.
    """
    label = 'n{} summary'.format(configuration)
    for part in path:
        label += '[{}]'.format(part) if isinstance(part, int) else '.'+part
    if type(archived) is not type(recomputed):
        raise ValueError('{}: type mismatch (archived={!r}, replay={!r})'.format(label, archived, recomputed))
    if isinstance(recomputed, dict):
        if set(archived) != set(recomputed):
            raise ValueError(label+': field inventory mismatch')
        for name in recomputed:
            compare_case_summary(archived[name], recomputed[name], configuration, path+(name,))
    elif isinstance(recomputed, list):
        if len(archived) != len(recomputed):
            raise ValueError(label+': list length mismatch')
        for i, (a, b) in enumerate(zip(archived, recomputed)):
            compare_case_summary(a, b, configuration, path+(i,))
    elif isinstance(recomputed, float):
        if not math.isfinite(archived) or not math.isfinite(recomputed):
            raise ValueError(label+': nonfinite diagnostic')
        derived_pivot = (len(path) == 3 and path[0] == 'scans'
                         and isinstance(path[1], int)
                         and path[2] == 'minimum_equal_time_relative_pivot')
        bound = 0.
        if derived_pivot:
            if min(archived, recomputed) < 0 or (archived > 0) != (recomputed > 0):
                raise ValueError(label+': pivot positivity class changed')
            bound = 32*sys.float_info.epsilon + 5e-12*max(abs(archived), abs(recomputed))
        if abs(archived-recomputed) > bound:
            raise ValueError('{}: archived={:.17g}, replay={:.17g}, difference={:.3g}, bound={:.3g}'.format(
                label, archived, recomputed, abs(archived-recomputed), bound))
    elif archived != recomputed:
        raise ValueError('{}: archived={!r}, replay={!r}'.format(label, archived, recomputed))


class CsvHash:
    def __init__(self):
        self.hash = hashlib.sha256()
        self.fields = None
        self.count = 0

    def add(self, rows):
        stream = io.StringIO(newline='')
        fields = list(rows[0])
        if self.fields is not None and fields != self.fields:
            raise ValueError('Changed CSV field order')
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        if self.fields is None:
            self.fields = fields
            writer.writeheader()
        writer.writerows(rows)
        self.hash.update(stream.getvalue().encode('utf-8'))
        self.count += len(rows)


def raw_replay(here=HERE, progress=False):
    source = read(here/'inputs/SOURCES.json')
    for name, receipt in source['dependencies'].items():
        check_file(ROOT/name, receipt)
    with extracted_review(here) as (root, inventory):
        base = checker(root, here)
        build = root/NATIVE/'build/MANIFEST.json'
        workflow = root/FULL/'FILE_MANIFEST.json'
        target = root/FULL/'campaigns/full71'
        pilot = root/NATIVE/'campaigns/pilot'
        if digest(build) != BUILD_HASH or digest(workflow) != WORKFLOW_HASH or digest(target/'PREPARED.json') != PREPARED_HASH:
            raise ValueError('Native build/workflow/preparation pins differ')
        omitted = {r['path']: r['sha256'] for r in inventory['omitted_native_inputs']}
        manifest = read(build)
        if manifest['inputs_sha256'] != digest(root/NATIVE/'INPUTS.json'):
            raise ValueError('Native build input specification differs')
        for name, sha in manifest['files'].items():
            actual = digest(root/name) if (root/name).is_file() else omitted.pop(name, None)
            if actual != sha:
                raise ValueError('Native compiled input receipt differs: ' + name)
        if omitted:
            raise ValueError('Unbound omitted native inputs')
        for group, folder in (('files', root/FULL), ('dependencies', root)):
            for name, receipt in read(workflow)[group].items():
                path = folder/name
                if not path.exists() and group == 'dependencies':
                    path = ROOT/name
                check_file(path, receipt)
        plan, pp, summary = read(target/'PREPARED.json'), read(pilot/'PREPARED.json'), read(target/'SUMMARY.json')
        new = [c for c in CONFIGS if c not in (170, 640, 650)]
        geom = source['geometry']
        if (plan['target_configurations'] != CONFIGS or plan['reused_configurations'] != [170, 640, 650]
                or plan['geometry'] != geom or plan['piter'] != 10 or plan['cases_per_task'] != 2
                or plan['slots'] != [list(range(i, i+2)) for i in range(0, 68, 2)]
                or [e['configuration'] for e in plan['cases']] != new
                or [e['index'] for e in plan['cases']] != list(range(68))
                or any(e['piters'] != [10] or e['input']['configuration'] != e['configuration'] for e in plan['cases'])):
            raise ValueError('Full71 measurement plan differs')
        if (pp['geometry'] != geom or [e['configuration'] for e in pp['cases']] != [170, 640, 650]
                or pp['build_manifest_sha256'] != BUILD_HASH or plan['native_build_manifest_sha256'] != BUILD_HASH
                or plan['workflow_manifest_sha256'] != WORKFLOW_HASH
                or plan['pilot']['prepared_sha256'] != digest(pilot/'PREPARED.json')
                or (target/'prepared.sha256').read_text().strip() != PREPARED_HASH
                or (pilot/'prepared.sha256').read_text().strip() != digest(pilot/'PREPARED.json')):
            raise ValueError('Pilot reuse or preparation receipts differ')
        if not summary['success'] or summary['validated_configurations'] != CONFIGS or summary['missing_configurations']:
            raise ValueError('Incomplete native full71 summary')
        reference = ROOT/'laplace_static/analysis/qx_static_extension_71configs/inputs/job5958012/campaigns/late_modes_full71/correlators.csv'
        if plan['reference_sha256'] != digest(reference):
            raise ValueError('Published diagonal reference differs')
        references = base.check.load_references(reference, CONFIGS)
        entries = {e['configuration']: (target, e) for e in plan['cases']}
        entries.update({e['configuration']: (pilot, e) for e in pp['cases']})
        C = np.full((71, 4, 33, 4, 4), np.nan + 1j*np.nan, dtype=np.complex128)
        text = np.full((71, 4, 33, 4), '', dtype='U32')
        matrix_hash, rank_hash = CsvHash(), CsvHash()
        checks, launches, reuse = [], [], []
        for ci, cfg in enumerate(CONFIGS):
            parent, entry = entries[cfg]
            if entry['piters'] != ([10, 40] if cfg == 650 else [10]):
                raise ValueError('Changed projection control inventory')
            out = base.complete_case(parent, entry)
            if out is None:
                raise ValueError('Missing complete case')
            rows, blocks, result = base.replay_case(out, entry, references)
            checked_csv(out/'matrix_correlators.csv', rows)
            checked_csv(out/'rank_blocks.csv', blocks)
            recorded = read(out/'SUMMARY.json')
            if not set(result).issubset(recorded):
                raise ValueError('n{} summary: missing native validation field'.format(cfg))
            compare_case_summary({k: recorded[k] for k in result}, result, cfg)
            if parent == pilot:
                reuse.append(dict(configuration=cfg, attempt=out.relative_to(root).as_posix(),
                    summary_sha256=digest(out/'SUMMARY.json'), case_sha256=digest(out/'CASE.json'),
                    launches_sha256=digest(out/'LAUNCHES.json')))
            selected = [r for r in rows if r['piter'] == '10']
            selected_blocks = [r for r in blocks if r['piter'] == '10']
            matrix_hash.add(selected)
            rank_hash.add(selected_blocks)
            for row in selected:
                if row['axis'] != '-1':
                    continue
                ri, t, sink, src = RADII.index(int(row['R'])), int(row['tau']), CUTOFFS.index(int(row['Nv_sink'])), CUTOFFS.index(int(row['Nv_src']))
                if np.isfinite(C[ci, ri, t, sink, src]):
                    raise ValueError('Duplicate primary matrix element')
                C[ci, ri, t, sink, src] = complex(float(row['L_RE']), float(row['L_IM']))
                if sink == src:
                    text[ci, ri, t, sink] = row['L_RE']
            checks.append(result)
            launches.extend(json.loads((out/'LAUNCHES.json').read_text()))
            if progress and (ci % 10 == 0 or ci == 70):
                print('Replayed {}/71 configurations'.format(ci+1), flush=True)
        if reuse != plan['pilot']['cases'] or not np.isfinite(C).all() or (text == '').any():
            raise ValueError('Incomplete primary array or pilot reuse receipt')
        for name, accumulated, expected_count in (('matrix_correlators.csv', matrix_hash, 599808), ('rank_blocks.csv', rank_hash, 872448)):
            count_key = 'matrix_rows' if name == 'matrix_correlators.csv' else 'rank_block_rows'
            if accumulated.count != expected_count or summary[count_key] != expected_count:
                raise ValueError('Combined numerical row count differs')
            if accumulated.hash.hexdigest() != summary['outputs_sha256'][name] or digest(target/name) != accumulated.hash.hexdigest():
                raise ValueError('Combined numerical CSV differs from complete stdout replay')
        scans = [s for r in checks for s in r['scans'] if s['piter'] == 10]
        audit = [r['audit'] for r in checks]
        scan_times = [r['seconds'] for r in launches if r['kind'] == 'scan' and r['piter'] == 10]
        report = dict(configuration_labels=CONFIGS, matrix_rows=matrix_hash.count, rank_block_rows=rank_hash.count,
            audit_slices=sum(r['slices'] for r in audit), audited_modes=sum(r['modes'] for r in audit), launches=len(launches),
            maximum_saved_mode_relative_residual=max(r['maximum_relative_residual'] for r in audit),
            exact_p10_diagonal_reference_rows=sum(r['diagonal_reference_rows'] for r in scans),
            diagonal_maximum_absolute_difference=max(r['diagonal_reference_maximum_absolute_difference'] for r in scans),
            direct_contraction_checks=sum(r['direct_tests'] for r in scans), legacy_diagonal_checks=sum(r['legacy_tests'] for r in scans),
            equal_time_gram_checks=sum(r['equal_time_gram_checks'] for r in scans), periodic_roundtrips=sum(r['periodic_roundtrips'] for r in scans),
            maximum_unitarity_defect=max(r['unit_max'] for r in scans), maximum_determinant_defect=max(r['det_max'] for r in scans),
            n650_p10_p40_maximum_component_difference=next(r['p10_p40_maximum_component_difference'] for r in checks if r['configuration'] == 650),
            mean_p10_scan_seconds=float(np.mean(scan_times)), minimum_p10_scan_seconds=min(scan_times), maximum_p10_scan_seconds=max(scan_times),
            archive_members=source['archive_members'], archive_sha256=source['archive_sha256'],
            omitted_native_binaries='Four binary/library hashes checked against native build manifest; binaries not executed.',
            matrix_orientation='configuration, R, tau, sink row, source column; axis=-1 already averaged over x,y,z')
        return C, text, report


if __name__ == '__main__':
    C, text, report = raw_replay(progress=True)
    np.savez_compressed(str(HERE/'inputs/primary_matrices.npz'), C=C, diagonal_text=text,
                        configurations=CONFIGS, radii=RADII, cutoffs=CUTOFFS)
    (HERE/'results/VALIDATION.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print('PASS full native replay:', report['matrix_rows'], 'matrix rows;', report['exact_p10_diagonal_reference_rows'], 'exact diagonal regressions.')
