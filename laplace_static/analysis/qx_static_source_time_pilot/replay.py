"""Cold native-log replay using unchanged pinned published checker sources."""
import importlib.util
import json
import math
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from common import HERE, CONFIGS, RADII, CUTOFFS, digest, read


def close_summary(actual, expected, path='summary'):
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError('Replay summary fields: ' + path)
        for key in expected:
            close_summary(actual[key], expected[key], path + '.' + key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError('Replay summary length: ' + path)
        for i, (a, b) in enumerate(zip(actual, expected)):
            close_summary(a, b, path + '[%d]' % i)
    elif path.endswith('.minimum_equal_time_relative_pivot'):
        if not (math.isfinite(actual) and actual >= 0 and (actual == 0) == (expected == 0)
                and math.isclose(actual, expected, rel_tol=5e-12,
                                 abs_tol=32 * sys.float_info.epsilon)):
            raise ValueError('Equal-time pivot replay: ' + path)
    elif actual != expected:
        raise ValueError('Native replay mismatch: ' + path)


def check_training(root, sources, reference):
    primary = root / 'laplace_static/analysis/qx_static_matrix_71configs/inputs/primary_matrices.npz'
    cache = np.load(str(primary), allow_pickle=False)
    if (list(cache['configurations']) != list(range(100, 801, 10))
            or list(cache['radii']) != list(RADII) or list(cache['cutoffs']) != list(CUTOFFS)
            or cache['C'].shape != (71, 4, 33, 4, 4)):
        raise ValueError('Published primary cache geometry')
    labels = list(map(int, cache['configurations']))
    training = [c for c in labels if c not in CONFIGS]
    if (training != sources['training_configuration_labels'] or len(training) != 68
            or sources['excluded_pilot_configuration_labels'] != list(CONFIGS)
            or sources['training_t0'] != 90):
        raise ValueError('Training/pilot leakage or changed source scope')
    C = cache['C']
    for (cfg, source, sink, tau, radius), row in reference.items():
        z = C[labels.index(cfg), RADII.index(radius), tau, CUTOFFS.index(sink), CUTOFFS.index(source)]
        if float(row['L_RE']) != z.real or float(row['L_IM']) != z.imag:
            raise ValueError('Baseline reference differs from published primary cache')
    location = root / 'laplace_static/analysis/qx_static_matrix_71configs'
    old_path = sys.path[:]
    old_stats = sys.modules.pop('statistics_tools', None)
    sys.path.insert(0, str(location))
    try:
        spec = importlib.util.spec_from_file_location('source_time_published_training', str(location / 'projections.py'))
        algorithm = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(algorithm)
        seen = set()
        for vector in sources['source_candidates']:
            radius, basis, anchor = vector['R'], vector['basis'], vector['anchor']
            key = radius, basis, anchor
            if key in seen:
                raise ValueError('Duplicate frozen source')
            seen.add(key)
            indices = [CUTOFFS.index(n) for n in vector['Nv']]
            if indices != algorithm.BASES[basis]:
                raise ValueError('Changed frozen source basis')
            mean = algorithm.submatrix(C[np.array([c in training for c in labels]), RADII.index(radius)], indices).mean(0)
            fresh, _ = algorithm.train(mean, anchor)
            frozen = np.array(vector['coefficients_re']) + 1j * np.array(vector['coefficients_im'])
            norm = np.vdot(frozen, mean[0].dot(frozen))
            overlap = abs(np.vdot(frozen, mean[0].dot(fresh)))
            H = (mean[anchor] + mean[anchor].conj().T) / 2
            eigenvalue = float(np.vdot(fresh, H.dot(fresh)).real)
            residual = np.linalg.norm(H.dot(frozen) - eigenvalue * mean[0].dot(frozen))
            scale = (np.linalg.norm(H) + abs(eigenvalue) * np.linalg.norm(mean[0])) * np.linalg.norm(frozen)
            if abs(norm - 1) > 1e-9 or abs(overlap - 1) > 1e-8 or residual / scale > 1e-9:
                raise ValueError('Frozen source does not reproduce the excluded-pilot training solution')
        expected = {(r, b, a) for r in RADII for b in ('pair', 'triple', 'full') for a in (4, 6, 8)}
        if seen != expected:
            raise ValueError('Frozen source inventory')
    finally:
        sys.path[:] = old_path
        sys.modules.pop('statistics_tools', None)
        if old_stats is not None:
            sys.modules['statistics_tools'] = old_stats
    return dict(training_configurations=68, excluded_pilot_configurations=list(CONFIGS),
                frozen_candidates=36, baseline_entries_exactly_matched=6336,
                phase_invariant_training_replay=True, sources_retuned_for_analysis=False)


def cold_replay(progress=False):
    inputs = read(HERE / 'inputs/SOURCES.json')
    archive = HERE / 'inputs/job6036903_review.zip'
    if archive.stat().st_size != inputs['archive']['bytes'] or digest(archive) != inputs['archive']['sha256']:
        raise ValueError('Changed native review archive')
    before = dict(sys.modules)
    old_path = sys.path[:]
    with tempfile.TemporaryDirectory(prefix='qx_source_time_replay_') as directory:
        root = Path(directory)
        try:
            with zipfile.ZipFile(str(archive)) as z:
                if len(z.namelist()) != 100 or len(set(z.namelist())) != 100:
                    raise ValueError('Native archive inventory')
                for info in z.infolist():
                    p = Path(info.filename)
                    if p.is_absolute() or '..' in p.parts or (info.external_attr >> 16) & 0o170000 == 0o120000:
                        raise ValueError('Unsafe native archive member')
                z.extractall(str(root))
            for name, record in inputs['additional_inputs'].items():
                source = HERE / 'inputs/published' / name
                if source.stat().st_size != record['bytes'] or digest(source) != record['sha256']:
                    raise ValueError('Changed replay dependency: ' + name)
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(str(source), str(target))
            native = root / 'laplace_static/validation/qx_static_source_time_native'
            design = root / 'laplace_static/validation/qx_static_source_time'
            for parent in (native, design):
                manifest = read(parent / 'FILE_MANIFEST.json')
                for name, record in manifest['files'].items():
                    path = parent / name
                    if path.stat().st_size != record['bytes'] or digest(path) != record['sha256']:
                        raise ValueError('Changed archived source: ' + name)
                for name, record in manifest['dependencies'].items():
                    path = root / name
                    if path.exists() and (path.stat().st_size != record['bytes'] or digest(path) != record['sha256']):
                        raise ValueError('Published helper differs from native pin: ' + name)
            # Generic names are private to this temporary replay environment.
            for name in ('campaign', 'parents', 'pilot_plan', 'receipts', 'check_outputs', 'fixed_projections'):
                sys.modules.pop(name, None)
            sys.path.insert(0, str(native))
            spec = importlib.util.spec_from_file_location('source_time_native_replay', str(native / 'campaign.py'))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            target = native / 'campaigns/pilot'
            prepared = read(target / 'PREPARED.json')
            archived = read(target / 'SUMMARY.json')
            if digest(target / 'PREPARED.json') != (target / 'prepared.sha256').read_text().strip():
                raise ValueError('Prepared campaign checksum')
            if (not archived['success'] or archived['validated_configurations'] != list(CONFIGS)
                    or archived['missing_configurations'] or [e['configuration'] for e in prepared['cases']] != list(CONFIGS)):
                raise ValueError('Incomplete pilot cohort')
            build = read(native / 'build/MANIFEST.json')
            if digest(native / 'build/MANIFEST.json') != prepared['build_manifest_sha256']:
                raise ValueError('Prepared build receipt changed')
            available_build, unavailable_build = [], []
            for name, expected in build['files'].items():
                if (root / name).exists():
                    if digest(root / name) != expected:
                        raise ValueError('Archived build input changed: ' + name)
                    available_build.append(name)
                else:
                    unavailable_build.append(name)
            reference = module.check.load_references(module.REFERENCE)
            sources = read(module.SOURCES)
            if (digest(module.SOURCES) != prepared['source_vectors_sha256']
                    or digest(module.REFERENCE) != prepared['reference_sha256']
                    or digest(design / 'PLAN.json') != prepared['design_plan_sha256']):
                raise ValueError('Prepared frozen source/reference/design changed')
            training = check_training(root, sources, reference)
            values, blocks, cases, launches = [], [], [], []
            for entry in prepared['cases']:
                out = module.complete_case(target, entry)
                if out is None:
                    raise ValueError('Missing native completion marker')
                rows, parts, case = module.replay_case(out, entry, reference)
                module.csv_equal(out / 'matrix_correlators.csv', rows)
                module.csv_equal(out / 'rank_blocks.csv', parts)
                close_summary(case, next(c for c in archived['cases'] if c['configuration'] == entry['configuration']))
                values += rows
                blocks += parts
                cases.append(case)
                launches += read(out / 'LAUNCHES.json')
                if progress:
                    print('Replayed native case n%d: %d scans' % (entry['configuration'], len(case['scans'])), flush=True)
            module.csv_equal(target / 'matrix_correlators.csv', values)
            module.csv_equal(target / 'rank_blocks.csv', blocks)
            projected, contrasts = module.derive(values, reference, sources)
            module.csv_equal(target / 'fixed_source_correlators.csv', projected)
            module.csv_equal(target / 'paired_fixed_source_contrasts.csv', contrasts)
            if (len(values), len(blocks), len(projected), len(contrasts)) != (67584, 147456, 13068, 9504):
                raise ValueError('Changed native export inventory')
            # Independent dense contraction catches an incorrect matrix orientation/conjugation.
            matrices = {}
            for row in values:
                if row['axis'] == '-1':
                    key = tuple(int(row[k]) for k in ('cfg', 't0', 'piter', 'R', 'tau'))
                    matrices.setdefault(key, np.zeros((4, 4), dtype=complex))[CUTOFFS.index(int(row['Nv_sink'])), CUTOFFS.index(int(row['Nv_src']))] = complex(float(row['L_RE']), float(row['L_IM']))
            for (cfg, src, sink, tau, radius), row in reference.items():
                matrices.setdefault((cfg, 90, 10, radius, tau), np.zeros((4, 4), dtype=complex))[CUTOFFS.index(sink), CUTOFFS.index(src)] = complex(float(row['L_RE']), float(row['L_IM']))
            vectors = {(v['R'], v['basis'], v['anchor']): v for v in sources['source_candidates']}
            for row in projected:
                vector = vectors[row['R'], row['basis'], row['anchor']]
                v = np.zeros(4, dtype=complex)
                for n, re, im in zip(vector['Nv'], vector['coefficients_re'], vector['coefficients_im']):
                    v[CUTOFFS.index(n)] = complex(re, im)
                C = matrices[row['configuration'], row['t0'], row['piter'], row['R'], row['tau']]
                expected = complex(row['Q_RE'], row['Q_IM'])
                terms = v.conj()[:, None] * C * v[None, :]
                bound = 128 * np.finfo(float).eps * np.abs(terms).sum() + 1e-30
                if abs(np.vdot(v, C.dot(v)) - expected) > bound:
                    raise ValueError('Independent frozen-source contraction')
            report = dict(success=True, configurations=list(CONFIGS), cases=cases,
                          matrix_rows=len(values), rank_block_rows=len(blocks), fixed_source_rows=len(projected),
                          paired_fixed_source_rows=len(contrasts), training=training,
                          archived_build_files_verified=available_build,
                          build_inputs_not_contained_in_review=unavailable_build,
                          replay_scope='Archived stdout, source pins, receipts, reductions and projections; no MPI execution or physical input reopening.')
            return dict(values=values, reference=reference, projected=projected, contrasts=contrasts,
                        sources=sources, cases=cases, launches=launches, report=report)
        finally:
            sys.path[:] = old_path
            for name in list(sys.modules):
                if name not in before:
                    del sys.modules[name]
            for name in ('campaign', 'parents', 'pilot_plan', 'receipts', 'check_outputs', 'fixed_projections'):
                if name in before:
                    sys.modules[name] = before[name]
