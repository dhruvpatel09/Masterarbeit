"""Paths, provenance and immutable run selection for this comparison."""
import hashlib
import json
import re
from pathlib import Path

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[2]
GAUGE_ROOT = Path('/hppfs/scratch/00/di98qes3/qX')
N100_RUN = Path('/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/'
                'n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4')
CONFIGS = (200, 210, 220, 230)
TAUS = (2, 4, 6)
NCASES = len(CONFIGS) * len(TAUS)
AUDIT = ROOT / 'mental/validation/qx_saved_modes_200_230/results_5921385/SUMMARY.json'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def verify_manifest(path):
    for line in Path(path).read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        target = Path(name)
        if not target.is_absolute():
            target = ROOT / target
        if digest(target) != expected:
            raise RuntimeError(f'Checksum mismatch: {target}')

def write_json(path, data):
    with Path(path).open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')

def check_file(path):
    path = Path(path)
    with path.open('rb') as stream:
        if not stream.read(1):
            raise RuntimeError(f'Empty input: {path}')
    status = path.stat()
    return dict(path=str(path), bytes=status.st_size, mtime_ns=status.st_mtime_ns)

def case_for(index):
    if not 0 <= index < NCASES:
        raise RuntimeError('Array index must be 0..11')
    cfg, tau = CONFIGS[index // 3], TAUS[index % 3]
    return dict(index=index, configuration=cfg, tau=tau,
                t0=98-tau//2, t1=98+tau//2, tmid=98)

def prepared():
    verify_manifest(CASE / 'package.sha256')
    verify_manifest(CASE / 'inputs.sha256')
    verify_manifest(CASE / 'build/manifest.sha256')
    verify_manifest(CASE / 'prepared.sha256')
    data = json.loads((CASE / 'PREPARED.json').read_text())
    if data['cases'] != [case_for(i) for i in range(NCASES)]:
        raise RuntimeError('Unexpected prepared cases')
    for run in data['runs']:
        marker = Path(run['run']) / 'COMPLETED.json'
        if digest(marker) != run['completion_marker_sha256']:
            raise RuntimeError(f'Completion marker changed: {marker}')
    for item in data['inputs']:
        status = Path(item['path']).stat()
        if (status.st_size,status.st_mtime_ns) != (item['bytes'],item['mtime_ns']):
            raise RuntimeError(f"Input changed after preparation: {item['path']}")
    return data

def array_id(value):
    if not re.fullmatch(r'\d+', value):
        raise RuntimeError('Expected a numeric Slurm array job ID')
    return value
