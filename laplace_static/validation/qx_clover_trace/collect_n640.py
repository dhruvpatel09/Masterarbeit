#!/usr/bin/env python3
"""Find and export n640 solver logs from repository and production storage."""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRODUCTION_RUNS = Path('/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs')
MATCH = re.compile(r'5924010|5925063|n640|(?:^|[/_.-])640(?:[/_.-]|$)')
CONTROL_MATCH = re.compile(r'5924010|5925063|/n640|qX_p2n640|"configuration"\s*:\s*640\b|\bfirst\s+640\b')
TEXT_SUFFIXES = {'.in', '.out', '.err', '.log', '.txt', '.json', '.sha256', '.sh', '.csv', '.dat'}
CONTROL_SUFFIXES = {'.in', '.json', '.sh'}
MAX_CONTROL = 1024*1024
MAX_FILE = 16*1024*1024
MAX_TOTAL = 64*1024*1024
PRUNE = {'eigenvectors', 'build', '.git', '__pycache__'}


def search_bases(root, extra_roots=()):
    # Controller folders need not contain n640 in their names.
    # Production runs are stored outside the repository.
    selected = [Path(root)/'mental', Path(root)/'laplace_static/validation']
    selected += [Path(p).expanduser().absolute() for p in extra_roots]
    unique = []
    for base in selected:
        if base not in unique:
            unique.append(base)
    return unique


def candidates(root, extra_roots=()):
    found = set()
    for base in search_bases(root, extra_roots):
        if not base.is_dir():
            continue
        for parent, dirs, files in os.walk(str(base), followlinks=False):
            dirs[:] = sorted(d for d in dirs if d not in PRUNE and
                             not (Path(parent)/d).is_symlink())
            for name in sorted(files):
                path = Path(parent)/name
                if path.is_symlink() or path.suffix.lower() not in TEXT_SUFFIXES:
                    continue
                if MATCH.search(str(path)):
                    found.add(path)
                elif path.suffix.lower() in CONTROL_SUFFIXES:
                    try:
                        if path.stat().st_size <= MAX_CONTROL:
                            text = path.read_text(encoding='utf-8')
                            if CONTROL_MATCH.search(text):
                                found.add(path)
                    except (OSError, UnicodeError):
                        continue
    return sorted(found)


def archive_name(path, root, bases):
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        for index, base in enumerate(bases):
            try:
                return 'external/search_root_%02d/%s' % (index, path.relative_to(base).as_posix())
            except ValueError:
                continue
    raise ValueError('File outside the recorded search roots: '+str(path))


def search_report(root, extra_roots=()):
    report = []
    for base in search_bases(root, extra_roots):
        row = dict(path=str(base), exists=base.is_dir())
        if row['exists']:
            try:
                entries = sorted(p.name for p in base.iterdir())
                row['top_level_entries'] = entries[:60]
                row['top_level_entry_count'] = len(entries)
            except OSError as exc:
                row['error'] = str(exc)
        report.append(row)
    return report


def export(root, output, extra_roots=()):
    root, output = Path(root), Path(output)
    bases = search_bases(root, extra_roots)
    paths = candidates(root, extra_roots)
    if not paths:
        raise ValueError('No matching n640 logs/inputs found; use --list to inspect searched locations')
    manifest = dict(scope='Existing n640 failure evidence; no job execution or input changes',
                    attempts=[5924010, 5925063], search_roots=search_report(root, extra_roots),
                    files=[], skipped=[])
    try:
        manifest['git_head'] = subprocess.check_output(
            ['git', '-C', str(root), 'rev-parse', 'HEAD'], universal_newlines=True,
            stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        manifest['git_head'] = None
    total, contents = 0, []
    provenance = root/'laplace_static/analysis/qx_eigenpair_completion/inputs/unresolved_n640.json'
    evidence_paths = set(paths)
    if provenance.is_file() and provenance not in evidence_paths:
        paths.append(provenance)
    for path in paths:
        name = archive_name(path, root, bases)
        try:
            size = path.stat().st_size
            if size > MAX_FILE or total+size > MAX_TOTAL:
                manifest['skipped'].append(dict(path=str(path), bytes=size, reason='size limit'))
                continue
            data = path.read_bytes()
            data.decode('utf-8')
            if b'\x00' in data:
                raise ValueError('NUL byte in text candidate')
        except (OSError, UnicodeError, ValueError) as exc:
            manifest['skipped'].append(dict(path=str(path), reason=str(exc)))
            continue
        manifest['files'].append(dict(path=str(path), archive_path=name, bytes=len(data),
                                      sha256=hashlib.sha256(data).hexdigest()))
        total += len(data)
        contents.append((path, name, data))
    if not any(path in evidence_paths for path, _, _ in contents):
        raise ValueError('No readable n640 evidence; no archive written')
    with zipfile.ZipFile(str(output), 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr('N640_EVIDENCE.json', json.dumps(manifest, indent=2)+'\n')
        for _, name, data in contents:
            archive.writestr(name, data)
    print('REVIEW_ARCHIVE', output.resolve())
    print('Included %d files; skipped %d (listed in N640_EVIDENCE.json)' %
          (len(contents), len(manifest['skipped'])))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='Show search locations and candidate paths')
    parser.add_argument('--search-root', type=Path, action='append', default=[],
                        help='Additional production/controller directory; repeatable')
    parser.add_argument('--output', type=Path, default=ROOT/'qx_n640_diagnostics_review.zip')
    args = parser.parse_args()
    extra = [PRODUCTION_RUNS]+args.search_root
    if args.list:
        print(json.dumps(dict(search_roots=search_report(ROOT, extra)), indent=2))
        paths = candidates(ROOT, extra)
        for path in paths:
            print(path)
        print('Matched text files:', len(paths))
    else:
        try:
            export(ROOT, args.output, extra)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            print(json.dumps(dict(search_roots=search_report(ROOT, extra)), indent=2))
            sys.exit(2)


if __name__ == '__main__':
    main()
