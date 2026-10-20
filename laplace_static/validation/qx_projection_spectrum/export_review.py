#!/usr/bin/env python3
"""Validate a completed campaign and atomically export a complete review ZIP."""
import argparse
import json
import os
import shutil
import tempfile
import zipfile
from pathlib import Path
import campaign as c

SMALL_LIMIT = 8*1024**2
COLLECTION_LIMIT = 64*1024**2


def check_files(files, large_exports):
    for path in files:
        if path.is_symlink() or path.suffix in ('.h5', '.a', '.o', '.zip') or '__pycache__' in path.parts:
            raise ValueError('Unexpected review file: '+str(path))
        limit = COLLECTION_LIMIT if path in large_exports else SMALL_LIMIT
        if path.stat().st_size > limit:
            raise ValueError('Review file exceeds its validated size limit: '+str(path))


def verify_existing_archive(archive, files, root):
    if not archive.exists():
        return False
    expected = {p.relative_to(root).as_posix(): p for p in files}
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        if len(names) != len(set(names)) or not set(names).issubset(expected):
            raise ValueError('Existing review archive contains unrelated entries')
        for name in names:
            if bundle.read(name) != expected[name].read_bytes():
                raise ValueError('Existing review archive differs from validated data: '+name)
    return set(names) == set(expected)


def export_archive(archive, files, root, large_exports):
    check_files(files, large_exports)
    if verify_existing_archive(archive, files, root):
        return 'already complete'
    descriptor, name = tempfile.mkstemp(prefix='.'+archive.stem+'_', suffix='.tmp', dir=str(archive.parent))
    os.close(descriptor)
    temporary = Path(name)
    try:
        with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(files):
                bundle.write(path, path.relative_to(root).as_posix())
        with zipfile.ZipFile(temporary) as bundle:
            if bundle.testzip() is not None or len(bundle.namelist()) != len(files):
                raise ValueError('New review archive failed verification')
        os.replace(str(temporary), str(archive))
    finally:
        if temporary.exists():
            temporary.unlink()
    return 'complete'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True)
    args = parser.parse_args()
    folder, data = c.prepared(args.name)
    raw, windows, summaries = [], [], []
    prepared_hash = c.digest(folder/'PREPARED.json')
    for entry in data['cases']:
        out = folder/('case_%03d' % entry['index'])
        if (out/'FAILED.json').exists():
            raise ValueError('Case has a failure record')
        identity = json.loads((out/'CASE.json').read_text())
        if any(identity[k] != v for k, v in entry.items()) or identity['prepared_sha256'] != prepared_hash:
            raise ValueError('Wrong case identity')
        launches = json.loads((out/'LAUNCHES.json').read_text())
        wanted = [(entry['configuration'], entry['tau'], p, 'p'+str(p)) for p in entry['piters']]
        if data['stage'] == 'pilot':
            wanted += [(100, 4, p, 'n100_p'+str(p)) for p in (4, 40)]
        if [(r['configuration'], r['tau'], r['piter'], r['label']) for r in launches] != wanted or any(r['exit_code'] for r in launches):
            raise ValueError('MPI launch inventory failed')
        result = c.read_case(out, entry, data)
        if json.loads((out/'SUMMARY.json').read_text()) != result:
            raise ValueError('Case summary differs from revalidated measurements')
        summaries.append(result)
        raw.extend(c.csv_rows(out/'profiles_raw.csv'))
        windows.extend(c.csv_rows(out/'su3_windows.csv'))
    check = c.checker()
    collection = folder/'collection'
    with tempfile.TemporaryDirectory(prefix='qx_review_collection_') as name:
        temporary = Path(name)
        check.write_csv(temporary/'profiles_raw.csv', check.RAW_COLUMNS, raw)
        check.write_csv(temporary/'su3_windows.csv', check.SU3_COLUMNS, windows)
        c.write_json(temporary/'SUMMARY.json', dict(
            success=True, stage=data['stage'], cases=summaries,
            raw_records=len(raw), window_records=len(windows)))
        collection.mkdir(exist_ok=True)
        for source in temporary.iterdir():
            target = collection/source.name
            if target.exists():
                if target.read_bytes() != source.read_bytes():
                    raise ValueError('Existing collection differs from validated exports: '+source.name)
            else:
                shutil.copyfile(str(source), str(target))
    files = [p for p in folder.rglob('*') if p.is_file()]
    files += [p for p in c.CASE.iterdir() if p.is_file() and p.suffix not in ('.out', '.err')]
    files = sorted(set(files))
    archive = c.ROOT/('qx_projection_spectrum_'+args.name+'_review.zip')
    state = export_archive(
        archive, files, c.ROOT,
        {collection/'profiles_raw.csv', collection/'su3_windows.csv'})
    print('PASS review export:', state, len(files), 'files')
    print('REVIEW_ARCHIVE', archive)


if __name__ == '__main__':
    main()
