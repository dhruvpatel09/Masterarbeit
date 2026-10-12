"""Revalidate archived campaigns without changing the supplied evidence."""
import contextlib
import csv
import hashlib
import io
import json
import shutil
import sys
import tempfile
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root, measurement_name, queue_name, first_job, out = sys.argv[1:]
    root, out = Path(root).resolve(), Path(out)
    case = root / 'laplace_static/validation' / measurement_name
    queue = case.with_name(queue_name)
    sys.path.insert(0, str(case))
    import common as c
    import check_task as check
    ledger = json.loads((queue / 'RUNS.json').read_text())
    collection = queue / ('collection_mixed_' + first_job)
    saved_collection = json.loads((collection / 'SUMMARY.json').read_text())
    assert ledger['first_job'] == first_job
    assert saved_collection['success'] is True
    entries = ledger['cases']
    ncases = len(c.CONFIGS) * len(c.TAUS)
    assert len(entries) == ncases
    assert len({e['array_job_id'] for e in entries}) == ncases
    remote_prefix = '/hppfs/work/pn29se/di24gus/Masterarbeit/'
    present, absent = [], []
    for manifest in (case/'package.sha256', case/'prepared.sha256',
                     case/'build/manifest.sha256', queue/'package.sha256', queue/'inputs.sha256'):
        for line in manifest.read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            relative = name[len(remote_prefix):] if name.startswith(remote_prefix) else name
            target = root / relative
            if target.is_file():
                assert sha(target) == digest, ('Checksum mismatch', manifest, name)
                present.append({'manifest': str(manifest.relative_to(root)), 'file': name})
            else:
                absent.append({'manifest': str(manifest.relative_to(root)), 'file': name})
    for name, digest in ledger['provenance'].items():
        assert sha(case/name) == digest
    results, launches_all = [], []
    collected = {name: [] for name in ('profiles_raw.csv', 'profiles_paired.csv',
                                      'su3_times.csv', 'su3_windows.csv')}
    for i, entry in enumerate(entries):
        selected = c.case_for(i)
        assert all(entry[k] == value for k, value in selected.items())
        assert entry['state'] == 'completed'
        job = entry['array_job_id']
        folder = case/f'results_{job}_{i}'
        assert not (folder/'FAILED.json').exists()
        assert json.loads((folder/'CASE.json').read_text()) == dict(selected, array_job_id=job)
        for name in ledger['provenance']:
            assert (folder/name.replace('/', '_')).read_bytes() == (case/name).read_bytes()
        launches = json.loads((folder/'LAUNCHES.json').read_text())
        wanted = [(selected['configuration'], selected['tau'], p, f'p{p}') for p in (4,40)]
        if i == 0:
            wanted += [(100,4,p,f'n100_p{p}') for p in (4,40)]
        assert [(r['configuration'],r['tau'],r['piter'],r['label']) for r in launches] == wanted
        assert all(r['exit_code'] == 0 and r['seconds'] > 0 for r in launches)
        launches_all.extend(launches)
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            for launch in launches:
                shutil.copyfile(folder/(launch['label']+'.out'), work/(launch['label']+'.out'))
            with contextlib.redirect_stdout(io.StringIO()):
                computed = check.check(work, selected['configuration'], selected['tau'], regression=(i==0))
            saved = json.loads((folder/'SUMMARY.json').read_text())
            assert computed == saved, ('Summary disagreement', i)
            for name in collected:
                assert (work/name).read_bytes() == (folder/name).read_bytes(), ('Export disagreement', i, name)
                with (work/name).open(newline='') as stream:
                    collected[name].extend(csv.DictReader(stream))
        computed.update(index=i, array_job_id=job, mpi_seconds=sum(r['seconds'] for r in launches))
        assert computed == saved_collection['cases'][i]
        results.append(computed)
    counts = {}
    for name, rows in collected.items():
        with (collection/name).open(newline='') as stream:
            assert rows == list(csv.DictReader(stream)), ('Collection differs', name)
        counts[name] = len(rows)
    assert counts == {'profiles_raw.csv': ncases*664, 'profiles_paired.csv': ncases*332,
                      'su3_times.csv': len(c.CONFIGS)*2*sum(c.TAUS), 'su3_windows.csv': ncases*2}
    result = dict(success=True, configurations=list(c.CONFIGS), cases=ncases,
                  local_manifest_entries_verified=len(present),
                  manifest_entries_unavailable_in_review=absent, cases_rechecked=results,
                  launch_count=len(launches_all), counts=counts,
                  archive_scope='Archived source hashes, matching manifests/ledger, actual stdout checks and byte-exact exports. Runtime binaries, gauge fields and eigenvector HDF5 bytes are not in these review archives.')
    out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases_rechecked','manifest_entries_unavailable_in_review')}))


if __name__ == '__main__':
    main()
