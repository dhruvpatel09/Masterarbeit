#!/usr/bin/env python3
"""Check existing audited runs and required inputs; create no eigenvectors."""
import json
from pathlib import Path
import common as c

def main():
    if (c.CASE / 'PREPARED.json').exists() or (c.CASE / 'prepared.sha256').exists():
        raise RuntimeError('Preparation already exists; no files overwritten')
    c.verify_manifest(c.CASE / 'package.sha256')
    c.verify_manifest(c.CASE / 'inputs.sha256')
    c.verify_manifest(c.CASE / 'build/manifest.sha256')
    audit = json.loads(c.AUDIT.read_text())
    if (audit.get('success') is not True or audit.get('mode_records') != 5200 or
            audit.get('slice_records') != 52 or
            sorted(r['configuration'] for r in audit['runs']) != list(c.CONFIGS) or
            any(r.get('success') is not True or r['mpi_exit_code'] != 0
                for r in audit['runs'])):
        raise RuntimeError('The saved-mode audit is not the expected successful audit')
    selected = json.loads((c.CASE / 'audited_runs.json').read_text())['runs']
    inputs = []
    for run in selected:
        cfg = run['configuration']
        marker = Path(run['run']) / 'COMPLETED.json'
        if c.digest(marker) != run['completion_marker_sha256']:
            raise RuntimeError(f'Audited completion marker changed: {marker}')
        done = json.loads(marker.read_text())
        if (done['configuration'], done['timeslices'], done['eigenvector_files']) != (cfg, 192, 192):
            raise RuntimeError(f'Unexpected inventory: {marker}')
        inputs.append(c.check_file(c.GAUGE_ROOT / f'qX_p2n{cfg}'))
        for t in (95, 96, 97, 99, 100, 101):
            inputs.append(c.check_file(Path(run['run']) / 'eigenvectors' /
                                       f'qX_p2n{cfg}_evec_t{t}.h5'))
    inputs.append(c.check_file(c.GAUGE_ROOT / 'qX_p2n100'))
    for t in (96, 100):
        inputs.append(c.check_file(c.N100_RUN / 'eigenvectors' / f'qX_p2n100_evec_t{t}.h5'))
    data = dict(runs=selected, cases=[c.case_for(i) for i in range(c.NCASES)],
                gauge_root=str(c.GAUGE_ROOT), n100_run=str(c.N100_RUN), inputs=inputs,
                audit_job=5920205, audit_scope=audit['scope'],
                settings=dict(Nv=100, midpoint=98, R=[1,2,4], piter=[4,40],
                              hyp_steps=1, alpha=[1,1,0.5], grid=[1,4,4,4], ranks=64),
                quality_scope='All tau2/4/6 endpoints and midpoint98 passed the saved-pair residual/Gram audit on 13 selected slices in job5920205. Other eigenvector times are outside that sampled numerical audit.')
    c.write_json(c.CASE / 'PREPARED.json', data)
    target = c.CASE / 'PREPARED.json'
    with (c.CASE / 'prepared.sha256').open('x') as stream:
        stream.write(c.digest(target) + '  ' + str(target.relative_to(c.ROOT)) + '\n')
    print('Prepared 12 paired cases for n160..n190, tau=2,4,6; task0 also checks n100 tau4.')
    print('All gauge/eigenvector inputs remain in their existing locations.')

if __name__ == '__main__':
    main()
