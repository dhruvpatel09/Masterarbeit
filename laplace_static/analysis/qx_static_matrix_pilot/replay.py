"""Replay the immutable MPI pilot archive without accessing cluster data paths."""
import csv
import hashlib
import importlib.util
import json
import math
import sys
import tempfile
import zipfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
VALIDATION=ROOT/'laplace_static/validation/qx_static_matrix'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def load_campaign():
    previous=sys.path[:];sys.path.insert(0,str(VALIDATION))
    try:
        spec=importlib.util.spec_from_file_location('qx_pilot_matrix_campaign',str(VALIDATION/'campaign.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module
    finally:sys.path[:]=previous


def compare(actual,expected,path='root',rtol=5e-11,atol=1e-20):
    if isinstance(expected,bool) or expected is None or isinstance(expected,(str,int)):
        if actual!=expected:raise ValueError('Replay mismatch at '+path)
    elif isinstance(expected,float):
        if not isinstance(actual,(int,float)) or not math.isfinite(actual) or abs(actual-expected)>atol+rtol*max(abs(actual),abs(expected)):
            raise ValueError('Numerical replay mismatch at '+path)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(actual)!=len(expected):raise ValueError('Replay inventory mismatch at '+path)
        for i,(a,b) in enumerate(zip(actual,expected)):compare(a,b,path+'[%d]'%i,rtol,atol)
    elif isinstance(expected,dict):
        if not isinstance(actual,dict) or set(actual)!=set(expected):raise ValueError('Replay fields mismatch at '+path)
        for key in expected:compare(actual[key],expected[key],path+'.'+key,rtol,atol)
    else:raise TypeError('Unsupported replay value at '+path)


def extract_checked(archive,destination):
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ValueError('Archive CRC failed')
        names=set()
        for info in z.infolist():
            path=Path(info.filename)
            if path.is_absolute() or '..' in path.parts or info.filename in names or not info.filename.startswith('laplace_static/validation/qx_static_matrix/'):
                raise ValueError('Unexpected or unsafe archive path')
            if ((info.external_attr>>16)&0o170000)==0o120000:raise ValueError('Linked archive entry')
            names.add(info.filename)
        z.extractall(destination)


def replay():
    sources=json.loads((HERE/'inputs/SOURCES.json').read_text());archive=HERE/sources['review_archive']['path']
    if archive.stat().st_size!=sources['review_archive']['bytes'] or digest(archive)!=sources['review_archive']['sha256']:raise ValueError('Original review archive changed')
    c=load_campaign();c.verify_sources()
    reference=c.check.load_references(c.REFERENCE,c.PILOT)
    values=[];blocks=[];cases=[];timing=[]
    with tempfile.TemporaryDirectory() as temporary:
        dest=Path(temporary);extract_checked(archive,dest)
        source=dest/'laplace_static/validation/qx_static_matrix';target=source/'campaigns/pilot'
        pins=json.loads((source/'INPUTS.json').read_text())
        if (source/'INPUTS.json').read_bytes()!=(VALIDATION/'INPUTS.json').read_bytes():raise ValueError('Pilot source manifest differs from the published validation workflow')
        for name,expected in pins['files'].items():
            local=dest/name
            if local.is_file() and digest(local)!=expected:raise ValueError('Archived source changed: '+name)
        plan=json.loads((target/'PREPARED.json').read_text());archived=json.loads((target/'SUMMARY.json').read_text())
        if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared checksum mismatch')
        if plan['geometry']!=c.geometry() or [r['configuration'] for r in plan['cases']]!=list(c.PILOT):raise ValueError('Changed pilot geometry')
        build=json.loads((source/'build/MANIFEST.json').read_text())
        if digest(source/'build/MANIFEST.json')!=plan['build_manifest_sha256'] or build['inputs_sha256']!=digest(source/'INPUTS.json'):raise ValueError('Build/source binding changed')
        archived_build=[];unavailable_build=[]
        for name,expected in build['files'].items():
            path=dest/name
            if path.is_file():
                if digest(path)!=expected:raise ValueError('Archived build input changed: '+name)
                archived_build.append(name)
            else:unavailable_build.append(name)
        for entry in plan['cases']:
            out=c.complete_case(target,entry)
            if out is None:raise ValueError('Missing complete pilot case')
            rows,parts,summary=c.replay_case(out,entry,reference)
            for name,expected in [('matrix_correlators.csv',rows),('rank_blocks.csv',parts)]:
                with (out/name).open(newline='') as f:
                    if list(csv.DictReader(f))!=expected:raise ValueError('Case table does not replay exactly: '+name)
            recorded=json.loads((out/'SUMMARY.json').read_text())
            compare({key:recorded[key] for key in summary},summary,'case_summary')
            values+=rows;blocks+=parts;cases.append(summary)
            launches=json.loads((out/'LAUNCHES.json').read_text());seconds=sum(r['seconds'] for r in launches)
            timing.append(dict(configuration=entry['configuration'],audit_seconds=sum(r['seconds'] for r in launches if r['kind']=='audit'),p10_seconds=sum(r['seconds'] for r in launches if r['piter']==10),p40_seconds=sum(r['seconds'] for r in launches if r['piter']==40),total_mpi_seconds=seconds,launches=len(launches)))
        for name,expected in [('matrix_correlators.csv',values),('rank_blocks.csv',blocks)]:
            with (target/name).open(newline='') as f:
                if list(csv.DictReader(f))!=expected:raise ValueError('Combined table does not replay exactly: '+name)
        compare(archived['cases'],cases,'combined_summary.cases')
        if not archived['success'] or archived['missing_configurations'] or archived['matrix_rows']!=len(values) or archived['rank_block_rows']!=len(blocks):raise ValueError('Incomplete or inconsistent combined summary')
        scans=[scan for case in cases for scan in case['scans']]
        validation=dict(success=True,configurations=list(c.PILOT),scans=4,mpi_launches=7,audited_modes=9900,
            matrix_rows=len(values),rank_block_rows=len(blocks),p10_diagonal_reference_rows=sum(s['diagonal_reference_rows'] for s in scans if s['piter']==10),
            p10_diagonal_maximum_absolute_difference=max(s['diagonal_reference_maximum_absolute_difference'] for s in scans if s['piter']==10),
            p40_diagonal_maximum_absolute_difference=max(s['diagonal_reference_maximum_absolute_difference'] for s in scans if s['piter']==40),
            p10_p40_full_matrix_maximum_component_difference=cases[-1]['p10_p40_maximum_component_difference'],
            direct_tests=sum(s['direct_tests'] for s in scans),legacy_tests=sum(s['legacy_tests'] for s in scans),
            equal_time_gram_checks=sum(s['equal_time_gram_checks'] for s in scans),periodic_roundtrips=sum(s['periodic_roundtrips'] for s in scans),
            maximum_saved_mode_relative_residual=max(case['audit']['maximum_relative_residual'] for case in cases),
            maximum_central_link_unitarity=max(s['unit_max'] for s in scans),maximum_central_link_determinant=max(s['det_max'] for s in scans),
            timing=timing,archived_build_files_verified=archived_build,build_artifacts_recorded_but_not_in_review=unavailable_build,
            scope='Raw-log and provenance replay. Physical gauge/eigenvector paths and unavailable native binaries/libraries are not re-opened or re-executed offline.')
    return values,blocks,validation
