#!/usr/bin/env python3
"""Standard-library replay of every admitted qX extension case."""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'inputs/job5958012'
CONFIGS=list(range(100,801,10))
CUTOFFS=[40,60,80,100]
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


class CSVHash:
    def __init__(self):self.hash=hashlib.sha256();self.writer=None;self.count=0
    def write(self,text):self.hash.update(text.encode('utf-8'));return len(text)
    def add(self,rows):
        if rows and self.writer is None:
            self.writer=csv.DictWriter(self,fieldnames=list(rows[0]),lineterminator='\n');self.writer.writeheader()
        self.writer.writerows(rows);self.count+=len(rows)
    def verify(self,path):
        if self.hash.hexdigest()!=digest(path):raise ValueError('Numerical export differs from raw replay: '+str(path))


def replay(progress=False):
    inventory=json.loads((HERE/'inputs/ARCHIVE_FILES.json').read_text());mapping={}
    for entry in inventory['files']:
        path=HERE/entry['stored']
        if path.stat().st_size!=entry['bytes'] or digest(path)!=entry['sha256']:
            raise ValueError('Archived input changed: '+entry['original'])
        mapping[entry['original']]=path
    if len(mapping)!=831:raise ValueError('Incomplete review archive')
    sources=json.loads((HERE/'inputs/SOURCES.json').read_text())
    if digest(SOURCE/'INPUTS.json')!=sources['measurement_inputs_sha256']:
        raise ValueError('Source receipt differs from published measurement')
    spec=importlib.util.spec_from_file_location('qx_extension_captured_campaign',str(SOURCE/'campaign.py'))
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(SOURCE))
    try:spec.loader.exec_module(module)
    finally:sys.path.pop(0)
    production=SOURCE/'campaigns/late_modes_full71';pilot=SOURCE/'campaigns/late_modes_pilot'
    plan=json.loads((production/'PREPARED.json').read_text());pilot_plan=json.loads((pilot/'PREPARED.json').read_text())
    build=json.loads((SOURCE/'build_provenance/MANIFEST.json').read_text())
    for target,data in ((production,plan),(pilot,pilot_plan)):
        if digest(target/'PREPARED.json')!=(target/'prepared.sha256').read_text().strip():raise ValueError('Prepared receipt changed')
        if digest(SOURCE/'build_provenance/MANIFEST.json')!=data['build_manifest_sha256']:raise ValueError('Build receipt changed')
    if build['inputs_sha256']!=digest(SOURCE/'INPUTS.json'):raise ValueError('Build source pins changed')
    if plan['stage']!='full' or plan['target_configurations']!=CONFIGS or len(plan['cases'])!=68:
        raise ValueError('Unexpected production inventory')
    if plan['pilot']['prepared_sha256']!=digest(pilot/'PREPARED.json') or plan['reused_configurations']!=[170,640,650]:
        raise ValueError('Reused pilot receipt changed')
    if [e['configuration'] for e in pilot_plan['cases']]!=[170,640,650]:raise ValueError('Wrong pilot inventory')
    g=plan['geometry']
    if g['taus']!=list(range(33)) or g['R']!=RADII or g['cutoffs']!=CUTOFFS or g['saved_Nv']!=100:
        raise ValueError('Unexpected geometry')
    checked={}
    for label,files in [('source',json.loads((SOURCE/'INPUTS.json').read_text())['files']),('build',build['files'])]:
        available=0
        for name,expected in files.items():
            if name in mapping:
                if digest(mapping[name])!=expected:raise ValueError('Captured source or build changed: '+name)
                available+=1
        checked[label]=dict(archived_files_verified=available,external_hashes_retained=len(files)-available)
    refs=module.references();data={};summaries=[];timing=[];groups=set();total_launches=0;validated_rows=0
    raw_hash=CSVHash();blocks_hash=CSVHash()
    for target,current in ((production,plan),(pilot,pilot_plan)):
        pilot_raw=CSVHash();pilot_blocks=CSVHash()
        for entry in current['cases']:
            cfg=entry['configuration']
            if cfg in groups:raise ValueError('Configuration measured twice in production cohort')
            groups.add(cfg);out=module.complete_case(target,entry)
            if out is None:raise ValueError('Missing validated case')
            if target==pilot:
                receipt=next(r for r in plan['pilot']['cases'] if r['configuration']==cfg)
                expected='/qx_static_extension/campaigns/late_modes_pilot/case_%03d/%s'%(entry['index'],out.name)
                if not receipt['attempt'].endswith(expected):raise ValueError('Reused attempt differs from production plan')
            rows,blocks,summary=module.replay_case(out,entry,refs)
            stored=json.loads((out/'SUMMARY.json').read_text())
            if any(stored.get(k)!=v for k,v in summary.items()):raise ValueError('Case summary differs from replay')
            for values,name in ((rows,'correlators.csv'),(blocks,'spatial_blocks.csv')):
                receipt=CSVHash();receipt.add(values);receipt.verify(out/name)
            if target==pilot:pilot_raw.add(rows);pilot_blocks.add(blocks)
            production_rows=[r for r in rows if r['piter']=='10'];production_blocks=[r for r in blocks if r['piter']=='10']
            raw_hash.add(production_rows);blocks_hash.add(production_blocks)
            for row in production_rows:
                if row['scheme']=='fixed' and row['axis']=='-1':
                    key=tuple(int(row[k]) for k in ('Nv','tau','R'));bucket=data.setdefault(key,{})
                    if cfg in bucket:raise ValueError('Duplicate combined-axis primary')
                    bucket[cfg]=row['L_RE']
            launches=json.loads((out/'LAUNCHES.json').read_text());total_launches+=len(launches);validated_rows+=len(rows)
            times={r['piter']:r['seconds'] for r in launches}
            timing.append(dict(configuration=cfg,audit_seconds=times[0],p10_seconds=times[10],
                p40_control_seconds=times.get(40),production_case_seconds=times[0]+times[10],
                origin_job=int(launches[0]['job_id'])))
            summaries.append(summary)
            if progress and len(summaries)%10==0:print('REPLAY %d/71 configurations'%len(summaries),flush=True)
        if target==pilot:
            pilot_raw.verify(pilot/'correlators.csv');pilot_blocks.verify(pilot/'spatial_blocks.csv')
    raw_hash.verify(production/'correlators.csv');blocks_hash.verify(production/'spatial_blocks.csv')
    summary=json.loads((production/'SUMMARY.json').read_text())
    if not summary['success'] or summary['cases']!=summaries or summary['missing_configurations'] or summary['validated_configurations']!=CONFIGS:
        raise ValueError('Collected summary differs from replay')
    if (raw_hash.count,blocks_hash.count,total_launches,validated_rows)!=(452412,340800,145,457272):
        raise ValueError('Wrong measurement counts')
    if groups!=set(CONFIGS) or len(data)!=1584 or any(set(v)!=set(CONFIGS) for v in data.values()):
        raise ValueError('Unmatched full-ensemble configuration samples')
    diagnostics=dict(success=True,configuration_count=71,configuration_labels=CONFIGS,
        raw_production_rows=raw_hash.count,regional_production_rows=blocks_hash.count,
        all_validated_rows_including_p40_controls=validated_rows,mpi_launches=total_launches,
        audited_modes=71*3300,archived_reference_rows=sum(s['archived_reference_rows'] for c in summaries for s in c['scans']),
        archived_reference_maximum_difference=max(s['archived_reference_maximum_absolute_difference'] for c in summaries for s in c['scans']),
        maximum_relative_mode_residual=max(c['quality']['maximum_relative_residual'] for c in summaries),
        maximum_unitarity_defect=max(float(s['unit_max']) for c in summaries for s in c['scans']),
        maximum_determinant_defect=max(float(s['det_max']) for c in summaries for s in c['scans']),
        maximum_p10_p40_component_difference=max(c['p10_p40_maximum_absolute_difference'] or 0 for c in summaries),
        direct_checks=sum(s['direct_tests'] for c in summaries for s in c['scans']),
        periodic_roundtrips=sum(s['periodic_roundtrips'] for c in summaries for s in c['scans']),provenance=checked)
    return data,diagnostics,sorted(timing,key=lambda r:r['configuration'])


if __name__=='__main__':
    _,diagnostics,_=replay(progress=True)
    print(json.dumps(diagnostics,indent=2))
