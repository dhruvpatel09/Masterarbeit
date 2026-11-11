#!/usr/bin/env python3
"""Geometry and fixed technical scope of a source-time validation pilot."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE_COMMIT='8b8b69863cd802fde5ea0d2b044aea5c98f91aa8'
CONFIGS=tuple(range(100,801,10))
PILOT=(170,640,650)
SOURCES=(58,90,122)
TAUS=tuple(range(33))
RADII=(4,12,24,32)
CUTOFFS=(40,60,80,100)
BASES={'pair':(40,100),'triple':(40,60,100),'full':CUTOFFS}
ANCHORS=(4,6,8)


def interval(t0,tau):
    if t0 not in SOURCES or tau not in TAUS:raise ValueError('Unknown source time or separation')
    t1=t0+tau
    if not 0<t0<=t1<191:raise ValueError('Interval crosses or touches an open temporal boundary')
    return dict(t0=t0,t1=t1,tau=tau,left_boundary_distance=t0,right_boundary_distance=191-t1)


def scans(cfg):
    if cfg not in PILOT:raise ValueError('Unknown pilot configuration')
    result=[dict(label='shift58_p10',t0=58,piter=10,purpose='source_time_shift'),
            dict(label='shift122_p10',t0=122,piter=10,purpose='source_time_shift')]
    if cfg==650:result += [dict(label='reference90_p10',t0=90,piter=10,purpose='fresh_reference_regression'),
                          dict(label='shift122_p40',t0=122,piter=40,purpose='projection_control')]
    return result


def audit_times(cfg):return sorted({s['t0']+tau for s in scans(cfg) for tau in TAUS})


def matrix_keys():
    return {(source,sink,t,r,a) for source in CUTOFFS for sink in CUTOFFS
            for t in TAUS for r in RADII for a in (0,1,2,-1)}


def build_plan():
    cases=[dict(configuration=cfg,new_scans=scans(cfg),saved_mode_audit_times=audit_times(cfg),
                reused_t0=90,reused_piter=10) for cfg in PILOT]
    return dict(base_commit=BASE_COMMIT,stage='design_and_input_preflight',
        native_extension_implemented=False,submission_implemented=False,
        lattice=[192,64,64,64],boundary_conditions=dict(time='open',space='periodic'),
        source_times=list(SOURCES),separations=list(TAUS),radii=list(RADII),
        main_target_radii=[24,32],short_distance_controls=[4,12],
        source_cutoffs=list(CUTOFFS),sink_cutoffs=list(CUTOFFS),
        matrix_orientation='sink row, source column; retain all 16 ordered pairs and complex values',
        pilot_configurations=list(PILOT),future_production_target=list(CONFIGS),
        training_configurations=[c for c in CONFIGS if c not in PILOT],
        fixed_source_bases={name:list(v) for name,v in BASES.items()},training_times=list(ANCHORS),
        primary_fixed_source_basis=['pair','full'],primary_training_time=4,
        fixed_source_policy='Train only at t0=90 on the 68 configurations outside the pilot. Freeze the same vectors for every source position. No pilot retuning.',
        reference_policy='Reuse validated t0=90 data on all three pilot configurations; require a fresh n650 t0=90 p10 regression before interpreting shifted results.',
        per_scan_matrix_records=len(matrix_keys()),planned_scan_launches=sum(len(c['new_scans']) for c in cases),
        planned_audit_launches=len(cases),planned_total_mpi_launches=sum(len(c['new_scans'])+1 for c in cases),
        planned_saved_mode_slices=sum(len(c['saved_mode_audit_times']) for c in cases),saved_modes_per_slice=100,
        cases=cases,proposed_allocation=dict(nodes_per_task=8,ranks_per_task=64,threads_per_rank=1,maximum_concurrent_tasks=3),
        statistical_unit='configuration, not spatial site, axis or source time',
        interpretation='Three-configuration technical pilot; no seven-point correlated fit, final potential, proof of boundary effects or automatic full71 submission.')


if __name__=='__main__':print(json.dumps(build_plan(),indent=2))
