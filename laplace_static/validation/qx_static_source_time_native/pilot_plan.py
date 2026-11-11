"""Native execution order and additional rank-block coverage of the fixed pilot."""
from parents import design
PILOT=design.PILOT
CUTOFFS=design.CUTOFFS
RADII=design.RADII
TAUS=design.TAUS
SOURCES=design.SOURCES
AXES=(0,1,2,-1)
DIRECT_TAUS=(0,4,16,32)
BLOCK_TAUS=(0,32)
BLOCK_RADII=(4,24,32)


def execution_scans(cfg):
    return sorted(design.scans(cfg),key=lambda s:s['purpose']!='fresh_reference_regression')


def expected_keys():return design.matrix_keys()
def direct_keys():return {(s,k,t,r,a,x) for s in CUTOFFS for k in CUTOFFS for t in DIRECT_TAUS for r in RADII for a in (0,1,2) for x in (0,4095)}
def legacy_keys():return {key for key in direct_keys() if key[0]==key[1]}
def block_keys():return {(s,k,t,r,a) for s in CUTOFFS for k in CUTOFFS for t in BLOCK_TAUS for r in BLOCK_RADII for a in (0,1,2)}


def geometry():
    return dict(protocol='source_time_v1',source_times=list(SOURCES),taus=list(TAUS),radii=list(RADII),
        source_cutoffs=list(CUTOFFS),sink_cutoffs=list(CUTOFFS),block_taus=list(BLOCK_TAUS),
        block_radii=list(BLOCK_RADII),axis_choices=list(AXES),saved_modes=100,
        nodes=8,ranks=64,threads=1,pilot_configurations=list(PILOT),
        source_plan=design.build_plan(),new_matrix_records=8*len(expected_keys()),
        new_rank_block_records=8*len(block_keys())*64)
