"""Explicit source/sink and validation geometry of the technical pilot."""
CUTOFFS=(40,60,80,100)
RADII=(4,12,24,32)
TAUS=tuple(range(33))
AUDIT_TIMES=tuple(range(90,123))
DIRECT_TAUS=(0,4,16,32)
BLOCK_TAUS=(0,32)
BLOCK_RADII=(4,32)
PILOT=(170,640,650)
AXES=(0,1,2,-1)
SITES=(0,4095)


def expected_keys():
    return {(s,k,t,r,a) for s in CUTOFFS for k in CUTOFFS for t in TAUS for r in RADII for a in AXES}


def direct_keys():
    return {(s,k,t,r,a,x) for s in CUTOFFS for k in CUTOFFS for t in DIRECT_TAUS for r in RADII for a in (0,1,2) for x in SITES}


def legacy_keys():return {key for key in direct_keys() if key[0]==key[1]}


def block_keys():
    return {(s,k,t,r,a) for s in CUTOFFS for k in CUTOFFS for t in BLOCK_TAUS for r in BLOCK_RADII for a in (0,1,2)}
