"""Geometry and cutoff inventory shared by generation and output validation."""
RADII=[1,2,4,8,12,16,20,23,24,25,28,32]
CUTOFFS=[40,60,80,100]
TAUS=list(range(33))
AUDIT_TIMES=list(range(90,123))
DIRECT_TAUS=[0,4,16,32]
DIRECT_RADII=[1,4,16,32]
BLOCK_TAUS=[0,8,16,24,32]
BLOCK_RADII=[4,8,16,24,32]


def cutoffs(piter):
    if piter not in (10,40):raise ValueError('Expected piter 10 or 40')
    return CUTOFFS if piter==10 else [100]


def expected_keys(piter):
    result={('fixed',n,t,r,a) for n in cutoffs(piter) for t in TAUS
            for r in RADII for a in (-1,0,1,2)}
    result|={('regression',100,t,r,a) for t in (2,4,6)
             for r in (1,2,4) for a in (-1,0,1,2)}
    return result


def direct_keys(piter):
    result={('fixed',n,t,r,a,s) for n in cutoffs(piter) for t in DIRECT_TAUS
            for r in DIRECT_RADII for a in range(3) for s in (0,4095)}
    result|={('regression',100,4,r,a,s) for r in (1,4)
             for a in range(3) for s in (0,4095)}
    return result
