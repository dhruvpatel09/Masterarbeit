#!/usr/bin/env python3
"""Generate the complete offline stability package from pinned measurements."""
import sys,platform
import numpy as np
import scipy
from common import HERE,read,write_csv,write_json
from source_splits import make_assignments
from analysis import compute

if __name__=='__main__':
    plan=read(HERE/'ANALYSIS_PLAN.json');target=HERE/'inputs/resampling_plan.npz'
    if not target.exists():
        rng=np.random.RandomState(plan['bootstrap_seed'])
        np.savez_compressed(str(target),configurations=plan['configuration_labels'],
            bootstrap_indices=rng.randint(0,71,size=(plan['bootstrap_draws'],71)).astype(np.uint8),
            split_assignments=make_assignments(plan))
    tables,arrays,summary=compute(True)
    for name,values in tables.items():write_csv(HERE/'results'/(name+'.csv'),values)
    np.savez_compressed(str(HERE/'results/stability_resamples.npz'),**arrays)
    write_json(HERE/'results/SUMMARY.json',summary)
    write_json(HERE/'inputs/ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()))
    print('PASS offline matrix stability:',summary['total_curves'],'curves;',summary['correlator_fit_attempts'],'correlated fit attempts.')
