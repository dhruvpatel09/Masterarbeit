#!/usr/bin/env python3
"""Recompute derived tables from the native-replayed primary matrix cache."""
import platform
import sys
import numpy as np
import scipy
from replay import HERE
from analysis import analyze
from statistics_tools import write_csv, write_json

if __name__ == '__main__':
    data=np.load(str(HERE/'inputs/primary_matrices.npz'),allow_pickle=False)
    tables,arrays,summary=analyze(data['C'],data['diagonal_text'],True)
    for name,values in tables.items():
        write_csv(HERE/'results'/(name+'.csv'),values)
    np.savez_compressed(str(HERE/'results/projection_resamples.npz'),**arrays)
    write_json(HERE/'results/SUMMARY.json',summary)
    write_json(HERE/'inputs/ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()))
    print('PASS derived full71 analysis:',len(tables),'tables;',summary['projected_effective_energy_rows'],'projected energy rows.')
