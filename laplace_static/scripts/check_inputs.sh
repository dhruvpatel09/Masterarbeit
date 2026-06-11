#!/bin/bash
set -euo pipefail

H5ROOT=/home/m2130292/Masterarbeit/extern/hdf5-1.14.6-openmpi410
BASE=/home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full

echo "Checking representative eigenvector/eigenvalue files..."
for CFG in 1 50 100; do
  echo
  echo "=== n${CFG} ==="
  echo "evec count:"
  ls "${BASE}/n${CFG}/eigenvectors"/*.h5 | wc -l

  echo "evec layout:"
  "${H5ROOT}/bin/h5ls" -r "${BASE}/n${CFG}/eigenvectors/Em1p4n${CFG}_evec_t0.h5" | head -25

  echo "eval layout:"
  "${H5ROOT}/bin/h5ls" -r "${BASE}/n${CFG}/eigenvalues"/*.h5 | head -25
done
