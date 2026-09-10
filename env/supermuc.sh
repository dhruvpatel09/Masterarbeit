#!/usr/bin/env bash
#
# Canonical SuperMUC-NG software environment for Masterarbeit.
#
# Usage:
#   source env/supermuc.sh
#

module purge || true

module load stack/24.6.0
module load gcc/14.3.0
module load intel-mpi/2021.17.0
module load hdf5/1.14.5-gcc14-impi
module load intel-mkl/2025.3.0

# The project uses MPI parallelism. Keep oneMKL sequential within each
# MPI rank to avoid introducing an additional threaded parallel layer.
export MKL_THREADING_LAYER=SEQUENTIAL
