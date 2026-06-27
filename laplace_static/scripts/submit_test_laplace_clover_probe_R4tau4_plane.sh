#!/bin/bash
#SBATCH --job-name=lap_clov_R4t4
#SBATCH --partition=compute2011
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

cd /home/m2130292/Masterarbeit/laplace_static

module purge
module load compiler/gcc/10.2.1
module load mpi/openmpi/4.1.0-no_ucx
module load mkl/latest

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_THREADING_LAYER=sequential

export NCFG="${NCFG:-${SLURM_ARRAY_TASK_ID:-1}}"
export RSEP="${RSEP:-4}"
export TAU="${TAU:-4}"
export T0_MODE="${T0_MODE:-average}"
export T0_START="${T0_START:-0}"
export T0_STRIDE="${T0_STRIDE:-1}"
export DX_MIN="${DX_MIN:--2}"
export DX_MAX="${DX_MAX:-6}"
export DY_MIN="${DY_MIN:--6}"
export DY_MAX="${DY_MAX:-6}"
export DZ_FIXED="${DZ_FIXED:-0}"

echo "HOST=$(hostname)"
echo "DATE=$(date)"
echo "JOB_ID=${SLURM_JOB_ID:-none}"
echo "ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-none}"
echo "NCFG=${NCFG}"
echo "RSEP=${RSEP}"
echo "TAU=${TAU}"
echo "T0_MODE=${T0_MODE}"
echo "T0_START=${T0_START}"
echo "T0_STRIDE=${T0_STRIDE}"
echo "DX_MIN=${DX_MIN}"
echo "DX_MAX=${DX_MAX}"
echo "DY_MIN=${DY_MIN}"
echo "DY_MAX=${DY_MAX}"
echo "DZ_FIXED=${DZ_FIXED}"

mpirun \
  -x NCFG \
  -x RSEP \
  -x TAU \
  -x T0_MODE \
  -x T0_START \
  -x T0_STRIDE \
  -x DX_MIN \
  -x DX_MAX \
  -x DY_MIN \
  -x DY_MAX \
  -x DZ_FIXED \
  -np 1 \
  ./build/test_laplace_clover_probe_R4tau4_plane
