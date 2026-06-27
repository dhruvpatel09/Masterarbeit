#!/bin/bash
#SBATCH --job-name=laplace_RT
#SBATCH --partition=compute2011
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

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
export RMAX="${RMAX:-3}"
export TMAX="${TMAX:-4}"
export T0_MODE="${T0_MODE:-fixed}"
export T0_FIXED="${T0_FIXED:-0}"
export T0_START="${T0_START:-0}"
export T0_STRIDE="${T0_STRIDE:-1}"

echo "HOST=$(hostname)"
echo "DATE=$(date)"
echo "JOB_ID=${SLURM_JOB_ID:-none}"
echo "ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-none}"
echo "NCFG=${NCFG}"
echo "RMAX=${RMAX}"
echo "TMAX=${TMAX}"
echo "T0_MODE=${T0_MODE}"
echo "T0_FIXED=${T0_FIXED}"
echo "T0_START=${T0_START}"
echo "T0_STRIDE=${T0_STRIDE}"

mpirun \
  -x NCFG \
  -x RMAX \
  -x TMAX \
  -x T0_MODE \
  -x T0_FIXED \
  -x T0_START \
  -x T0_STRIDE \
  -np 1 \
  ./build/test_laplace_RT_scan
