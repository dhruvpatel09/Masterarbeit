#!/bin/bash
#SBATCH --job-name=clover_EB2
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

echo "HOST=$(hostname)"
echo "DATE=$(date)"
echo "JOB_ID=${SLURM_JOB_ID:-none}"
echo "ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-none}"
echo "NCFG=${NCFG}"

mpirun \
  -x NCFG \
  -np 1 \
  ./build/test_clover_EB2_one_config
