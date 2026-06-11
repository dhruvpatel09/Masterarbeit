#!/bin/bash
#SBATCH --job-name=laplace_pair_R1T1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=compute2011
#SBATCH --mem=8G
#SBATCH --chdir=/home/m2130292/Masterarbeit/laplace_static
#SBATCH --output=logs/laplace_pair_R1T1_%j.out
#SBATCH --error=logs/laplace_pair_R1T1_%j.err

set -euo pipefail

module purge
module load compiler/gcc/10.2.1
module load mpi/openmpi/4.1.0-no_ucx
module load mkl/latest

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_THREADING_LAYER=sequential

echo "START: $(date)"
echo "HOSTNAME=$(hostname)"
echo "PWD=$(pwd)"

export NCFG="${NCFG:-${SLURM_ARRAY_TASK_ID:-1}}"
echo "NCFG=${NCFG}"
mpirun -x NCFG -np 1 ./build/test_laplace_pair_R1T1

echo "END: $(date)"
