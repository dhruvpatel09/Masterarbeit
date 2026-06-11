#!/bin/bash
#SBATCH --job-name=tau_one_point
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=compute2011
#SBATCH --mem=8G
#SBATCH --chdir=/home/m2130292/Masterarbeit/laplace_static
#SBATCH --output=logs/tau_one_point_%j.out
#SBATCH --error=logs/tau_one_point_%j.err

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

mpirun -np 1 ./build/test_tau_one_point

echo "END: $(date)"
