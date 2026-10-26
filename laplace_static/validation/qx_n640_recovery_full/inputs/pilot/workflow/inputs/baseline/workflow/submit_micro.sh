#!/usr/bin/env bash
#SBATCH --job-name=qx_n640_diagnostic
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --qos=micro
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/qx_n640_diagnostic/diagnostic_%j.out
#SBATCH --error=laplace_static/validation/qx_n640_diagnostic/diagnostic_%j.err
#SBATCH --no-requeue
set -euo pipefail
test "${SLURM_JOB_NUM_NODES:?}" -eq 8
test "${SLURM_NTASKS:?}" -eq 64
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
export I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core
python3 laplace_static/validation/qx_n640_diagnostic/campaign.py run --name "${1:?campaign name required}"
