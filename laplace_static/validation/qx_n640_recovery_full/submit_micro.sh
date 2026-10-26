#!/usr/bin/env bash
#SBATCH --job-name=qx_n640_full192
#SBATCH --partition=micro
#SBATCH --qos=micro
#SBATCH --account=pn29se
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --export=NONE
#SBATCH --no-requeue
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/qx_n640_recovery_full/slots_%A_%a.out
#SBATCH --error=laplace_static/validation/qx_n640_recovery_full/slots_%A_%a.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MKL_THREADING_LAYER=SEQUENTIAL I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core
unset MENTAL_EVAL_ONLY
python3 laplace_static/validation/qx_n640_recovery_full/campaign.py run-slot --name "$1" --slot "$SLURM_ARRAY_TASK_ID"
