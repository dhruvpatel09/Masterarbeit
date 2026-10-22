#!/usr/bin/env bash
#SBATCH --job-name=qx_clover_trace
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/qx_clover_trace/pilot_%j.out
#SBATCH --error=laplace_static/validation/qx_clover_trace/pilot_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
python3 laplace_static/validation/qx_clover_trace/campaign.py run --name "$1"
