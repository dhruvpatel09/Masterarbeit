#!/usr/bin/env bash
#SBATCH --job-name=qx_n100_mode_audit
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
#SBATCH --output=laplace_static/validation/qx_saved_modes_n100/audit_%j.out
#SBATCH --error=laplace_static/validation/qx_saved_modes_n100/audit_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
python3 laplace_static/validation/qx_saved_modes_n100/audit.py run
