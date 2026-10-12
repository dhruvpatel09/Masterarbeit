#!/usr/bin/env bash
#SBATCH --job-name=qx_hyp200230
#SBATCH --account=pn29se
#SBATCH --partition=test
#SBATCH --qos=test
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --array=0
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_projection_configs_qx_200_230/profiles_%A_%a.out
#SBATCH --error=laplace_static/validation/hyp_projection_configs_qx_200_230/profiles_%A_%a.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
python3 laplace_static/validation/hyp_projection_configs_qx_200_230/run_task.py
