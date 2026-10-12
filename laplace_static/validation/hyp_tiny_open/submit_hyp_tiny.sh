#!/usr/bin/env bash
#SBATCH --job-name=hyp_tiny
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --qos=micro
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_tiny_open/hyp_tiny_%j.out
#SBATCH --error=laplace_static/validation/hyp_tiny_open/hyp_tiny_%j.err

set -euo pipefail
source env/supermuc.sh
module load slurm_setup

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

mpiexec -n "$SLURM_NTASKS" \
    solver/build_supermuc_tiny/test_hyp_tiny \
    laplace_static/validation/tiny_open8x4x4x4/data/cnfg/open8x4x4x4n1
