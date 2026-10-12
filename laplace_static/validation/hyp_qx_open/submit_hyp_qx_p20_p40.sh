#!/usr/bin/env bash
#SBATCH --job-name=hyp_qx_2040
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --qos=micro
#SBATCH --time=00:15:00
#SBATCH --nodes=8
#SBATCH --ntasks=384
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_qx_open/hyp_qx_p20_p40_%j.out
#SBATCH --error=laplace_static/validation/hyp_qx_open/hyp_qx_p20_p40_%j.err

set -euo pipefail
source env/supermuc.sh
module load slurm_setup

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

gauge_file=/hppfs/scratch/00/di98qes3/qX/qX_p2n100
executable=solver/build_supermuc_qx_p6x4x4x4/test_hyp_qx_p20_p40

test -r "$gauge_file"
test -x "$executable"

echo "START $(date -Is)"
echo "JOB=$SLURM_JOB_ID RANKS=$SLURM_NTASKS"
echo "GAUGE=$gauge_file"
echo "GLOBAL=192x64x64x64 GRID=6x4x4x4 LOCAL=32x16x16x16"
echo "HYP_STEPS=1 ALPHA=1,1,0.5 PITER=20,40"

mpiexec -n "$SLURM_NTASKS" "$executable" "$gauge_file"

echo "END $(date -Is)"
