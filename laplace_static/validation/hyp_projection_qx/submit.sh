#!/usr/bin/env bash
#SBATCH --job-name=qx_hypproj4
#SBATCH --account=pn29se
#SBATCH --partition=test
#SBATCH --qos=test
#SBATCH --nodes=8
#SBATCH --ntasks=384
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_projection_qx/projection_%j.out
#SBATCH --error=laplace_static/validation/hyp_projection_qx/projection_%j.err
set -euo pipefail
case_dir="$PWD/laplace_static/validation/hyp_projection_qx"
sha256sum -c qx_projection_scaling_package.sha256
sha256sum -c "$case_dir/inputs.sha256"
sha256sum -c "$case_dir/build/manifest.sha256"
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
gauge=/hppfs/scratch/00/di98qes3/qX/qX_p2n100
test -r "$gauge"
echo "START $(date -Is) JOB=$SLURM_JOB_ID CFG=100 RANKS=$SLURM_NTASKS HYP_STEPS=1 ALPHA=1,1,0.5 PITER=4,40 REPEAT4=YES"
mpiexec -n "$SLURM_NTASKS" "$case_dir/build/test_projection_qx" "$gauge"
echo "END $(date -Is)"
