#!/usr/bin/env bash
#SBATCH --job-name=qx_hyp_time
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
#SBATCH --output=laplace_static/validation/hyp_projection_timeslices_qx/projection_time_%j.out
#SBATCH --error=laplace_static/validation/hyp_projection_timeslices_qx/projection_time_%j.err
set -euo pipefail
case_dir="$PWD/laplace_static/validation/hyp_projection_timeslices_qx"
test "${SLURM_JOB_NUM_NODES:?}" -eq 8
test "${SLURM_NTASKS:?}" -eq 384
sha256sum -c "$case_dir/package.sha256"
sha256sum -c "$case_dir/inputs.sha256"
sha256sum -c "$case_dir/build/manifest.sha256"
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
gauge=/hppfs/scratch/00/di98qes3/qX/qX_p2n100
test -r "$gauge"
result="$case_dir/results_${SLURM_JOB_ID}"
mkdir "$result"
cp "$case_dir/build/manifest.sha256" "$result/manifest.sha256"
echo "START $(date -Is) JOB=$SLURM_JOB_ID CFG=100 GRID=6x4x4x4 RANKS=$SLURM_NTASKS HYP_STEPS=1 ALPHA=1,1,0.5 PITER=4,40 REPEAT4=YES"
mpiexec -n "$SLURM_NTASKS" \
    "$case_dir/build/test_projection_timeslices_qx" "$gauge" \
    | tee "$result/diagnostics.out"
python3 "$case_dir/check_export.py" "$result/diagnostics.out" "$result"
echo "END $(date -Is) RESULTS=$result"
