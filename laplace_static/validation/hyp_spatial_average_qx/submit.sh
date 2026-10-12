#!/usr/bin/env bash
#SBATCH --job-name=qx_hypavg
#SBATCH --account=pn29se
#SBATCH --partition=test
#SBATCH --qos=test
#SBATCH --time=00:30:00
#SBATCH --nodes=8
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_spatial_average_qx/qx_hypavg_%j.out
#SBATCH --error=laplace_static/validation/hyp_spatial_average_qx/qx_hypavg_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export QX_ROOT=/hppfs/scratch/00/di98qes3/qX
export EVROOT=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/eigenvectors
export T0=96 TAU_LENGTH=4 GX=0 GY=0 GZ=0 PRINT_ALL=1
case_dir=laplace_static/validation/hyp_spatial_average_qx
executable="$case_dir/build/test_spatial_average_qx"
reference=laplace_static/validation/hyp_multilink_qx/results_5903809
for name in baseline off zero hyp40 hyp80; do
    test -r "$reference/$name.out"
done
results="$case_dir/results_$SLURM_JOB_ID"
mkdir "$results"
test -r "$QX_ROOT/qX_p2n100"
for t in 96 100; do
    file="$EVROOT/qX_p2n100_evec_t$t.h5"
    if [[ ! -r "$file" ]]; then echo "Missing eigenvector file: $file" >&2; exit 1; fi
done
test -x "$executable"
sha256sum -c "$case_dir/build/manifest.sha256"
cp "$case_dir/build/manifest.sha256" "$results/manifest.sha256"
sha256sum "$reference/baseline.out" "$reference/off.out" "$reference/zero.out" \
    "$reference/hyp40.out" "$reference/hyp80.out" > "$results/reference.sha256"
echo "START $(date -Is) JOB=$SLURM_JOB_ID GRID=1x4x4x4 RANKS=$SLURM_NTASKS"
export HYP_MODE=3 HYP_PITER=40
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/baseline.out"
echo "BASELINE_DONE $(date -Is)"
export HYP_MODE=0 HYP_PITER=40
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/off.out"
echo "OFF_DONE $(date -Is)"
export HYP_MODE=1
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/zero.out"
echo "ZERO_DONE $(date -Is)"
export HYP_MODE=2
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/hyp40.out"
echo "HYP40_DONE $(date -Is)"
export HYP_PITER=80
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/hyp80.out"
echo "HYP80_DONE $(date -Is)"
python3 "$case_dir/compare.py" "$results" "$reference"
echo "END $(date -Is)"
