#!/usr/bin/env bash
#SBATCH --job-name=qx_p4obs
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
#SBATCH --output=laplace_static/validation/hyp_projection_observable_qx/p4obs_%j.out
#SBATCH --error=laplace_static/validation/hyp_projection_observable_qx/p4obs_%j.err
set -euo pipefail
case_dir=laplace_static/validation/hyp_projection_observable_qx
original=laplace_static/validation/hyp_profile_scan_qx
reference=laplace_static/validation/hyp_tau_scan_qx/results_5904081
executable="$original/build/test_profile_scan_qx"
test "${SLURM_NTASKS:?}" -eq 64
sha256sum -c "$case_dir/package.sha256"
sha256sum -c "$case_dir/inputs.sha256"
sha256sum -c "$original/build/manifest.sha256"
test -x "$executable"
for tau in 2 4 6; do test -s "$reference/tau$tau.out"; done
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export QX_ROOT=/hppfs/scratch/00/di98qes3/qX
export EVROOT=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/eigenvectors
export GX=0 GY=0 GZ=0 PRINT_ALL=1 HYP_MODE=2
test -r "$QX_ROOT/qX_p2n100"
for t in 95 96 97 99 100 101; do test -s "$EVROOT/qX_p2n100_evec_t$t.h5"; done
results="$case_dir/results_$SLURM_JOB_ID"
mkdir "$results"
cp "$original/build/manifest.sha256" "$results/manifest.sha256"
cp "$case_dir/package.sha256" "$results/package.sha256"
sha256sum "$reference/tau2.out" "$reference/tau4.out" "$reference/tau6.out" > "$results/reference.sha256"
echo "START $(date -Is) JOB=$SLURM_JOB_ID CFG=100 MID=98 TAU=2,4,6 R=1,2,4 PITER=4 REF_PITER=40 RANKS=$SLURM_NTASKS"
export HYP_PITER=4
for tau in 2 4 6; do
    export TAU_LENGTH="$tau" T0="$((98-tau/2))"
    mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/p4_tau$tau.out"
    echo "P4_TAU_${tau}_DONE $(date -Is)"
done
# Fresh control confirms the current executable reproduces the saved reference.
export HYP_PITER=40 TAU_LENGTH=4 T0=96
mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/p40_tau4_control.out"
echo "P40_TAU4_CONTROL_DONE $(date -Is)"
python3 "$case_dir/check_export.py" "$results" "$reference"
echo "END $(date -Is)"
