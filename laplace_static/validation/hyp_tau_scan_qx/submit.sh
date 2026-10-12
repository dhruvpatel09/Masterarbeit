#!/usr/bin/env bash
#SBATCH --job-name=qx_hyptau
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
#SBATCH --output=laplace_static/validation/hyp_tau_scan_qx/qx_hyptau_%j.out
#SBATCH --error=laplace_static/validation/hyp_tau_scan_qx/qx_hyptau_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export QX_ROOT=/hppfs/scratch/00/di98qes3/qX
export EVROOT=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/eigenvectors
export GX=0 GY=0 GZ=0 PRINT_ALL=1 HYP_MODE=2 HYP_PITER=40
case_dir=laplace_static/validation/hyp_tau_scan_qx
original=laplace_static/validation/hyp_profile_scan_qx
executable="$original/build/test_profile_scan_qx"
reference="$original/results_5903817/hyp40.out"
sha256sum -c "$case_dir/package.sha256"
sha256sum -c "$case_dir/inputs.sha256"
sha256sum -c "$original/build/manifest.sha256"
test -x "$executable"
test -r "$reference"
test -r "$QX_ROOT/qX_p2n100"
for t in 95 96 97 99 100 101; do test -s "$EVROOT/qX_p2n100_evec_t$t.h5"; done
results="$case_dir/results_$SLURM_JOB_ID"
mkdir "$results"
cp "$original/build/manifest.sha256" "$results/manifest.sha256"
cp "$case_dir/package.sha256" "$results/package.sha256"
sha256sum "$reference" > "$results/reference.sha256"
echo "START $(date -Is) JOB=$SLURM_JOB_ID MID=98 TAU=2,4,6 R=1,2,4 RANKS=$SLURM_NTASKS"
for tau in 2 4 6; do
    export TAU_LENGTH="$tau" T0="$((98-tau/2))"
    mpiexec -n "$SLURM_NTASKS" "$executable" > "$results/tau$tau.out"
    echo "TAU_${tau}_DONE $(date -Is)"
done
python3 "$case_dir/check_export.py" "$results" "$reference"
echo "END $(date -Is)"
