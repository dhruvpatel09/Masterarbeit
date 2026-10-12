#!/usr/bin/env bash
#SBATCH --job-name=qx_hypobs
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
#SBATCH --output=laplace_static/validation/hyp_observable_qx/qx_hypobs_%j.out
#SBATCH --error=laplace_static/validation/hyp_observable_qx/qx_hypobs_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export QX_ROOT=/hppfs/scratch/00/di98qes3/qX
export EVROOT=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/eigenvectors
export T0=96 GX=0 GY=0 GZ=0 PRINT_ALL=1
case_dir=laplace_static/validation/hyp_observable_qx
build=solver/build_supermuc_qx_p1x4x4x4
results="$case_dir/results_$SLURM_JOB_ID"
mkdir "$results"
test -r "$QX_ROOT/qX_p2n100"
for t in 96 97; do test -r "$EVROOT/qX_p2n100_evec_t$t.h5"; done
echo "START $(date -Is) JOB=$SLURM_JOB_ID GRID=1x4x4x4 RANKS=$SLURM_NTASKS"
mpiexec -n "$SLURM_NTASKS" "$build/tau_baseline" > "$results/baseline.out"
echo "BASELINE_DONE $(date -Is)"
export HYP_MODE=0 HYP_PITER=40
mpiexec -n "$SLURM_NTASKS" "$build/tau_hyp" > "$results/off.out"
echo "OFF_DONE $(date -Is)"
export HYP_MODE=1
mpiexec -n "$SLURM_NTASKS" "$build/tau_hyp" > "$results/zero.out"
echo "ZERO_DONE $(date -Is)"
export HYP_MODE=2
mpiexec -n "$SLURM_NTASKS" "$build/tau_hyp" > "$results/hyp40.out"
echo "HYP40_DONE $(date -Is)"
export HYP_PITER=80
mpiexec -n "$SLURM_NTASKS" "$build/tau_hyp" > "$results/hyp80.out"
echo "HYP80_DONE $(date -Is)"
python3 "$case_dir/compare.py" "$results"
echo "END $(date -Is)"
