#!/usr/bin/env bash
#SBATCH --job-name=hyp_obs
#SBATCH --account=pn29se
#SBATCH --partition=test
#SBATCH --qos=test
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit
#SBATCH --output=laplace_static/validation/hyp_observable_tiny/hyp_obs_%j.out
#SBATCH --error=laplace_static/validation/hyp_observable_tiny/hyp_obs_%j.err
set -euo pipefail
source env/supermuc.sh
module load slurm_setup
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
export TINY_DIST_ROOT="$PWD/laplace_static/validation/tiny_open8x4x4x4/data"
export NVECS=16 T_SRC=2 T_SINK=6 R_SEP=2
export POINT_X=0 POINT_Y=1 POINT_Z=2
export PROBE_LONG=0 PROBE_TR1=1 PROBE_TR2=2 PRINT_SITE_VALUES=1
case_dir=laplace_static/validation/hyp_observable_tiny
build=solver/build_supermuc_tiny_serial
results="$case_dir/results_$SLURM_JOB_ID"
mkdir "$results"
test -r "$TINY_DIST_ROOT/cnfg/open8x4x4x4n1"
for t in 2 6; do
    test -r "$TINY_DIST_ROOT/disteigvecs/disteigvecs_open8x4x4x4n1_t$t.h5"
done
echo "START $(date -Is) JOB=$SLURM_JOB_ID"
mpiexec -n 1 "$build/clover_baseline" > "$results/baseline.out"
export HYP_MODE=0 HYP_PITER=40
mpiexec -n 1 "$build/clover_hyp" > "$results/off.out"
export HYP_MODE=1
mpiexec -n 1 "$build/clover_hyp" > "$results/zero.out"
export HYP_MODE=2
mpiexec -n 1 "$build/clover_hyp" > "$results/hyp40.out"
export HYP_PITER=80
mpiexec -n 1 "$build/clover_hyp" > "$results/hyp80.out"
python3 "$case_dir/compare.py" "$results"
echo "END $(date -Is)"
