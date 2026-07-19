#!/bin/bash
#SBATCH --job-name=laplace_clover_tiny
#SBATCH --partition=compute2011
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail

cd /home/m2130292/Masterarbeit/laplace_static

module purge
module load compiler/gcc/10.2.1
module load mpi/openmpi/4.1.0-no_ucx
module load mkl/latest

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_THREADING_LAYER=sequential

export TINY_DIST_ROOT="${TINY_DIST_ROOT:-/home/m2130292/Masterarbeit/laplace_static/validation/tiny_open8x4x4x4/data}"
export NVECS_LIST="${NVECS_LIST:-16}"
export T_SRC="${T_SRC:-2}"
export T_SINK="${T_SINK:-6}"
export R_SEP="${R_SEP:-2}"
export POINT_X="${POINT_X:-0}"
export POINT_Y="${POINT_Y:-1}"
export POINT_Z="${POINT_Z:-2}"
export PROBE_LONG="${PROBE_LONG:-0}"
export PROBE_TR1="${PROBE_TR1:-1}"
export PROBE_TR2="${PROBE_TR2:-2}"
export PRINT_SITE_VALUES="${PRINT_SITE_VALUES:-0}"

results_dir="validation/tiny_open8x4x4x4/results"
mkdir -p "${results_dir}"

echo "HOST=$(hostname)"
echo "DATE=$(date)"
echo "JOB_ID=${SLURM_JOB_ID:-none}"
echo "TINY_DIST_ROOT=${TINY_DIST_ROOT}"
echo "NVECS_LIST=${NVECS_LIST}"
echo "T_SRC=${T_SRC}"
echo "T_SINK=${T_SINK}"
echo "R_SEP=${R_SEP}"
echo "POINT=${POINT_X},${POINT_Y},${POINT_Z}"
echo "PROBE_LOCAL=${PROBE_LONG},${PROBE_TR1},${PROBE_TR2}"
echo "PRINT_SITE_VALUES=${PRINT_SITE_VALUES}"

for nvecs in ${NVECS_LIST}; do
  export NVECS="${nvecs}"
  output="${results_dir}/laplace_clover_tiny_open_cfg1_t${T_SRC}-${T_SINK}_r${R_SEP}_Nv${NVECS}_probe${PROBE_LONG}_${PROBE_TR1}_${PROBE_TR2}.txt"

  echo "Running NVECS=${NVECS} -> ${output}"

  mpirun \
    -x TINY_DIST_ROOT \
    -x NVECS \
    -x T_SRC \
    -x T_SINK \
    -x R_SEP \
    -x POINT_X \
    -x POINT_Y \
    -x POINT_Z \
    -x PROBE_LONG \
    -x PROBE_TR1 \
    -x PROBE_TR2 \
    -x PRINT_SITE_VALUES \
    -np 1 \
    ./build/test_laplace_clover_tiny_open \
    | tee "${output}"
done
