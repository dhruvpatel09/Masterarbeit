#!/usr/bin/env bash
#SBATCH --job-name=qx_tau1
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks=64
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit/laplace_static
#SBATCH --output=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/laplace_static_smoke_20260915/tau_t96_xyz000/slurm/slurm_%j.out
#SBATCH --error=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/laplace_static_smoke_20260915/tau_t96_xyz000/slurm/slurm_%j.err

set -euo pipefail

ROOT="/hppfs/work/pn29se/di24gus/Masterarbeit"

source "${ROOT}/env/supermuc.sh"
module load slurm_setup

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export OMP_DYNAMIC=FALSE
export MKL_THREADING_LAYER=SEQUENTIAL

export QX_ROOT="/hppfs/scratch/00/di98qes3/qX"
export EVROOT="/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/eigenvectors"

export T0=96
export GX=0
export GY=0
export GZ=0
export PRINT_ALL=0

EXE="${ROOT}/laplace_static/build/test_tau_qX_one_point"

echo "START: $(date)"
echo "JOB_ID=${SLURM_JOB_ID:-unknown}"
echo "NODES=${SLURM_JOB_NUM_NODES:-unknown}"
echo "NTASKS=${SLURM_NTASKS:-unknown}"
echo "HOST=$(hostname)"
echo "EXE=${EXE}"
echo "QX_ROOT=${QX_ROOT}"
echo "EVROOT=${EVROOT}"
echo "T0=${T0}"
echo "XYZ=${GX},${GY},${GZ}"

test -x "${EXE}"
test -r "${QX_ROOT}/qX_p2n100"
test -r "${EVROOT}/qX_p2n100_evec_t96.h5"
test -r "${EVROOT}/qX_p2n100_evec_t97.h5"

mpiexec -n "${SLURM_NTASKS}" \
    "${EXE}"

echo "END: $(date)"
