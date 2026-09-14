#!/usr/bin/env bash
#SBATCH --job-name=qx_n100_ev100
#SBATCH --account=pn29se
#SBATCH --partition=micro
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=48
#SBATCH --ntasks=768
#SBATCH --time=02:00:00
#SBATCH --export=NONE
#SBATCH --chdir=/hppfs/work/pn29se/di24gus/Masterarbeit/mental
#SBATCH --output=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/slurm/slurm_%j.out
#SBATCH --error=/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/slurm/slurm_%j.err

set -euo pipefail

source "/hppfs/work/pn29se/di24gus/Masterarbeit/env/supermuc.sh"
module load slurm_setup

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export OMP_DYNAMIC=FALSE
export MKL_THREADING_LAYER=SEQUENTIAL


INPUT="/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/mental_qX_n100_Nv100_evec.in"
MAINLOG="/hppfs/scratch/0E/di24gus/Masterarbeit/qX/mental_runs/n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/log/mental_qX_n100_Nv100_evec.log"

echo "START: $(date)"
echo "JOB_ID=${SLURM_JOB_ID:-unknown}"
echo "NODES=${SLURM_JOB_NUM_NODES:-unknown}"
echo "NTASKS=${SLURM_NTASKS:-unknown}"
echo "HOST=$(hostname)"
echo "INPUT=$INPUT"
echo "MAINLOG=$MAINLOG"

mpiexec -n "${SLURM_NTASKS}"     ./mental     -i "$INPUT"     -o "$MAINLOG"

echo "END: $(date)"
