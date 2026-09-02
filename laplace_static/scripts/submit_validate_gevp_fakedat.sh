#!/bin/bash
#SBATCH --job-name=fake_gevp_val
#SBATCH --partition=compute2011
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=laplace_static/validation/fake_gevp/slurm-%j.out
#SBATCH --error=laplace_static/validation/fake_gevp/slurm-%j.err

set -euo pipefail

repo_root="${SLURM_SUBMIT_DIR:-$(pwd)}"
input_path="${FAKE_GEVP_INPUT:-${repo_root}/laplace_static/validation/fake_gevp/data/fakedat.mat}"
output_dir="${FAKE_GEVP_OUTPUT_DIR:-${repo_root}/laplace_static/validation/fake_gevp/results}"

if test ! -f "${repo_root}/laplace_static/scripts/validate_gevp_fakedat.py"; then
  echo "ERROR: submit this job from the Masterarbeit worktree root" >&2
  exit 2
fi

cd "${repo_root}"
mkdir -p "${output_dir}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export PYTHONDONTWRITEBYTECODE=1

python3 -B laplace_static/scripts/analyze_laplace_gaussian_gevp.py --self-test
python3 -B laplace_static/scripts/validate_gevp_fakedat.py \
  --input "${input_path}" \
  --output-dir "${output_dir}"
