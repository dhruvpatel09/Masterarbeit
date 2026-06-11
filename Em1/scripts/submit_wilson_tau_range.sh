#!/bin/bash
set -euo pipefail

FIRST="$1"
LAST="$2"
TAUS="$3"
MODE="${4:-new}"   # new or append

mkdir -p inputs slurm slurm_logs

prev_jobid=""

for tau in $TAUS; do
    tag=$(printf "%02d" "$tau")

    infile="inputs/wilson_tau_tau${tag}_${FIRST}_${LAST}.in"
    jobfile="slurm/submit_wilson_tau_tau${tag}_${FIRST}_${LAST}.sh"

    cp wilson_tau.in "$infile"
    sed -i "s/^output.*/output  Em1p4_wilson_mu0_tau${tag}/" "$infile"
    sed -i "s/^tau.*/tau         ${tau}/" "$infile"
    sed -i "s/^first.*/first       ${FIRST}/" "$infile"
    sed -i "s/^last.*/last        ${LAST}/" "$infile"

    runline="mpirun ./wilson_tau -i $infile"
    if [[ "$MODE" == "append" ]]; then
        runline="mpirun ./wilson_tau -i $infile -a"
    fi

    cat > "$jobfile" <<EOF
#!/bin/bash
#SBATCH --job-name=EM1_t${tag}_${FIRST}_${LAST}
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=compute2011
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

module load compiler/gcc/10.2.1
module load mpi/openmpi/4.1.0-no_ucx

$runline
EOF

    chmod +x "$jobfile"

    if [[ -z "$prev_jobid" ]]; then
        jobid=$(sbatch "$jobfile" | awk '{print $4}')
    else
        jobid=$(sbatch --dependency=afterok:${prev_jobid} "$jobfile" | awk '{print $4}')
    fi

    echo "Submitted tau=${tau}, configs ${FIRST}-${LAST}, mode=${MODE}, job=${jobid}"
    prev_jobid="$jobid"
done