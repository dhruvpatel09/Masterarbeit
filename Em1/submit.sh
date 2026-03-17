#!/bin/bash
# submit.sh
#SBATCH --job-name=EM1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=compute2011
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# compile openQCD
module load compiler/gcc/10.2.1 
module load mpi/openmpi/4.1.0-no_ucx

# run the program
# mpirun ./polya -i polya.in
mpirun ./polya -i polya.in -a