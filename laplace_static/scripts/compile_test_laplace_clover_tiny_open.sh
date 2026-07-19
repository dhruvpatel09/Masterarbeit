#!/bin/bash
set -euo pipefail

cd /home/m2130292/Masterarbeit/laplace_static

module purge
module load compiler/gcc/10.2.1
module load mpi/openmpi/4.1.0-no_ucx
module load mkl/latest

HOMEBASE=/home/m2130292/Masterarbeit
H5ROOT=${HOMEBASE}/extern/hdf5-1.14.6-openmpi410
MKLROOT_LOCAL=${MKLROOT}

QCDINC=${HOMEBASE}/qcd/include
QCDLIB=${HOMEBASE}/qcd/lib

CFLAGS="-O2 -std=c11 -Wall -Wextra -Wpedantic -I${QCDINC} -I${H5ROOT}/include"
LDFLAGS="-L${QCDLIB} -Wl,-rpath,${QCDLIB} -lqcd"
H5LDFLAGS="-L${H5ROOT}/lib -Wl,-rpath,${H5ROOT}/lib -lhdf5_hl -lhdf5"
MKLLDFLAGS="-L${MKLROOT_LOCAL}/lib/intel64 -Wl,-rpath,${MKLROOT_LOCAL}/lib/intel64 -lmkl_rt -lpthread -ldl"

mkdir -p build

rm -f \
  build/test_laplace_clover_tiny_open \
  build/test_laplace_clover_tiny_open.o

mpicc ${CFLAGS} \
  -c src/test_laplace_clover_tiny_open.c \
  -o build/test_laplace_clover_tiny_open.o

mpicxx -o build/test_laplace_clover_tiny_open \
  build/test_laplace_clover_tiny_open.o \
  ${LDFLAGS} \
  ${H5LDFLAGS} \
  ${MKLLDFLAGS} \
  -lz -lm

echo "Built build/test_laplace_clover_tiny_open"
ldd build/test_laplace_clover_tiny_open \
  | egrep "not found|hdf5|mpi|mkl|qcd|z" || true
