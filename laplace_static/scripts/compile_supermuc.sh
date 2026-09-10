#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
LAPROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT="$(cd "${LAPROOT}/.." && pwd)"

source "${ROOT}/env/supermuc.sh"

H5PCC="$(command -v h5pcc)"
H5ROOT="$(cd "$(dirname "${H5PCC}")/.." && pwd)"

if [[ -d "${H5ROOT}/lib" ]]; then
    H5LIB="${H5ROOT}/lib"
elif [[ -d "${H5ROOT}/lib64" ]]; then
    H5LIB="${H5ROOT}/lib64"
else
    echo "Could not locate HDF5 library directory below ${H5ROOT}" >&2
    exit 1
fi

QCDINC="${ROOT}/qcd/include"
QCDLIB="${ROOT}/qcd/lib"
EVECH5="${ROOT}/dist/disteigvecshdf5/c"

MKLLIB="${MKLROOT}/lib/intel64"

mkdir -p "${LAPROOT}/build"

COMMON_INC=(
    "-I${QCDINC}"
    "-I${H5ROOT}/include"
)

COMMON_LINK=(
    "-L${QCDLIB}"
    "-Wl,-rpath,${QCDLIB}"
    -lqcd
    "-L${H5LIB}"
    "-Wl,-rpath,${H5LIB}"
    -lhdf5_hl
    -lhdf5
    "-L${MKLLIB}"
    "-Wl,-rpath,${MKLLIB}"
    -lmkl_rt
    -lpthread
    -ldl
    -lz
    -lm
)

build_reader_target() {
    local src="$1"
    local exe="$2"
    local tag="$3"

    local reader_obj="${LAPROOT}/build/DistEigvecsHdf5Reader_${tag}.o"
    local src_obj="${LAPROOT}/build/${exe}.o"
    local binary="${LAPROOT}/build/${exe}"

    rm -f "${reader_obj}" "${src_obj}" "${binary}"

    mpicc -O2 -Wall \
        "${COMMON_INC[@]}" \
        "-I${EVECH5}/include" \
        -c "${EVECH5}/src/DistEigvecsHdf5Reader.c" \
        -o "${reader_obj}"

    mpicc -O2 -Wall \
        "${COMMON_INC[@]}" \
        "-I${EVECH5}/include" \
        -c "${LAPROOT}/src/${src}" \
        -o "${src_obj}"

    mpicxx -o "${binary}" \
        "${src_obj}" \
        "${reader_obj}" \
        "${COMMON_LINK[@]}"

    echo "Built ${binary}"
}

build_tiny() {
    local exe="test_laplace_clover_tiny_open"
    local src_obj="${LAPROOT}/build/${exe}.o"
    local binary="${LAPROOT}/build/${exe}"

    rm -f "${src_obj}" "${binary}"

    mpicc -O2 -std=c11 -Wall -Wextra -Wpedantic \
        "${COMMON_INC[@]}" \
        -c "${LAPROOT}/src/${exe}.c" \
        -o "${src_obj}"

    mpicxx -o "${binary}" \
        "${src_obj}" \
        "${COMMON_LINK[@]}"

    echo "Built ${binary}"
}

echo "=== SuperMUC-NG laplace_static build ==="
echo "ROOT=${ROOT}"
echo "H5ROOT=${H5ROOT}"
echo "MKLROOT=${MKLROOT}"
echo "MKL_THREADING_LAYER=${MKL_THREADING_LAYER}"
echo

build_reader_target \
    test_laplace_RT_weighted_scan.c \
    test_laplace_RT_weighted_scan \
    RT_weighted

build_reader_target \
    test_laplace_clover_probe_R4tau4_plane.c \
    test_laplace_clover_probe_R4tau4_plane \
    clover

build_tiny

echo
echo "=== Runtime-library check ==="

for exe in \
    test_laplace_RT_weighted_scan \
    test_laplace_clover_probe_R4tau4_plane \
    test_laplace_clover_tiny_open
do
    echo "--- ${exe} ---"
    ldd "${LAPROOT}/build/${exe}" \
        | grep -E 'not found|hdf5|mpi|mkl|z' || true
done
