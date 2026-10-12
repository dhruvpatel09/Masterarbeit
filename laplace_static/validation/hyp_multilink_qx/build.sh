#!/usr/bin/env bash
set -euo pipefail
cd /hppfs/work/pn29se/di24gus/Masterarbeit
source env/supermuc.sh
root="$PWD"
case_dir="$root/laplace_static/validation/hyp_multilink_qx"
build="$root/solver/build_supermuc_qx_p1x4x4x4"
sha256sum -c "$case_dir/inputs.sha256"
test -s "$build/libsolver.a"
h5pcc_path="$(command -v h5pcc)"
h5root="$(cd "$(dirname "$h5pcc_path")/.." && pwd)"
if [[ -d "$h5root/lib" ]]; then h5lib="$h5root/lib"; else h5lib="$h5root/lib64"; fi
test -f "$h5root/include/hdf5.h"
test -d "$h5lib"
common=( -std=c11 -O2 -fopenmp -Wall -Wextra
    "-I$root/solver/include" "-I$root/qcd/include" "-I$h5root/include" )
libs=( -Wl,--start-group "$build/libsolver.a" "$root/qcd/lib/libqcd.a" -Wl,--end-group
    "-L$h5lib" "-Wl,-rpath,$h5lib" -lhdf5_hl -lhdf5
    "-L$MKLROOT/lib/intel64" "-Wl,-rpath,$MKLROOT/lib/intel64"
    -lmkl_rt -lpthread -ldl -lz -lm )
reader="$root/dist/disteigvecshdf5/c"
mkdir -p "$case_dir/build"
mpicc "${common[@]}" "-I$reader/include" \
    -c "$reader/src/DistEigvecsHdf5Reader.c" -o "$case_dir/build/reader.o"
mpicc "${common[@]}" "-I$reader/include" \
    "$case_dir/test_multilink_qx.c" "$case_dir/build/reader.o" \
    "${libs[@]}" -o "$case_dir/build/test_multilink_qx"
sha256sum "$case_dir/test_multilink_qx.c" "$case_dir/build/test_multilink_qx" \
    "$build/libsolver.a" "$root/qcd/lib/libqcd.a" \
    "$build/include/global.h" > "$case_dir/build/manifest.sha256"
echo "Built $case_dir/build/test_multilink_qx"
