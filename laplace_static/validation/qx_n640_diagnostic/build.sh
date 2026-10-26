#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$case_dir/../../.." && pwd)"
cd "$root"
python3 "$case_dir/campaign.py" check-inputs
source "$root/env/supermuc.sh"
export MKL_THREADING_LAYER=SEQUENTIAL
build="$case_dir/build"
test ! -e "$build"
mkdir "$build"
python3 "$case_dir/instrument_lanczos.py" "$build/qcd_lanczos_trace.c"
python3 "$case_dir/campaign.py" config-header > "$build/diagnostic_config.h"
cd "$build"
flags=( -O2 -std=c11 -Wall -Wextra "-I$build" "-I$root/qcd/include"
        "-I$root/mental/include_mkl_lapacke" "-I$MKLROOT/include" )
mpicc "${flags[@]}" -c "$build/qcd_lanczos_trace.c" "$case_dir/test_n640.c"
mpicc "${flags[@]}" '-DQX_VARIANT="library"' -c "$case_dir/test_n640.c" -o test_n640_reference.o
mpicxx -o test_n640_reference test_n640_reference.o "$root/qcd/lib/libqcd.a" \
    "-L$MKLROOT/lib/intel64" "-Wl,-rpath,$MKLROOT/lib/intel64" \
    -lmkl_rt -lpthread -ldl -lz -lm
mpicxx -o test_n640 test_n640.o qcd_lanczos_trace.o "$root/qcd/lib/libqcd.a" \
    "-Wl,-Map,$build/link.map" \
    "-L$MKLROOT/lib/intel64" "-Wl,-rpath,$MKLROOT/lib/intel64" \
    -lmkl_rt -lpthread -ldl -lz -lm
ldd ./test_n640 > linked_libraries.txt
if grep -q 'not found' linked_libraries.txt; then cat linked_libraries.txt; exit 1; fi
module list > modules.txt 2>&1
mpicc -show > compiler.txt
python3 "$case_dir/campaign.py" record-build
printf 'Built %s\n' "$build/test_n640"
