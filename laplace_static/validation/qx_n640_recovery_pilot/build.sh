#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
cd "$root"
python3 "$here/campaign.py" check-inputs
test ! -e "$here/build"
source "$root/env/supermuc.sh"
export MKL_THREADING_LAYER=SEQUENTIAL
mkdir "$here/build"
python3 "$here/patch_sources.py" "$here/build"
cd "$here/build"
hp="$(command -v h5pcc)"
h5root="$(cd "$(dirname "$hp")/.." && pwd)"
if [[ -d "$h5root/lib" ]]; then h5lib="$h5root/lib"; else h5lib="$h5root/lib64"; fi
wr="$root/dist/distredelemshdf5/c"
ev="$root/dist/disteigvecshdf5/c"
el="$root/dist/disteigvalshdf5/c"
flags=( -O2 -std=gnu11 -Wall -Wextra "-I$here/build" "-I$here" "-I$root/qcd/include"
  "-I$h5root/include" "-I$root/mental/include_mkl_lapacke" "-I$MKLROOT/include"
  "-I$wr/include" "-I$ev/include" "-I$el/include" )
mpicc "${flags[@]}" -c "$wr/src/DistRedelemsHdf5Writer.c" "$ev/src/DistEigvecsHdf5Writer.c" \
  "$ev/src/DistEigvecsHdf5Reader.c" "$el/src/DistEigvalsHdf5Writer.c" mental_recovery.c check_saved.c
libs=( "$root/qcd/lib/libqcd.a" "-L$h5lib" "-Wl,-rpath,$h5lib" -lhdf5_hl -lhdf5
  "-L$MKLROOT/lib/intel64" "-Wl,-rpath,$MKLROOT/lib/intel64" -lmkl_rt -lpthread -ldl -lz -lm )
mpicxx -o mental_recovery mental_recovery.o DistRedelemsHdf5Writer.o DistEigvecsHdf5Writer.o \
  DistEigvecsHdf5Reader.o DistEigvalsHdf5Writer.o "${libs[@]}"
mpicxx -o check_saved check_saved.o DistEigvecsHdf5Reader.o "${libs[@]}"
ldd mental_recovery > linked_generator.txt
ldd check_saved > linked_checker.txt
if grep -q 'not found' linked_generator.txt linked_checker.txt; then
  cat linked_generator.txt linked_checker.txt
  exit 1
fi
module list > modules.txt 2>&1
mpicc -show > compiler.txt
python3 "$here/campaign.py" record-build
printf 'Built isolated n640 selected-time writer and saved-mode checker\n'
