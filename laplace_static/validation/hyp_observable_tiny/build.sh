#!/usr/bin/env bash
set -euo pipefail
cd /hppfs/work/pn29se/di24gus/Masterarbeit
source env/supermuc.sh
root="$PWD"
case_dir="$root/laplace_static/validation/hyp_observable_tiny"
sha256sum -c "$case_dir/baseline.sha256"
build="$root/solver/build_supermuc_tiny_serial"
if [[ -e "$build" ]]; then
    echo "Build directory already exists: $build; stopping." >&2
    exit 1
fi
python3 - <<'PY'
from pathlib import Path
import shutil,re
src=Path('solver/build_supermuc_tiny')
dst=Path('solver/build_supermuc_tiny_serial')
text=(src/'include/global.h').read_text()
values={f'NPROC{i}':1 for i in range(4)}
values.update({f'L{i}':(8 if i==0 else 4) for i in range(4)})
values.update({f'L{i}_TRD':(8 if i==0 else 4) for i in range(4)})
for key,val in values.items():
    text,n=re.subn(r'^#define\s+'+key+r'\s+\d+\s*$',f'#define {key} {val}',text,flags=re.M)
    if n!=1: raise SystemExit('Unexpected header definition: '+key)
dst.mkdir()
shutil.copy2(src/'Makefile',dst/'Makefile')
shutil.copytree(src/'include',dst/'include')
(dst/'include/global.h').write_text(text)
PY
(
    cd "$build"
    make -j2 mklib GCC=mpicc CC=mpicc \
        "MDIR=$root/openQCD-2.4.2/modules" "MAINDIR=$root/solver" \
        "QCD_HOME=$root/qcd" \
        "INCPATH=./include $root/solver/include $root/solver/smear $root/qcd/include" \
        'CFLAGS=-std=c99 -O2 -g -fopenmp -Wall' 2>&1 | tee build.log
)
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
mpicc "${common[@]}" "$root/laplace_static/src/test_laplace_clover_tiny_open.c" \
    "${libs[@]}" -o "$build/clover_baseline"
mpicc "${common[@]}" "$case_dir/test_laplace_clover_hyp_tiny.c" \
    "${libs[@]}" -o "$build/clover_hyp"
echo "Built baseline and HYP integration executables in $build"
