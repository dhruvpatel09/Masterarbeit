#!/usr/bin/env bash
set -euo pipefail
cd /hppfs/work/pn29se/di24gus/Masterarbeit
case_dir="$PWD/laplace_static/validation/qx_saved_modes_n100"
python3 "$case_dir/audit.py" verify-sources
python3 - <<'PY'
import importlib.util
from pathlib import Path
p=Path('laplace_static/validation/qx_projection_spectrum/campaign.py').resolve()
s=importlib.util.spec_from_file_location('qx_audit_build_guard',str(p))
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
m.build_manifest()
PY
test ! -e "$case_dir/build"
source env/supermuc.sh
root="$PWD"
solver_build="$root/solver/build_supermuc_qx_p1x4x4x4"
reader="$root/dist/disteigvecshdf5/c"
preserved="$root/laplace_static/analysis/qx_eigenpair_completion/inputs"
h5pcc_path="$(command -v h5pcc)"
h5root="$(cd "$(dirname "$h5pcc_path")/.." && pwd)"
if [[ -d "$h5root/lib" ]]; then h5lib="$h5root/lib"; else h5lib="$h5root/lib64"; fi
test -f "$h5root/include/hdf5.h"
test -s "$solver_build/libsolver.a"
mkdir "$case_dir/build"
module list > "$case_dir/build/modules.txt" 2>&1
mpicc -show > "$case_dir/build/compiler.txt"
common=( -std=c11 -O2 -fopenmp -Wall -Wextra
    "-I$root/solver/include" "-I$root/qcd/include" "-I$h5root/include"
    "-I$reader/include" "-I$preserved" )
mpicc "${common[@]}" -c "$reader/src/DistEigvecsHdf5Reader.c" -o "$case_dir/build/reader.o"
mpicc "${common[@]}" -c "$case_dir/test_saved_modes.c" -o "$case_dir/build/audit.o"
mpicc -fopenmp "$case_dir/build/audit.o" "$case_dir/build/reader.o" \
    -Wl,--start-group "$solver_build/libsolver.a" "$root/qcd/lib/libqcd.a" -Wl,--end-group \
    -L"$h5lib" -Wl,-rpath,"$h5lib" -lhdf5_hl -lhdf5 \
    -L"$MKLROOT/lib/intel64" -Wl,-rpath,"$MKLROOT/lib/intel64" \
    -lmkl_rt -lpthread -ldl -lz -lm -o "$case_dir/build/test_saved_modes"
ldd "$case_dir/build/test_saved_modes" > "$case_dir/build/ldd.txt"
if grep -q 'not found' "$case_dir/build/ldd.txt"; then
    echo 'A shared library is missing; stopping.' >&2
    exit 1
fi
sha256sum "$case_dir/build/test_saved_modes" "$case_dir/test_saved_modes.c" \
    "$preserved/metrics.h" "$reader/src/DistEigvecsHdf5Reader.c" \
    "$solver_build/libsolver.a" "$root/qcd/lib/libqcd.a" \
    "$solver_build/include/global.h" "$root/solver/include/solver.h" \
    "$root/env/supermuc.sh" > "$case_dir/build/manifest.sha256"
echo "Built $case_dir/build/test_saved_modes"
