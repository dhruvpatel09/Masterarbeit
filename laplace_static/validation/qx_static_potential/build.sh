#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
cd "$root"
python3 "$here/campaign.py" check-inputs
test ! -e "$here/build"
# Preserve the library and environment used by the validated qX profile build.
sha256sum -c laplace_static/validation/qx_clover_trace/build/manifest.sha256
source env/supermuc.sh
mkdir "$here/build"
python3 "$here/generate_helpers.py" "$here/build"
reader="$root/dist/disteigvecshdf5/c"
solver_build="$root/solver/build_supermuc_qx_p1x4x4x4"
h5pcc_path="$(command -v h5pcc)"
h5root="$(cd "$(dirname "$h5pcc_path")/.." && pwd)"
if [[ -d "$h5root/lib" ]]; then h5lib="$h5root/lib"; else h5lib="$h5root/lib64"; fi
old="$root/laplace_static/validation/hyp_projection_configs_qx_200_230"
module list > "$here/build/modules.txt" 2>&1
mpicc -show > "$here/build/compiler.txt"
flags=( -std=c11 -O2 -fopenmp -Wall -Wextra
    "-I$here/build" "-I$root/solver/include" "-I$root/qcd/include"
    "-I$reader/include" "-I$h5root/include" "-I$old" )
mpicc "${flags[@]}" -c "$reader/src/DistEigvecsHdf5Reader.c" -o "$here/build/reader.o"
for target in scan audit_selected; do
    if [[ "$target" == scan ]]; then src="$here/scan.c"; else src="$here/build/audit_selected.c"; fi
    mpicc "${flags[@]}" -c "$src" -o "$here/build/$target.o"
    mpicc -fopenmp "$here/build/$target.o" "$here/build/reader.o" \
        -Wl,--start-group "$solver_build/libsolver.a" "$root/qcd/lib/libqcd.a" -Wl,--end-group \
        -L"$h5lib" -Wl,-rpath,"$h5lib" -lhdf5_hl -lhdf5 \
        -L"$MKLROOT/lib/intel64" -Wl,-rpath,"$MKLROOT/lib/intel64" \
        -lmkl_rt -lpthread -ldl -lz -lm -o "$here/build/$target"
    ldd "$here/build/$target" > "$here/build/linked_$target.txt"
    if grep -q 'not found' "$here/build/linked_$target.txt"; then
        cat "$here/build/linked_$target.txt"
        exit 1
    fi
done
python3 "$here/campaign.py" record-build
printf 'Built fixed-source static correlator scan and selected-time saved-mode audit\n'
