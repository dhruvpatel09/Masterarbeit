#!/usr/bin/env bash
set -euo pipefail
root=/hppfs/work/pn29se/di24gus/Masterarbeit
cd "$root"
case_dir="$root/laplace_static/validation/hyp_projection_timeslices_qx"
sha256sum -c "$case_dir/package.sha256"
sha256sum -c "$case_dir/inputs.sha256"
test ! -e "$case_dir/build"
source env/supermuc.sh
mkdir "$case_dir/build"
mpicc -std=c99 -O2 -fopenmp -Wall -Wextra \
    -I"$root/solver/include" -I"$root/qcd/include" \
    "$case_dir/test_projection_timeslices_qx.c" \
    -Wl,--start-group \
    "$root/solver/build_supermuc_qx_p6x4x4x4/libsolver.a" "$root/qcd/lib/libqcd.a" \
    -Wl,--end-group \
    -L"$MKLROOT/lib/intel64" -Wl,-rpath,"$MKLROOT/lib/intel64" \
    -lmkl_rt -lpthread -ldl -lm -o "$case_dir/build/test_projection_timeslices_qx"
sha256sum "$case_dir/build/test_projection_timeslices_qx" \
    "$case_dir/test_projection_timeslices_qx.c" "$case_dir/check_export.py" \
    "$case_dir/build.sh" "$case_dir/submit.sh" \
    "$root/solver/build_supermuc_qx_p6x4x4x4/libsolver.a" "$root/qcd/lib/libqcd.a" \
    "$root/solver/build_supermuc_qx_p6x4x4x4/include/global.h" \
    > "$case_dir/build/manifest.sha256"
echo "Built $case_dir/build/test_projection_timeslices_qx"
