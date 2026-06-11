#!/bin/bash
set -euo pipefail

cd /home/m2130292/Masterarbeit

OUT=laplace_static/logs/qcd_types_$(date +%Y%m%d_%H%M%S).log

{
  echo "=== qcd_complex_16 ==="
  sed -n '110,130p' qcd/include/qcd.h

  echo
  echo "=== qcd_geometry ==="
  sed -n '135,255p' qcd/include/qcd.h

  echo
  echo "=== qcd_spinorComponent3d ==="
  sed -n '292,306p' qcd/include/qcd.h

  echo
  echo "=== qcd_gaugeField ==="
  sed -n '327,341p' qcd/include/qcd.h

  echo
  echo "=== qcd_gaugeField3d ==="
  sed -n '360,373p' qcd/include/qcd.h

  echo
  echo "=== qcd_3x3Ops useful macros/functions ==="
  grep -n -m 80 "MUL\|ADJOINT\|trace\|copy\|zero\|unit" qcd/include/qcd_3x3Ops.h || true

  echo
  echo "=== qcd complex ops/macros ==="
  grep -n -m 80 "qcd_CMUL\|qcd_CADD\|qcd_CONJ\|qcd_CSCALE\|qcd_CSUB" qcd/include/qcd_COps.h || true
} > "$OUT" 2>&1

echo "Wrote: $OUT"
tail -80 "$OUT"
