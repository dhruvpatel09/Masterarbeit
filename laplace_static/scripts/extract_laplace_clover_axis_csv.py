#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd

if len(sys.argv) != 5:
    raise SystemExit("Usage: extract_laplace_clover_axis_csv.py AXIS TAU JOBID OUTCSV")

axis = int(sys.argv[1])
tau = int(sys.argv[2])
jobid = sys.argv[3]
out = Path(sys.argv[4])

logs = sorted(Path("logs").glob(f"lap_ax{axis}_R4t{tau}_n1-100_{jobid}_*.out"))

cols = [
    "cfg","r","tau","axis","dx","dy","dz","Nsrc",
    "L_Re","L_Im",
    "E2_vac","B2_vac","S_vac",
    "LE2_Re","LE2_Im","LB2_Re","LB2_Im","LS_Re","LS_Im",
    "rho_E2","rho_B2","rho_S",
]

rows = []

for log in logs:
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("PROBE "):
            parts = line.split()[1:]
            if len(parts) != len(cols):
                raise RuntimeError(f"Bad PROBE column count in {log}: {len(parts)}")
            rows.append(parts)

df = pd.DataFrame(rows, columns=cols)

for c in cols:
    df[c] = pd.to_numeric(df[c], errors="raise")

df = df.sort_values(["cfg", "axis", "dx", "dy"]).reset_index(drop=True)

assert df.shape == (100 * 117, len(cols)), df.shape
assert sorted(df["cfg"].unique()) == list(range(1, 101))
assert sorted(df["axis"].unique()) == [axis]
assert sorted(df["dx"].unique()) == list(range(-2, 7))
assert sorted(df["dy"].unique()) == list(range(-6, 7))
assert (df["Nsrc"] == 48 * 24**3).all()
assert df.notna().all().all()

df.to_csv(out, index=False, lineterminator="\n")

print("PASS wrote:", out)
print("shape =", df.shape)
print(df[["rho_E2","rho_B2","rho_S"]].describe())
