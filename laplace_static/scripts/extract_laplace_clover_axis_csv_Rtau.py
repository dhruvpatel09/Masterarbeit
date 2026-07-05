#!/usr/bin/env python3

from pathlib import Path
import sys
import pandas as pd

if len(sys.argv) != 6:
    raise SystemExit("Usage: extract_laplace_clover_axis_csv_Rtau.py R TAU AXIS JOBID OUTCSV")

R = int(sys.argv[1])
tau = int(sys.argv[2])
axis = int(sys.argv[3])
jobid = sys.argv[4]
out = Path(sys.argv[5])

logs = sorted(Path("logs").glob(f"lap_ax{axis}_R{R}t{tau}_n1-100_{jobid}_*.out"))

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

nprobe = df.groupby("cfg").size().iloc[0]

assert len(logs) == 100, len(logs)
assert df.shape == (100 * nprobe, len(cols)), df.shape
assert sorted(df["cfg"].unique()) == list(range(1, 101))
assert sorted(df["r"].unique()) == [R]
assert sorted(df["tau"].unique()) == [tau]
assert sorted(df["axis"].unique()) == [axis]
assert (df.groupby("cfg").size() == nprobe).all()
assert (df["Nsrc"] == 48 * 24**3).all()
assert df.notna().all().all()

df.to_csv(out, index=False, lineterminator="\n")

print("PASS wrote:", out)
print("shape =", df.shape)
print("nprobe per cfg =", nprobe)
print("dx range =", df["dx"].min(), df["dx"].max())
print("dy range =", df["dy"].min(), df["dy"].max())
print()
print(df[["rho_E2","rho_B2","rho_S"]].describe())
