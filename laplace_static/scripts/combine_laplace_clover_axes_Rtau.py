#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path

if len(sys.argv) != 3:
    raise SystemExit("Usage: combine_laplace_clover_axes_Rtau.py R TAU")

R = int(sys.argv[1])
tau = int(sys.argv[2])

files = [
    f"results/laplace_clover_probe_R{R}tau{tau}_axis0_t0avg48_n1-100.csv",
    f"results/laplace_clover_probe_R{R}tau{tau}_axis1_t0avg48_n1-100.csv",
    f"results/laplace_clover_probe_R{R}tau{tau}_axis2_t0avg48_n1-100.csv",
]

out = Path(f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100.csv")

dfs = [pd.read_csv(f) for f in files]
df = pd.concat(dfs, ignore_index=True)
df = df.sort_values(["cfg", "axis", "dx", "dy"]).reset_index(drop=True)

nprobe = df.groupby(["cfg", "axis"]).size().iloc[0]

assert df.shape[0] == 3 * 100 * nprobe, df.shape
assert sorted(df["cfg"].unique()) == list(range(1, 101))
assert sorted(df["r"].unique()) == [R]
assert sorted(df["tau"].unique()) == [tau]
assert sorted(df["axis"].unique()) == [0, 1, 2]
assert (df["Nsrc"] == 48 * 24**3).all()
assert df.notna().all().all()

check = df.groupby(["cfg", "r", "tau", "dx", "dy", "dz"])["axis"].nunique()
assert check.min() == 3 and check.max() == 3, check.describe()

df.to_csv(out, index=False, lineterminator="\n")

print("PASS wrote:", out)
print("shape =", df.shape)
print("nprobe per cfg/axis =", nprobe)
print()
print(df.groupby("axis")[["rho_E2","rho_B2","rho_S"]].mean())
