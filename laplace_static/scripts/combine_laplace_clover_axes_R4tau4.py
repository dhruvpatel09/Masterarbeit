#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

files = [
    "results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100.csv",
    "results/laplace_clover_probe_R4tau4_axis1_t0avg48_n1-100.csv",
    "results/laplace_clover_probe_R4tau4_axis2_t0avg48_n1-100.csv",
]

out = Path("results/laplace_clover_probe_R4tau4_axes0-2_t0avg48_n1-100.csv")

dfs = [pd.read_csv(f) for f in files]
df = pd.concat(dfs, ignore_index=True)
df = df.sort_values(["cfg", "axis", "dx", "dy"]).reset_index(drop=True)

assert df.shape[0] == 3 * 100 * 117, df.shape
assert sorted(df["cfg"].unique()) == list(range(1, 101))
assert sorted(df["axis"].unique()) == [0, 1, 2]
assert sorted(df["dx"].unique()) == list(range(-2, 7))
assert sorted(df["dy"].unique()) == list(range(-6, 7))
assert (df["Nsrc"] == 48 * 24**3).all()
assert df.notna().all().all()

# Check every cfg/probe has all 3 axes.
check = df.groupby(["cfg", "r", "tau", "dx", "dy", "dz"])["axis"].nunique()
assert check.min() == 3 and check.max() == 3, check.describe()

df.to_csv(out, index=False, lineterminator="\n")

print("PASS wrote:", out)
print("shape =", df.shape)
print()
print(df.groupby("axis")[["rho_E2","rho_B2","rho_S"]].mean())
