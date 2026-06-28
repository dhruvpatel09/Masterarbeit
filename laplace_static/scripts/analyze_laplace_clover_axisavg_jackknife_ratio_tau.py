#!/usr/bin/env python3

from pathlib import Path
import sys
import numpy as np
import pandas as pd

if len(sys.argv) != 2:
    raise SystemExit("Usage: analyze_laplace_clover_axisavg_jackknife_ratio_tau.py TAU")

tau = int(sys.argv[1])

OBS = {
    "E2": ("LE2_Re", "E2_vac"),
    "B2": ("LB2_Re", "B2_vac"),
    "S":  ("LS_Re",  "S_vac"),
}

inp = Path(f"results/laplace_clover_probe_R4tau{tau}_axes0-2_t0avg48_n1-100.csv")
out = Path(f"results/laplace_clover_probe_R4tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")

df = pd.read_csv(inp)

axisavg = (
    df.groupby(["cfg", "r", "tau", "dx", "dy", "dz"], as_index=False)
      .agg(
          axis_count=("axis", "nunique"),
          L_Re=("L_Re", "mean"),
          LE2_Re=("LE2_Re", "mean"),
          LB2_Re=("LB2_Re", "mean"),
          LS_Re=("LS_Re", "mean"),
          E2_vac=("E2_vac", "mean"),
          B2_vac=("B2_vac", "mean"),
          S_vac=("S_vac", "mean"),
      )
)

assert (axisavg["axis_count"] == 3).all()

rows = []
group_cols = ["r", "tau", "dx", "dy", "dz"]

for key, g in axisavg.groupby(group_cols, sort=True):
    g = g.sort_values("cfg")
    n = len(g)
    assert n == 100, (key, n)

    L = g["L_Re"].to_numpy()

    row = dict(zip(group_cols, key))
    row["ncfg"] = n
    row["naxis"] = 3
    row["L_Re_mean"] = L.mean()

    jk_L = np.array([np.delete(L, i).mean() for i in range(n)])
    row["L_Re_err"] = np.sqrt((n - 1) / n * np.sum((jk_L - jk_L.mean())**2))

    for obs, (LO_col, O_col) in OBS.items():
        LO = g[LO_col].to_numpy()
        O = g[O_col].to_numpy()

        jk = []
        for i in range(n):
            mask = np.ones(n, dtype=bool)
            mask[i] = False
            rho = LO[mask].mean() / L[mask].mean() - O[mask].mean()
            jk.append(rho)

        jk = np.array(jk)
        row[f"rho_{obs}_mean"] = jk.mean()
        row[f"rho_{obs}_err"] = np.sqrt((n - 1) / n * np.sum((jk - jk.mean())**2))

    rows.append(row)

pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")

print("PASS wrote:", out)
print("rows =", len(rows))
