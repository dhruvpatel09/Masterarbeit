#!/usr/bin/env python3

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pyerrors as pe

if len(sys.argv) != 2:
    raise SystemExit("Usage: analyze_laplace_clover_axisavg_pyerrors_axial_R.py R")

R = int(sys.argv[1])

def make_axisavg_cfg_table(R, tau):
    inp = Path(f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100.csv")
    df = pd.read_csv(inp)

    axisavg = (
        df.groupby(["cfg", "r", "tau", "dx", "dy", "dz"], as_index=False)
          .agg(
              axis_count=("axis", "nunique"),
              L_Re=("L_Re", "mean"),
              LS_Re=("LS_Re", "mean"),
              S_vac=("S_vac", "mean"),
          )
    )

    assert (axisavg["axis_count"] == 3).all()
    return axisavg

rows = []

for tau in [4, 6]:
    df = make_axisavg_cfg_table(R, tau)

    for dx in sorted(df["dx"].unique()):
        g = df[(df["dy"] == 0) & (df["dx"] == dx)].sort_values("cfg")
        assert len(g) == 100, (R, tau, dx, len(g))

        # Important: same ensemble label for L, LS, S so pyerrors keeps correlations.
        ens = f"Em1_R{R}_tau{tau}_dx{dx}_dy0_axisavg"

        L  = pe.Obs([g["L_Re"].to_numpy()],  [ens])
        LS = pe.Obs([g["LS_Re"].to_numpy()], [ens])
        S  = pe.Obs([g["S_vac"].to_numpy()], [ens])

        rho = LS / L - S
        rho.gamma_method()

        rows.append({
            "R": R,
            "tau": tau,
            "dx": dx,
            "dy": 0,
            "rho_S_pyerrors": rho.value,
            "err_pyerrors": rho.dvalue,
        })

out = Path(f"results/laplace_clover_R{R}_axisavg_rho_S_axial_tau4_tau6_pyerrors_corr.csv")
pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")

print("Wrote", out)
print()
print(pd.DataFrame(rows).to_string(index=False))
