#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import pyerrors as pe

def make_axisavg_cfg_table(tau):
    inp = Path(f"results/laplace_clover_probe_R4tau{tau}_axes0-2_t0avg48_n1-100.csv")
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
    df = make_axisavg_cfg_table(tau)

    for dx in sorted(df["dx"].unique()):
        g = df[(df["dy"] == 0) & (df["dx"] == dx)].sort_values("cfg")
        assert len(g) == 100, (tau, dx, len(g))

        # One ensemble name is enough here; cfg order is the Markov-chain order.
        L  = pe.Obs([g["L_Re"].to_numpy()],  [f"R4tau{tau}_dx{dx}_L"])
        LS = pe.Obs([g["LS_Re"].to_numpy()], [f"R4tau{tau}_dx{dx}_LS"])
        S  = pe.Obs([g["S_vac"].to_numpy()], [f"R4tau{tau}_dx{dx}_S"])

        rho = LS / L - S
        rho.gamma_method()

        rows.append({
            "tau": tau,
            "dx": dx,
            "dy": 0,
            "rho_S_pyerrors": rho.value,
            "err_pyerrors": rho.dvalue,
        })

out = Path("results/laplace_clover_R4_axisavg_rho_S_axial_tau4_tau6_pyerrors.csv")
pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")

print("Wrote", out)
print()
print(pd.DataFrame(rows).to_string(index=False))
