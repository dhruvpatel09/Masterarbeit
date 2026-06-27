#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

OBS = {
    "E2": ("LE2_Re", "E2_vac"),
    "B2": ("LB2_Re", "B2_vac"),
    "S":  ("LS_Re",  "S_vac"),
}

def jackknife_ratio_for_file(inp, out):
    df = pd.read_csv(inp)

    rows = []
    group_cols = ["r", "tau", "axis", "dx", "dy", "dz"]

    for key, g in df.groupby(group_cols, sort=True):
        cfgs = np.array(sorted(g["cfg"].unique()))
        n = len(cfgs)
        assert n == 100, (inp, key, n)

        L = g.sort_values("cfg")["L_Re"].to_numpy()

        row = dict(zip(group_cols, key))
        row["ncfg"] = n
        row["L_Re_mean"] = L.mean()
        row["L_Re_err"] = np.sqrt((n - 1) / n * np.sum((np.array([
            np.delete(L, i).mean() for i in range(n)
        ]) - L.mean())**2))

        for obs, (LO_col, O_col) in OBS.items():
            LO = g.sort_values("cfg")[LO_col].to_numpy()
            O  = g.sort_values("cfg")[O_col].to_numpy()

            jk = []
            for i in range(n):
                mask = np.ones(n, dtype=bool)
                mask[i] = False
                rho = LO[mask].mean() / L[mask].mean() - O[mask].mean()
                jk.append(rho)

            jk = np.array(jk)
            mean = jk.mean()
            err = np.sqrt((n - 1) / n * np.sum((jk - mean)**2))

            row[f"rho_{obs}_mean"] = mean
            row[f"rho_{obs}_err"] = err

        rows.append(row)

    out = Path(out)
    pd.DataFrame(rows).to_csv(out, index=False, lineterminator="\n")
    print("Wrote", out)


jackknife_ratio_for_file(
    "results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100.csv",
    "results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100_jk_ratio_summary.csv",
)

jackknife_ratio_for_file(
    "results/laplace_clover_probe_R4tau6_plane_t0avg48_n1-100.csv",
    "results/laplace_clover_probe_R4tau6_plane_t0avg48_n1-100_jk_ratio_summary.csv",
)
