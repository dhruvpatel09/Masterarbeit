#!/usr/bin/env python3

from pathlib import Path
import math
import numpy as np
import pandas as pd

def load_line(R, tau):
    df = pd.read_csv(
        f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv"
    )
    return df[df["dy"] == 0].sort_values("dx").copy()

def val_at(line, dx):
    row = line[line["dx"] == dx]
    if len(row) != 1:
        return np.nan, np.nan
    r = row.iloc[0]
    return float(r["rho_S_mean"]), float(r["rho_S_err"])

def midpoint_estimate(line, R):
    if R % 2 == 0:
        return val_at(line, R // 2)
    dx1 = R // 2
    dx2 = dx1 + 1
    v1, e1 = val_at(line, dx1)
    v2, e2 = val_at(line, dx2)
    return 0.5 * (v1 + v2), 0.5 * math.sqrt(e1**2 + e2**2)

rows = []

for R in [4, 5, 6]:
    line4 = load_line(R, 4)

    inside = line4[(line4["dx"] >= 0) & (line4["dx"] <= R)].copy()
    left = inside[inside["dx"] <= R / 2.0]
    right = inside[inside["dx"] >= R / 2.0]

    left_min = left.loc[left["rho_S_mean"].idxmin()]
    right_min = right.loc[right["rho_S_mean"].idxmin()]
    mid_v, mid_e = midpoint_estimate(line4, R)

    tau_table = pd.read_csv(
        f"results/laplace_clover_R{R}_axisavg_tau4_tau6_pyerrors_advisor_table.txt",
        comment="#",
        delim_whitespace=True,
        names=["dx", "rhoS_tau4", "err4", "rhoS_tau6", "err6", "delta", "pull"],
    )

    max_pull_all = tau_table["pull"].abs().max()
    max_pull_inside = tau_table[(tau_table["dx"] >= 0) & (tau_table["dx"] <= R)]["pull"].abs().max()

    rows.append({
        "R": R,
        "tau_main": 4,
        "left_min_dx": int(left_min["dx"]),
        "left_min_xmid": float(left_min["dx"] - R/2),
        "left_min_rho": float(left_min["rho_S_mean"]),
        "left_min_err": float(left_min["rho_S_err"]),
        "right_min_dx": int(right_min["dx"]),
        "right_min_xmid": float(right_min["dx"] - R/2),
        "right_min_rho": float(right_min["rho_S_mean"]),
        "right_min_err": float(right_min["rho_S_err"]),
        "midpoint_rho": mid_v,
        "midpoint_err_diag": mid_e,
        "dip_asymmetry_right_minus_left": float(right_min["rho_S_mean"] - left_min["rho_S_mean"]),
        "max_pull_tau4_vs_tau6_all": max_pull_all,
        "max_pull_tau4_vs_tau6_inside": max_pull_inside,
    })

out = pd.DataFrame(rows)

csv_path = "results/laplace_clover_R456_axisavg_advisor_summary.csv"
txt_path = "results/laplace_clover_R456_axisavg_advisor_summary.txt"

out.to_csv(csv_path, index=False)

with open(txt_path, "w") as f:
    f.write("# Axis-averaged Laplace-clover rho_S summary, Ncfg=100\n")
    f.write("# Main profile uses tau/a=4; tau/a=6 used as stability check.\n")
    f.write("# rho_S = <LS>/<L> - <S>\n")
    f.write("# Errors in profile columns are jackknife; tau-stability pulls use pyerrors/Gamma errors.\n")
    f.write("\n")
    f.write(out.to_string(index=False))
    f.write("\n")

print("Wrote", csv_path)
print("Wrote", txt_path)
print()
print(out.to_string(index=False))
