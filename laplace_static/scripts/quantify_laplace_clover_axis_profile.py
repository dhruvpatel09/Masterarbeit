#!/usr/bin/env python3

from pathlib import Path
import math
import pandas as pd

def load_line(R, tau):
    df = pd.read_csv(
        f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv"
    )
    return df[df["dy"] == 0].sort_values("dx").copy()

def val_at(line, dx):
    row = line[line["dx"] == dx]
    if len(row) != 1:
        return None, None
    r = row.iloc[0]
    return float(r["rho_S_mean"]), float(r["rho_S_err"])

def midpoint_estimate(line, R):
    if R % 2 == 0:
        return val_at(line, R // 2)

    dx1 = R // 2
    dx2 = dx1 + 1
    v1, e1 = val_at(line, dx1)
    v2, e2 = val_at(line, dx2)

    # Diagnostic error only, ignoring covariance.
    return 0.5 * (v1 + v2), 0.5 * math.sqrt(e1**2 + e2**2)

rows = []

for R in [4, 5]:
    for tau in [4, 6]:
        path = Path(f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")
        if not path.exists():
            continue

        line = load_line(R, tau)
        inside = line[(line["dx"] >= 0) & (line["dx"] <= R)].copy()

        left = inside[inside["dx"] <= R / 2.0]
        right = inside[inside["dx"] >= R / 2.0]

        left_min = left.loc[left["rho_S_mean"].idxmin()]
        right_min = right.loc[right["rho_S_mean"].idxmin()]

        mid_v, mid_e = midpoint_estimate(line, R)
        source0_v, source0_e = val_at(line, 0)
        sourceR_v, sourceR_e = val_at(line, R)
        outerL_v, outerL_e = val_at(line, -2)
        outerR_v, outerR_e = val_at(line, R + 2)

        rows.append({
            "R": R,
            "tau": tau,
            "left_min_dx": int(left_min["dx"]),
            "left_min_rho": float(left_min["rho_S_mean"]),
            "left_min_err": float(left_min["rho_S_err"]),
            "right_min_dx": int(right_min["dx"]),
            "right_min_rho": float(right_min["rho_S_mean"]),
            "right_min_err": float(right_min["rho_S_err"]),
            "midpoint_rho": mid_v,
            "midpoint_err_diag": mid_e,
            "source0_rho": source0_v,
            "source0_err": source0_e,
            "sourceR_rho": sourceR_v,
            "sourceR_err": sourceR_e,
            "outer_left_rho": outerL_v,
            "outer_left_err": outerL_e,
            "outer_right_rho": outerR_v,
            "outer_right_err": outerR_e,
            "dip_asymmetry_right_minus_left": float(right_min["rho_S_mean"] - left_min["rho_S_mean"]),
        })

out = pd.DataFrame(rows)

path = "results/laplace_clover_axis_profile_shape_observables_R4_R5.csv"
out.to_csv(path, index=False)

print("Wrote", path)
print()
print(out.to_string(index=False))
