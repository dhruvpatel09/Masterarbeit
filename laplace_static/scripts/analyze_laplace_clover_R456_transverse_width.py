#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

rows = []
profile_rows = []

for R in [4, 5, 6]:
    tau = 4
    path = Path(f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")
    if not path.exists():
        print("Missing", path)
        continue

    df = pd.read_csv(path)
    df = df.copy()
    df["x_mid"] = df["dx"] - R / 2.0
    df["weight"] = np.maximum(-df["rho_S_mean"], 0.0)

    for dx, g in df.groupby("dx"):
        g = g.sort_values("dy")
        denom = g["weight"].sum()

        if denom <= 0:
            width = np.nan
            integrated_depletion = 0.0
        else:
            width = np.sqrt(((g["dy"] ** 2) * g["weight"]).sum() / denom)
            integrated_depletion = denom

        rows.append({
            "R": R,
            "tau": tau,
            "dx": int(dx),
            "x_mid": float(dx - R / 2.0),
            "width_rms_dy": float(width),
            "integrated_depletion_minus_rho": float(integrated_depletion),
            "rho_axis_dy0": float(g[g["dy"] == 0]["rho_S_mean"].iloc[0]) if (g["dy"] == 0).any() else np.nan,
        })

    # Save selected transverse profiles: left dip, midpoint, right dip.
    selected_dx = [1, R // 2, R - 1]
    if R % 2 == 1:
        selected_dx = [1, R // 2, R // 2 + 1, R - 1]

    for dx in selected_dx:
        gg = df[df["dx"] == dx].sort_values("dy").copy()
        if len(gg) == 0:
            continue
        label = "selected"
        if dx == 1:
            label = "left_dip"
        elif dx == R - 1:
            label = "right_dip"
        elif abs(dx - R/2.0) <= 0.5:
            label = "midpoint_region"

        for _, r in gg.iterrows():
            profile_rows.append({
                "R": R,
                "tau": tau,
                "dx": int(r["dx"]),
                "x_mid": float(r["x_mid"]),
                "dy": int(r["dy"]),
                "profile_label": label,
                "rho_S_mean": float(r["rho_S_mean"]),
                "rho_S_err": float(r["rho_S_err"]),
                "minus_rho_S": float(max(-r["rho_S_mean"], 0.0)),
            })

width = pd.DataFrame(rows).sort_values(["R", "dx"])
profiles = pd.DataFrame(profile_rows).sort_values(["R", "dx", "dy"])

width_path = "results/laplace_clover_R456_tau4_transverse_width_summary.csv"
profile_path = "results/laplace_clover_R456_tau4_selected_transverse_profiles.csv"

width.to_csv(width_path, index=False)
profiles.to_csv(profile_path, index=False)

print("Wrote", width_path)
print("Wrote", profile_path)
print()
print(width.to_string(index=False))

# Plot width vs midpoint coordinate.
fig, ax = plt.subplots(figsize=(7.4, 4.8))

for R in [4, 5, 6]:
    g = width[width["R"] == R].sort_values("x_mid")
    ax.plot(g["x_mid"], g["width_rms_dy"], marker="o", linestyle="-", label=rf"$R/a={R}$")

ax.axvline(0, color="black", linestyle=":", linewidth=1.0, alpha=0.7)
ax.set_xlabel(r"$(\Delta x - R/2)/a$")
ax.set_ylabel(r"diagnostic RMS transverse width in $\Delta y/a$")
ax.set_title(r"Transverse width of $-\rho_S$: $R/a=4,5,6$, $\tau/a=4$")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

fig_path = outdir / "laplace_clover_R456_tau4_transverse_width_rms_dy.png"
fig.savefig(fig_path, dpi=300)
plt.close(fig)
print("Wrote", fig_path)

# Plot integrated depletion vs midpoint coordinate.
fig, ax = plt.subplots(figsize=(7.4, 4.8))

for R in [4, 5, 6]:
    g = width[width["R"] == R].sort_values("x_mid")
    ax.plot(g["x_mid"], g["integrated_depletion_minus_rho"], marker="o", linestyle="-", label=rf"$R/a={R}$")

ax.axvline(0, color="black", linestyle=":", linewidth=1.0, alpha=0.7)
ax.set_xlabel(r"$(\Delta x - R/2)/a$")
ax.set_ylabel(r"$\sum_{\Delta y}[-\rho_S]_+$")
ax.set_title(r"Integrated transverse depletion: $R/a=4,5,6$, $\tau/a=4$")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()

fig_path = outdir / "laplace_clover_R456_tau4_integrated_transverse_depletion.png"
fig.savefig(fig_path, dpi=300)
plt.close(fig)
print("Wrote", fig_path)
