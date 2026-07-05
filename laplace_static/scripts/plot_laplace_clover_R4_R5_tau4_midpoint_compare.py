#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

def load_line(R, tau):
    df = pd.read_csv(
        f"results/laplace_clover_probe_R{R}tau{tau}_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv"
    )
    line = df[df["dy"] == 0].sort_values("dx").copy()
    line["x_mid"] = line["dx"] - R / 2.0
    line["R"] = R
    line["tau"] = tau
    return line

l4 = load_line(4, 4)
l5 = load_line(5, 4)

fig, ax = plt.subplots(figsize=(7.4, 4.8))

ax.errorbar(
    l4["x_mid"], l4["rho_S_mean"], yerr=l4["rho_S_err"],
    marker="o", linestyle="-", capsize=3, label=r"$R/a=4,\ \tau/a=4$"
)

ax.errorbar(
    l5["x_mid"], l5["rho_S_mean"], yerr=l5["rho_S_err"],
    marker="s", linestyle="-", capsize=3, label=r"$R/a=5,\ \tau/a=4$"
)

ax.axhline(0, linewidth=1.0, color="black", alpha=0.6)
ax.axvline(0, linestyle=":", linewidth=1.0, color="black", alpha=0.75)

# Source markers relative to midpoint.
ax.axvline(-2.0, linestyle="--", linewidth=1.0, color="black", alpha=0.45)
ax.axvline(+2.0, linestyle="--", linewidth=1.0, color="black", alpha=0.45)
ax.axvline(-2.5, linestyle="--", linewidth=1.0, color="black", alpha=0.75)
ax.axvline(+2.5, linestyle="--", linewidth=1.0, color="black", alpha=0.75)

ax.set_xlabel(r"$(\Delta x - R/2)/a$ at $\Delta y=0$")
ax.set_ylabel(r"$\rho_S$")
ax.set_title(r"Midpoint-centered axial $\rho_S$: $R/a=4$ vs $R/a=5$, $\tau/a=4$")
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
out = outdir / "laplace_clover_R4_R5_tau4_rho_S_axis_profile_midpoint_axisavg_jk_ratio_compare.png"
fig.savefig(out, dpi=300)
plt.close(fig)

combined = pd.concat([l4, l5], ignore_index=True)
combined.to_csv(
    "results/laplace_clover_R4_R5_tau4_rho_S_axis_profile_midpoint_axisavg_jk_ratio_compare.csv",
    index=False,
)

print("Wrote", out)
print("Wrote results/laplace_clover_R4_R5_tau4_rho_S_axis_profile_midpoint_axisavg_jk_ratio_compare.csv")
print()
print(combined[["R", "tau", "dx", "x_mid", "rho_S_mean", "rho_S_err"]].to_string(index=False))
