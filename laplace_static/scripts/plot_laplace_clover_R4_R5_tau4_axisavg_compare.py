#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

s4 = pd.read_csv("results/laplace_clover_probe_R4tau4_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")
s5 = pd.read_csv("results/laplace_clover_probe_R5tau4_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")

l4 = s4[s4["dy"] == 0].sort_values("dx")
l5 = s5[s5["dy"] == 0].sort_values("dx")

fig, ax = plt.subplots(figsize=(7.2, 4.8))

ax.errorbar(
    l4["dx"], l4["rho_S_mean"], yerr=l4["rho_S_err"],
    marker="o", linestyle="-", capsize=3, label=r"$R/a=4,\ \tau/a=4$"
)

ax.errorbar(
    l5["dx"], l5["rho_S_mean"], yerr=l5["rho_S_err"],
    marker="s", linestyle="-", capsize=3, label=r"$R/a=5,\ \tau/a=4$"
)

ax.axhline(0, linewidth=1.0, color="black", alpha=0.6)
ax.axvline(0, linestyle="--", linewidth=1.0, color="black", alpha=0.75)
ax.axvline(4, linestyle="--", linewidth=1.0, color="black", alpha=0.45)
ax.axvline(5, linestyle="--", linewidth=1.0, color="black", alpha=0.75)

ax.set_xlabel(r"$\Delta x/a$ at $\Delta y=0$")
ax.set_ylabel(r"$\rho_S$")
ax.set_title(r"Axis-averaged axial $\rho_S$: $R/a=4$ vs $R/a=5$, $\tau/a=4$")
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
out = outdir / "laplace_clover_R4_R5_tau4_rho_S_axis_profile_axisavg_jk_ratio_compare.png"
fig.savefig(out, dpi=300)
plt.close(fig)

print("Wrote", out)
print()
print("R=4")
print(l4[["dx", "rho_S_mean", "rho_S_err"]].to_string(index=False))
print()
print("R=5")
print(l5[["dx", "rho_S_mean", "rho_S_err"]].to_string(index=False))
