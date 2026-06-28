#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

s4 = pd.read_csv("results/laplace_clover_probe_R4tau4_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")
s6 = pd.read_csv("results/laplace_clover_probe_R4tau6_axes0-2_t0avg48_n1-100_jk_ratio_summary.csv")

l4 = s4[s4["dy"] == 0].sort_values("dx")
l6 = s6[s6["dy"] == 0].sort_values("dx")

comp = pd.DataFrame({
    "dx": l4["dx"].to_numpy(),
    "rho_S_tau4": l4["rho_S_mean"].to_numpy(),
    "err_tau4": l4["rho_S_err"].to_numpy(),
    "rho_S_tau6": l6["rho_S_mean"].to_numpy(),
    "err_tau6": l6["rho_S_err"].to_numpy(),
})

comp["delta_tau6_minus_tau4"] = comp["rho_S_tau6"] - comp["rho_S_tau4"]
comp["pull"] = comp["delta_tau6_minus_tau4"] / np.sqrt(comp["err_tau4"]**2 + comp["err_tau6"]**2)

fig, ax = plt.subplots(figsize=(7.2, 4.8))

ax.errorbar(
    comp["dx"], comp["rho_S_tau4"], yerr=comp["err_tau4"],
    marker="o", linestyle="-", capsize=3, label=r"$\tau/a=4$"
)

ax.errorbar(
    comp["dx"], comp["rho_S_tau6"], yerr=comp["err_tau6"],
    marker="s", linestyle="-", capsize=3, label=r"$\tau/a=6$"
)

ax.axhline(0, linewidth=1.0, color="black", alpha=0.6)
ax.axvline(0, linestyle="--", linewidth=1.0, color="black", alpha=0.75)
ax.axvline(4, linestyle="--", linewidth=1.0, color="black", alpha=0.75)

ax.set_xlabel(r"$\Delta x/a$ at $\Delta y=0$")
ax.set_ylabel(r"$\rho_S$")
ax.set_title(r"Axis-averaged tau stability of $\rho_S$, $R/a=4$, $N_{\rm cfg}=100$")
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
out = outdir / "laplace_clover_R4_rho_S_axis_profile_tau4_tau6_t0avg48_n100_axisavg_jk_ratio_compare.png"
fig.savefig(out, dpi=300)
plt.close(fig)

comp_out = "results/laplace_clover_R4_rho_S_axis_tau4_tau6_axisavg_jk_ratio_compare_t0avg48_n100.csv"
comp.to_csv(comp_out, index=False)

print("Wrote", out)
print("Wrote", comp_out)
print()
print(comp.to_string(index=False))
print()
print("max |pull| =", comp["pull"].abs().max())
