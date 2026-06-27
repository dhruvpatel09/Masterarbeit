#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

s4 = pd.read_csv("results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100_summary.csv")
s6 = pd.read_csv("results/laplace_clover_probe_R4tau6_plane_t0avg48_n1-100_summary.csv")

l4 = s4[s4["dy"] == 0].sort_values("dx")
l6 = s6[s6["dy"] == 0].sort_values("dx")

fig, ax = plt.subplots(figsize=(7.2, 4.8))

ax.errorbar(
    l4["dx"], l4["rho_S_mean"], yerr=l4["rho_S_err"],
    marker="o", linestyle="-", capsize=3, label=r"$\tau/a=4$"
)

ax.errorbar(
    l6["dx"], l6["rho_S_mean"], yerr=l6["rho_S_err"],
    marker="s", linestyle="-", capsize=3, label=r"$\tau/a=6$"
)

ax.axhline(0, linewidth=1.0, color="black", alpha=0.6)
ax.axvline(0, linestyle="--", linewidth=1.0, color="black", alpha=0.75)
ax.axvline(4, linestyle="--", linewidth=1.0, color="black", alpha=0.75)

ax.set_xlabel(r"$\Delta x/a$ at $\Delta y=0$")
ax.set_ylabel(r"$\rho_S$")
ax.set_title(r"Tau stability of axial $\rho_S$, $R/a=4$, $N_{\rm cfg}=100$")
ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
out = outdir / "laplace_clover_R4_rho_S_axis_profile_tau4_tau6_t0avg48_n100_compare.png"
fig.savefig(out, dpi=300)
plt.close(fig)
print("Wrote", out)

comp = pd.DataFrame({
    "dx": l4["dx"].to_numpy(),
    "rho_S_tau4": l4["rho_S_mean"].to_numpy(),
    "err_tau4": l4["rho_S_err"].to_numpy(),
    "rho_S_tau6": l6["rho_S_mean"].to_numpy(),
    "err_tau6": l6["rho_S_err"].to_numpy(),
})
comp["delta_tau6_minus_tau4"] = comp["rho_S_tau6"] - comp["rho_S_tau4"]
comp.to_csv("results/laplace_clover_R4_rho_S_axis_tau4_tau6_compare_t0avg48_n100.csv", index=False)

print()
print(comp.to_string(index=False))
