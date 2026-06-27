#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

inp = Path("results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100_summary.csv")
outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inp)

# ----------------------------------------------------------------------
# 1. Polished rho_S heatmap
# ----------------------------------------------------------------------
obs = "rho_S_mean"
err = "rho_S_err"

grid = df.pivot(index="dy", columns="dx", values=obs).sort_index(ascending=True)
vmax = abs(grid.values).max()
vmin = -vmax

fig, ax = plt.subplots(figsize=(7.2, 5.2))

im = ax.imshow(
    grid.values,
    origin="lower",
    aspect="auto",
    vmin=vmin,
    vmax=vmax,
    cmap="RdBu_r",
    extent=[
        grid.columns.min() - 0.5,
        grid.columns.max() + 0.5,
        grid.index.min() - 0.5,
        grid.index.max() + 0.5,
    ],
)

ax.axvline(0, linestyle="--", linewidth=1.1, color="black", alpha=0.75)
ax.axvline(4, linestyle="--", linewidth=1.1, color="black", alpha=0.75)
ax.axhline(0, linestyle=":", linewidth=1.0, color="black", alpha=0.65)

ax.set_xlabel(r"$\Delta x/a$")
ax.set_ylabel(r"$\Delta y/a$")
ax.set_title(r"$\rho_S(\Delta x,\Delta y)$, $R/a=4$, $\tau/a=4$, $N_{\rm cfg}=100$")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$\rho_S = \langle LS\rangle/\langle L\rangle - \langle S\rangle$")

fig.tight_layout()
out = outdir / "laplace_clover_R4tau4_rho_S_heatmap_t0avg48_n100_publication.png"
fig.savefig(out, dpi=300)
plt.close(fig)
print("Wrote", out)

# ----------------------------------------------------------------------
# 2. rho_S line profile along the static-source axis dy=0
# ----------------------------------------------------------------------
line = df[df["dy"] == 0].sort_values("dx")

fig, ax = plt.subplots(figsize=(7.0, 4.6))

ax.errorbar(
    line["dx"],
    line[obs],
    yerr=line[err],
    marker="o",
    linestyle="-",
    capsize=3,
)

ax.axhline(0, linewidth=1.0, color="black", alpha=0.6)
ax.axvline(0, linestyle="--", linewidth=1.0, color="black", alpha=0.75)
ax.axvline(4, linestyle="--", linewidth=1.0, color="black", alpha=0.75)

ax.set_xlabel(r"$\Delta x/a$ at $\Delta y=0$")
ax.set_ylabel(r"$\rho_S$")
ax.set_title(r"Axial profile of $\rho_S$, $R/a=4$, $\tau/a=4$, $N_{\rm cfg}=100$")
ax.grid(True, alpha=0.3)

fig.tight_layout()
out = outdir / "laplace_clover_R4tau4_rho_S_axis_profile_t0avg48_n100_publication.png"
fig.savefig(out, dpi=300)
plt.close(fig)
print("Wrote", out)

# ----------------------------------------------------------------------
# 3. Text table for advisor notes
# ----------------------------------------------------------------------
tab = line[["dx", "rho_S_mean", "rho_S_err"]].copy()
out = Path("results/laplace_clover_R4tau4_rho_S_axis_profile_t0avg48_n100.txt")
tab.to_csv(out, sep=" ", index=False)
print("Wrote", out)

print()
print(tab.to_string(index=False))
