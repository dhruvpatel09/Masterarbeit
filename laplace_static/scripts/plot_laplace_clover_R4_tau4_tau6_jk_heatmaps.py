#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = Path("results/figures")
outdir.mkdir(parents=True, exist_ok=True)

inputs = [
    (
        4,
        "results/laplace_clover_probe_R4tau4_plane_t0avg48_n1-100_jk_ratio_summary.csv",
    ),
    (
        6,
        "results/laplace_clover_probe_R4tau6_plane_t0avg48_n1-100_jk_ratio_summary.csv",
    ),
]

# Use a common signed color scale for fair tau comparison.
all_vals = []
for tau, path in inputs:
    df = pd.read_csv(path)
    all_vals.append(df["rho_S_mean"].abs().max())

vmax = max(all_vals)
vmin = -vmax

for tau, path in inputs:
    df = pd.read_csv(path)
    grid = df.pivot(index="dy", columns="dx", values="rho_S_mean").sort_index(ascending=True)

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
    ax.set_title(rf"Jackknife ratio $\rho_S(\Delta x,\Delta y)$, $R/a=4$, $\tau/a={tau}$, $N_{{\rm cfg}}=100$")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\rho_S = \langle LS\rangle/\langle L\rangle - \langle S\rangle$")

    fig.tight_layout()
    out = outdir / f"laplace_clover_R4tau{tau}_rho_S_heatmap_t0avg48_n100_jk_ratio_publication.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Wrote", out)
