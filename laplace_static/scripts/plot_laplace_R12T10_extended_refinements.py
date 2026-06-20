#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

VEFF_PATH = RESULTS / "laplace_Veff_R0-12_T1-9_t0avg6_n1-100_pyerrors.csv"
POT_R6_PATH = RESULTS / "laplace_static_potential_R1-6_T3-5_t0avg6_correlated.csv"
FIT_STABILITY_PATH = RESULTS / "laplace_static_potential_fit_stability_t0avg6_R6.csv"

GOOD_MARGINAL = ["good", "marginal"]


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Wrote: {path}")


def load_inputs():
    veff = pd.read_csv(VEFF_PATH)
    pot_r6 = pd.read_csv(POT_R6_PATH)
    fit = pd.read_csv(FIT_STABILITY_PATH)
    return veff, pot_r6, fit


def get_cornell_fit(fit):
    row = fit[
        (fit["model"] == "cornell")
        & (fit["fit_label"] == "R1-6")
        & (fit["status"] == "ok")
    ]

    if row.empty:
        raise RuntimeError("Could not find Cornell R1-6 fit row")

    row = row.iloc[0]
    return float(row["V0"]), float(row["sigma"]), float(row["alpha"])


def write_common_t3_points(veff):
    common_t3 = veff[
        (veff["R"].between(7, 12))
        & (veff["T"] == 3)
        & (veff["status"] == "ok")
        & (veff["reliability"].isin(GOOD_MARGINAL))
    ].copy()

    common_t3 = common_t3[
        [
            "R",
            "T",
            "Veff",
            "err_pyerrors",
            "relative_err_percent",
            "reliability",
        ]
    ]

    common_t3["selection_rule"] = "common T=3 with reliability in {good,marginal}"

    out = RESULTS / "laplace_static_potential_diagnostic_points_R7-12_commonT3_t0avg6.csv"
    common_t3.to_csv(out, index=False, lineterminator="\n")
    print(f"Wrote: {out}")

    return common_t3


def plot_fit_plus_common_t3(veff, pot_r6, fit):
    V0, sigma, alpha = get_cornell_fit(fit)

    common_t3 = write_common_t3_points(veff)

    r_fit = np.linspace(1.0, 6.0, 300)
    y_fit = V0 + sigma * r_fit - alpha / r_fit

    r_extra = np.linspace(6.0, 12.0, 300)
    y_extra = V0 + sigma * r_extra - alpha / r_extra

    plt.figure(figsize=(8.4, 5.6))

    plt.plot(
        r_fit,
        y_fit,
        linewidth=1.8,
        label="Cornell fit, R/a=1..6",
    )

    plt.plot(
        r_extra,
        y_extra,
        linestyle="--",
        linewidth=1.4,
        alpha=0.8,
        label="Cornell extrapolation",
    )

    plt.errorbar(
        pot_r6["R"],
        pot_r6["plateau_value"],
        yerr=pot_r6["err_plateau_jackknife"],
        fmt="o",
        capsize=3,
        label="validated plateau points, R/a=1..6",
    )

    plt.errorbar(
        common_t3["R"],
        common_t3["Veff"],
        yerr=common_t3["err_pyerrors"],
        fmt="s",
        capsize=3,
        label="diagnostic common T/a=3 points, R/a=7..12",
    )

    for _, row in common_t3.iterrows():
        plt.annotate(
            row["reliability"],
            (row["R"], row["Veff"]),
            textcoords="offset points",
            xytext=(3, 6),
            fontsize=8,
        )

    plt.axvline(6, linestyle=":", linewidth=1.0, alpha=0.7)
    plt.axvline(12, linestyle=":", linewidth=1.0, alpha=0.7)

    plt.xlabel("R/a")
    plt.ylabel("aV(R)")
    plt.title("Validated R1-6 static potential plus common-T diagnostic extension")
    plt.xticks(range(1, 13))
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)

    savefig(
        FIGDIR
        / "laplace_static_potential_R1-6_fit_plus_R7-12_commonT3_diagnostic_t0avg6.png"
    )


def plot_R7_R8_good_marginal_zoom(veff):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.8, 4.9),
        sharey=True,
    )

    marker_for = {
        "good": "o",
        "marginal": "s",
    }

    for ax, R in zip(axes, [7, 8]):
        group = veff[
            (veff["R"] == R)
            & (veff["status"] == "ok")
            & (veff["reliability"].isin(GOOD_MARGINAL))
        ].sort_values("T")

        ax.axvspan(
            2.0,
            4.0,
            alpha=0.08,
            label="candidate T=2..4",
        )
        ax.axvspan(
            3.0,
            5.0,
            alpha=0.08,
            label="candidate T=3..5",
        )

        for reliability, marker in marker_for.items():
            sub = group[group["reliability"] == reliability]

            if sub.empty:
                continue

            ax.errorbar(
                sub["T"],
                sub["Veff"],
                yerr=sub["err_pyerrors"],
                fmt=marker,
                capsize=3,
                linestyle="none",
                label=reliability,
            )

        ax.plot(
            group["T"],
            group["Veff"],
            linewidth=1.0,
            alpha=0.5,
        )

        ax.set_title(f"R/a={R}")
        ax.set_xlabel("T/a")
        ax.set_xticks(range(1, 7))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("aVeff(R,T)")
    axes[0].set_ylim(0.9, 1.7)
    fig.suptitle("Good/marginal plateau-window zoom for R/a=7,8")

    savefig(
        FIGDIR
        / "laplace_plateau_window_comparison_R7_R8_good_marginal_zoom_t0avg6.png"
    )


def main():
    veff, pot_r6, fit = load_inputs()

    plot_fit_plus_common_t3(veff, pot_r6, fit)
    plot_R7_R8_good_marginal_zoom(veff)

    print()
    print("Extended refinement plots complete")


if __name__ == "__main__":
    main()
