#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

VEFF_PATH = RESULTS / "laplace_Veff_R0-12_T1-9_t0avg6_n1-100_pyerrors.csv"
STATUS_PATH = RESULTS / "laplace_extended_R0-12_T1-10_t0avg6_n1-100_status.csv"
RAW_EXT_PATH = RESULTS / "spatial_avg_RT_R0-12_T1-10_t0avg6_n1-100.csv"
RAW_R6_PATH = RESULTS / "spatial_avg_RT_R0-6_T1-8_t0avg6_n1-100.csv"
POT_R6_PATH = RESULTS / "laplace_static_potential_R1-6_T3-5_t0avg6_correlated.csv"
FIT_STABILITY_PATH = RESULTS / "laplace_static_potential_fit_stability_t0avg6_R6.csv"

OUT_SUMMARY = RESULTS / "laplace_extended_summary_t0avg6_R12T10.csv"
OUT_DIAG_POINTS = RESULTS / "laplace_static_potential_diagnostic_points_R7-12_t0avg6.csv"
OUT_OVERLAP = RESULTS / "laplace_R6_vs_R12T10_overlap_check_t0avg6.csv"

GOOD_MARGINAL = ["good", "marginal"]


def require(path):
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Wrote: {path}")


def load_inputs():
    veff = pd.read_csv(require(VEFF_PATH))
    status = pd.read_csv(require(STATUS_PATH))
    raw_ext = pd.read_csv(require(RAW_EXT_PATH))
    raw_r6 = pd.read_csv(require(RAW_R6_PATH))
    pot_r6 = pd.read_csv(require(POT_R6_PATH))
    fit = pd.read_csv(require(FIT_STABILITY_PATH))
    return veff, status, raw_ext, raw_r6, pot_r6, fit


def build_summary(veff, status):
    summary = veff.merge(
        status[
            [
                "R",
                "T",
                "mean_Re",
                "sem_Re",
                "relative_err_Re_percent",
                "nonpositive_configs",
            ]
        ],
        on=["R", "T"],
        how="left",
    )

    def region(row):
        R = int(row["R"])
        if R <= 0:
            return "R0_normalization"
        if R <= 6:
            return "validated_R1-6"
        if R <= 8:
            return "borderline_R7-8"
        return "diagnostic_R9-12"

    summary["region"] = summary.apply(region, axis=1)
    summary.to_csv(OUT_SUMMARY, index=False, lineterminator="\n")
    print(f"Wrote: {OUT_SUMMARY}")
    return summary


def choose_diagnostic_points(veff):
    rows = []

    for R in range(7, 13):
        group = veff[
            (veff["R"] == R)
            & (veff["status"] == "ok")
            & (veff["reliability"].isin(GOOD_MARGINAL))
        ].copy()

        if group.empty:
            rows.append(
                {
                    "R": R,
                    "chosen_T": np.nan,
                    "Veff": np.nan,
                    "err_pyerrors": np.nan,
                    "relative_err_percent": np.nan,
                    "reliability": "none",
                    "selection_rule": "no good/marginal point available",
                }
            )
            continue

        group = group.sort_values(
            ["T", "relative_err_percent"],
            ascending=[False, True],
        )
        chosen = group.iloc[0]

        rows.append(
            {
                "R": int(chosen["R"]),
                "chosen_T": int(chosen["T"]),
                "Veff": float(chosen["Veff"]),
                "err_pyerrors": float(chosen["err_pyerrors"]),
                "relative_err_percent": float(chosen["relative_err_percent"]),
                "reliability": str(chosen["reliability"]),
                "selection_rule": "latest T with reliability in {good,marginal}",
            }
        )

    diag = pd.DataFrame(rows)
    diag.to_csv(OUT_DIAG_POINTS, index=False, lineterminator="\n")
    print(f"Wrote: {OUT_DIAG_POINTS}")
    return diag


def plot_plateau_candidates(veff):
    plot_df = veff[
        (veff["R"].between(1, 12))
        & (veff["status"] == "ok")
        & (veff["reliability"].isin(GOOD_MARGINAL))
    ].copy()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.4, 5.2),
        sharey=True,
    )

    for ax, r_values, title in [
        (axes[0], range(1, 7), "validated region: R/a = 1..6"),
        (axes[1], range(7, 13), "extended diagnostic region: R/a = 7..12"),
    ]:
        for R in r_values:
            group = plot_df[plot_df["R"] == R].sort_values("T")

            if group.empty:
                continue

            ax.errorbar(
                group["T"],
                group["Veff"],
                yerr=group["err_pyerrors"],
                marker="o",
                capsize=3,
                linewidth=1.0,
                label=f"R/a={R}",
            )

        ax.set_title(title)
        ax.set_xlabel("T/a")
        ax.set_xticks(range(1, 10))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)

    axes[0].set_ylabel("a Veff(R,T)")
    fig.suptitle("Extended plateau-candidate points, pyerrors/Gamma method")

    savefig(FIGDIR / "laplace_Veff_plateau_candidates_t0avg6_R1-12.png")


def plot_fit_plus_diagnostic(pot_r6, fit, diag):
    fit_row = fit[
        (fit["model"] == "cornell")
        & (fit["fit_label"] == "R1-6")
        & (fit["status"] == "ok")
    ]

    if fit_row.empty:
        raise RuntimeError("Could not find Cornell R1-6 fit row")

    fit_row = fit_row.iloc[0]
    V0 = float(fit_row["V0"])
    sigma = float(fit_row["sigma"])
    alpha = float(fit_row["alpha"])

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

    valid_diag = diag.dropna(subset=["Veff", "err_pyerrors"])

    plt.errorbar(
        valid_diag["R"],
        valid_diag["Veff"],
        yerr=valid_diag["err_pyerrors"],
        fmt="s",
        capsize=3,
        label="diagnostic candidate points, R/a=7..12",
    )

    for _, row in valid_diag.iterrows():
        plt.annotate(
            f"T={int(row['chosen_T'])}",
            (row["R"], row["Veff"]),
            textcoords="offset points",
            xytext=(3, 6),
            fontsize=8,
        )

    plt.axvline(6, linestyle=":", linewidth=1.0, alpha=0.7)
    plt.axvline(12, linestyle=":", linewidth=1.0, alpha=0.7)

    plt.xlabel("R/a")
    plt.ylabel("aV(R)")
    plt.title("Validated R1-6 static potential plus R7-12 diagnostic extension")
    plt.xticks(range(1, 13))
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)

    savefig(
        FIGDIR
        / "laplace_static_potential_R1-6_fit_plus_R7-12_diagnostic_t0avg6.png"
    )


def plot_early_time_potential(veff):
    plt.figure(figsize=(8.4, 5.6))

    for T in [2, 3]:
        group = veff[
            (veff["T"] == T)
            & (veff["R"].between(1, 12))
            & (veff["status"] == "ok")
        ].sort_values("R")

        plt.errorbar(
            group["R"],
            group["Veff"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=3,
            linewidth=1.0,
            label=f"T/a={T}",
        )

    plt.axvline(
        12,
        linestyle="--",
        linewidth=1.0,
        alpha=0.7,
        label="half lattice",
    )

    plt.xlabel("R/a")
    plt.ylabel("aVeff(R,T)")
    plt.title("Early-time extended potential, pyerrors/Gamma method")
    plt.xticks(range(1, 13))
    plt.grid(True, alpha=0.3)
    plt.legend()

    savefig(
        FIGDIR
        / "laplace_static_potential_early_time_T2_T3_R1-12_t0avg6.png"
    )


def plot_reliability_map(veff):
    label_to_code = {
        "invalid": 0,
        "diagnostic_only": 1,
        "marginal": 2,
        "good": 3,
    }

    table = veff[veff["R"].between(1, 12)].copy()
    table["code"] = table["reliability"].map(label_to_code).fillna(0)

    pivot = (
        table.pivot(
            index="R",
            columns="T",
            values="code",
        )
        .sort_index(ascending=False)
    )

    cmap = ListedColormap(
        [
            "#3b3b3b",
            "#d95f02",
            "#e6ab02",
            "#1b9e77",
        ]
    )
    norm = BoundaryNorm(
        [-0.5, 0.5, 1.5, 2.5, 3.5],
        cmap.N,
    )

    plt.figure(figsize=(8.2, 5.8))
    image = plt.imshow(
        pivot.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )

    cbar = plt.colorbar(
        image,
        ticks=[0, 1, 2, 3],
    )
    cbar.ax.set_yticklabels(
        [
            "invalid",
            "diagnostic",
            "marginal",
            "good",
        ]
    )

    plt.xticks(
        range(len(pivot.columns)),
        pivot.columns,
    )
    plt.yticks(
        range(len(pivot.index)),
        pivot.index,
    )

    plt.xlabel("T/a")
    plt.ylabel("R/a")
    plt.title("Reliability map for extended Veff, pyerrors/Gamma method")

    savefig(FIGDIR / "laplace_reliability_map_t0avg6_R12T10.png")


def plot_relative_error_heatmap(veff):
    table = veff[veff["R"].between(1, 12)].copy()
    table.loc[
        table["status"] != "ok",
        "relative_err_percent",
    ] = np.nan

    pivot = (
        table.pivot(
            index="R",
            columns="T",
            values="relative_err_percent",
        )
        .sort_index(ascending=False)
    )

    values = pivot.to_numpy(dtype=float)
    log_values = np.log10(values)

    plt.figure(figsize=(8.2, 5.8))
    image = plt.imshow(
        log_values,
        aspect="auto",
        interpolation="nearest",
    )

    cbar = plt.colorbar(image)
    cbar.set_label("log10(relative uncertainty [%])")

    plt.xticks(
        range(len(pivot.columns)),
        pivot.columns,
    )
    plt.yticks(
        range(len(pivot.index)),
        pivot.index,
    )

    plt.xlabel("T/a")
    plt.ylabel("R/a")
    plt.title("Relative-error heatmap for extended Veff")

    savefig(FIGDIR / "laplace_relative_error_heatmap_t0avg6_R12T10.png")


def plot_nonpositive_fraction_heatmap(status):
    table = status[status["R"].between(1, 12)].copy()
    table["nonpositive_fraction"] = (
        table["nonpositive_configs"] / table["Ncfg"]
    )

    pivot = (
        table.pivot(
            index="R",
            columns="T",
            values="nonpositive_fraction",
        )
        .sort_index(ascending=False)
    )

    plt.figure(figsize=(8.2, 5.8))
    image = plt.imshow(
        pivot.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        vmin=0.0,
        vmax=1.0,
    )

    cbar = plt.colorbar(image)
    cbar.set_label("nonpositive fraction")

    plt.xticks(
        range(len(pivot.columns)),
        pivot.columns,
    )
    plt.yticks(
        range(len(pivot.index)),
        pivot.index,
    )

    plt.xlabel("T/a")
    plt.ylabel("R/a")
    plt.title("Nonpositive correlator fraction in extended scan")

    savefig(FIGDIR / "laplace_nonpositive_fraction_heatmap_t0avg6_R12T10.png")


def plot_overlap_check(raw_r6, raw_ext):
    old = raw_r6.sort_values(["cfg", "T", "R"]).copy()

    new = raw_ext[
        (raw_ext["R"].between(0, 6))
        & (raw_ext["T"].between(1, 8))
    ].copy()
    new = new.sort_values(["cfg", "T", "R"])

    merged = old.merge(
        new,
        on=["cfg", "T", "R", "Nsrc"],
        how="inner",
        suffixes=("_old", "_new"),
    )

    if len(merged) != len(old):
        raise RuntimeError("Overlap merge did not reproduce all old rows")

    merged["abs_delta_Re"] = (
        merged["Re_new"] - merged["Re_old"]
    ).abs()
    merged["abs_delta_Im"] = (
        merged["Im_new"] - merged["Im_old"]
    ).abs()

    summary = merged.groupby(
        ["R", "T"],
        as_index=False,
    ).agg(
        max_abs_delta_Re=("abs_delta_Re", "max"),
        max_abs_delta_Im=("abs_delta_Im", "max"),
    )

    summary.to_csv(
        OUT_OVERLAP,
        index=False,
        lineterminator="\n",
    )
    print(f"Wrote: {OUT_OVERLAP}")

    pivot = (
        summary.pivot(
            index="R",
            columns="T",
            values="max_abs_delta_Re",
        )
        .sort_index(ascending=False)
    )

    values = np.log10(
        pivot.to_numpy(dtype=float) + 1.0e-300
    )

    plt.figure(figsize=(7.4, 4.8))
    image = plt.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
    )

    cbar = plt.colorbar(image)
    cbar.set_label("log10(max |Delta Re|)")

    plt.xticks(
        range(len(pivot.columns)),
        pivot.columns,
    )
    plt.yticks(
        range(len(pivot.index)),
        pivot.index,
    )

    plt.xlabel("T/a")
    plt.ylabel("R/a")
    plt.title("Overlap check: old R6T8 data vs extended R12T10 data")

    savefig(FIGDIR / "laplace_R6_vs_R12T10_overlap_check_t0avg6.png")


def plot_plateau_window_R7_R8(veff):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.8, 4.9),
        sharey=True,
    )

    marker_for = {
        "good": "o",
        "marginal": "s",
        "diagnostic_only": "^",
    }

    for ax, R in zip(axes, [7, 8]):
        group = veff[
            (veff["R"] == R)
            & (veff["status"] == "ok")
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
        ax.set_xticks(range(1, 10))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    axes[0].set_ylabel("aVeff(R,T)")
    fig.suptitle("Plateau-window comparison for borderline radii R/a=7,8")

    savefig(FIGDIR / "laplace_plateau_window_comparison_R7_R8_t0avg6.png")


def main():
    veff, status, raw_ext, raw_r6, pot_r6, fit = load_inputs()

    summary = build_summary(veff, status)
    diag = choose_diagnostic_points(veff)

    plot_reliability_map(veff)
    plot_relative_error_heatmap(veff)
    plot_nonpositive_fraction_heatmap(status)
    plot_overlap_check(raw_r6, raw_ext)
    plot_plateau_candidates(veff)
    plot_early_time_potential(veff)
    plot_plateau_window_R7_R8(veff)
    plot_fit_plus_diagnostic(pot_r6, fit, diag)

    print()
    print("Extended diagnostic suite complete")
    print()
    print("Diagnostic candidate points R/a=7..12:")
    print(
        diag.to_string(
            index=False,
            float_format=lambda value: f"{value:.6g}",
        )
    )
    print()
    print("Generated summary rows:", len(summary))
    print("Generated diagnostic rows:", len(diag))


if __name__ == "__main__":
    main()
