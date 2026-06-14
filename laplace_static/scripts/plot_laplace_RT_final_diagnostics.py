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

PYERRORS_FILE = (
    RESULTS / "laplace_Veff_R0-3_T1-3_n1-100_pyerrors.csv"
)
COMPARISON_FILE = (
    RESULTS / "laplace_Veff_R0-3_T1-3_n1-100_error_comparison.csv"
)
RAW_FILE = (
    RESULTS / "spatial_avg_RT_R0-3_T1-4_n1-100.csv"
)


def jackknife_veff(raw, R):
    cfgs = sorted(raw["cfg"].unique())
    ncfg = len(cfgs)

    arrays = {}
    sums = {}

    for T in range(1, 5):
        values = (
            raw[
                (raw["R"] == R) &
                (raw["T"] == T)
            ]
            .sort_values("cfg")["Re"]
            .to_numpy(float)
        )

        arrays[T] = values
        sums[T] = values.sum()

    samples = []

    for omit in range(ncfg):
        row = []

        for T in range(1, 4):
            mean_T = (
                sums[T] - arrays[T][omit]
            ) / (ncfg - 1)

            mean_Tp1 = (
                sums[T + 1] - arrays[T + 1][omit]
            ) / (ncfg - 1)

            row.append(np.log(mean_T / mean_Tp1))

        samples.append(row)

    return np.asarray(samples)


def jackknife_covariance(samples):
    ncfg = len(samples)
    average = samples.mean(axis=0)
    centered = samples - average

    return (
        (ncfg - 1) / ncfg
        * centered.T
        @ centered
    )


def main():
    pyerr = pd.read_csv(PYERRORS_FILE)
    comparison = pd.read_csv(COMPARISON_FILE)
    raw = pd.read_csv(RAW_FILE)

    # ---------------------------------------------------------
    # 1. Zoomed pyerrors-only plateau plots
    # ---------------------------------------------------------
    for R in range(1, 4):
        g = pyerr[pyerr["R"] == R].sort_values("T")

        plt.figure(figsize=(7.0, 4.8))

        plt.errorbar(
            g["T"],
            g["Veff"],
            yerr=g["err_pyerrors"],
            marker="o",
            linewidth=1.5,
            capsize=5,
            elinewidth=1.4,
        )

        plt.xlabel("T/a")
        plt.ylabel("a Veff(R,T)")
        plt.title(f"Laplace effective potential, R/a = {R}")
        plt.xticks([1, 2, 3])
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        plt.savefig(
            FIGDIR / f"laplace_Veff_pyerrors_zoom_R{R}.png",
            dpi=200,
        )
        plt.close()

    # ---------------------------------------------------------
    # 2. Relative pyerrors uncertainty
    # ---------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        g = pyerr[pyerr["R"] == R].sort_values("T")

        relative = (
            100.0
            * g["err_pyerrors"]
            / np.abs(g["Veff"])
        )

        plt.plot(
            g["T"],
            relative,
            marker="o",
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("Relative uncertainty [%]")
    plt.title("Relative pyerrors uncertainty of Veff")
    plt.xticks([1, 2, 3])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        FIGDIR / "laplace_Veff_relative_error_pyerrors.png",
        dpi=200,
    )
    plt.close()

    # ---------------------------------------------------------
    # 3. Correlated drift relative to T/a = 3
    # ---------------------------------------------------------
    drift_rows = []

    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        samples = jackknife_veff(raw, R)
        covariance = jackknife_covariance(samples)

        central = np.array([
            comparison[
                (comparison["R"] == R) &
                (comparison["T"] == T)
            ]["Veff"].iloc[0]
            for T in [1, 2, 3]
        ])

        drift = central - central[2]

        drift_error = []

        for index in range(3):
            variance = (
                covariance[index, index]
                + covariance[2, 2]
                - 2.0 * covariance[index, 2]
            )

            drift_error.append(
                np.sqrt(max(variance, 0.0))
            )

        drift_error = np.asarray(drift_error)

        for index, T in enumerate([1, 2, 3]):
            significance = (
                drift[index] / drift_error[index]
                if drift_error[index] > 0.0
                else 0.0
            )

            drift_rows.append({
                "R": R,
                "T": T,
                "delta_to_T3": drift[index],
                "err_delta_jackknife": drift_error[index],
                "significance_sigma": significance,
            })

        plt.errorbar(
            [1, 2, 3],
            drift,
            yerr=drift_error,
            marker="o",
            capsize=4,
            label=f"R/a = {R}",
        )

    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("T/a")
    plt.ylabel("a[Veff(R,T) - Veff(R,3)]")
    plt.title("Correlated effective-potential drift")
    plt.xticks([1, 2, 3])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        FIGDIR / "laplace_Veff_drift_to_T3_jackknife.png",
        dpi=200,
    )
    plt.close()

    drift_df = pd.DataFrame(drift_rows)
    drift_df.to_csv(
        RESULTS /
        "laplace_Veff_R1-3_T1-3_drift_to_T3_jackknife.csv",
        index=False,
    )

    print("Drift diagnostics:")
    print(drift_df.to_string(index=False))

    print()
    print("Wrote:")
    for R in range(1, 4):
        print(
            FIGDIR /
            f"laplace_Veff_pyerrors_zoom_R{R}.png"
        )

    print(
        FIGDIR /
        "laplace_Veff_relative_error_pyerrors.png"
    )
    print(
        FIGDIR /
        "laplace_Veff_drift_to_T3_jackknife.png"
    )
    print(
        RESULTS /
        "laplace_Veff_R1-3_T1-3_drift_to_T3_jackknife.csv"
    )


if __name__ == "__main__":
    main()
