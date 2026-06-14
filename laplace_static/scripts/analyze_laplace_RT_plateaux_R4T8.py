#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_distribution


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

RAW_FILE = (
    RESULTS /
    "spatial_avg_RT_R0-4_T1-8_n1-100.csv"
)

PYERRORS_FILE = (
    RESULTS /
    "laplace_Veff_R0-4_T1-7_n1-100_pyerrors.csv"
)

WINDOWS_OUT = (
    RESULTS /
    "laplace_Veff_R1-4_T1-7_correlated_plateau_windows.csv"
)

SELECTED_OUT = (
    RESULTS /
    "laplace_static_potential_R1-4_T3-5_correlated.csv"
)

PLATEAU_FIGURE = (
    FIGDIR /
    "laplace_Veff_correlated_plateau_T3-5_R1-4.png"
)

POTENTIAL_FIGURE = (
    FIGDIR /
    "laplace_static_potential_correlated_T3-5.png"
)

EXPECTED_CFGS = list(range(1, 101))
CORRELATOR_T_VALUES = list(range(1, 9))
VEFF_T_VALUES = list(range(1, 8))
R_VALUES = list(range(1, 5))

SELECTED_TMIN = 3
SELECTED_TMAX = 5

EIGENVALUE_CUTOFF = 1.0e-12


def jackknife_error(samples):
    samples = np.asarray(samples, dtype=float)
    n_samples = len(samples)

    if n_samples < 2:
        return np.nan

    center = samples.mean()

    return np.sqrt(
        (n_samples - 1) / n_samples
        * np.sum((samples - center) ** 2)
    )


def validate_raw(raw):
    required = {"cfg", "T", "R", "Nsrc", "Re", "Im"}
    missing = required - set(raw.columns)

    if missing:
        raise RuntimeError(
            f"Missing raw-data columns: {sorted(missing)}"
        )

    if len(raw) != 4000:
        raise RuntimeError(
            f"Expected 4000 raw rows, found {len(raw)}"
        )

    if raw.duplicated(["cfg", "T", "R"]).any():
        raise RuntimeError(
            "Duplicate (cfg,T,R) raw-data rows found"
        )

    cfgs = sorted(raw["cfg"].unique().tolist())
    t_values = sorted(raw["T"].unique().tolist())
    r_values = sorted(raw["R"].unique().tolist())

    if cfgs != EXPECTED_CFGS:
        raise RuntimeError(
            f"Unexpected configuration IDs: {cfgs}"
        )

    if t_values != CORRELATOR_T_VALUES:
        raise RuntimeError(
            f"Unexpected correlator T values: {t_values}"
        )

    if r_values != list(range(0, 5)):
        raise RuntimeError(
            f"Unexpected R values: {r_values}"
        )

    if not (raw["Nsrc"] == 24 ** 3).all():
        raise RuntimeError(
            "Unexpected Nsrc values found"
        )

    return cfgs


def covariance_inverse(covariance):
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    largest = eigenvalues.max()

    if not np.isfinite(largest) or largest <= 0.0:
        raise RuntimeError(
            "Covariance matrix has no positive eigenvalues"
        )

    cutoff = largest * EIGENVALUE_CUTOFF
    keep = eigenvalues > cutoff
    rank = int(np.count_nonzero(keep))

    if rank < 2:
        raise RuntimeError(
            f"Insufficient covariance rank: {rank}"
        )

    kept_values = eigenvalues[keep]
    kept_vectors = eigenvectors[:, keep]

    inverse = (
        kept_vectors
        @ np.diag(1.0 / kept_values)
        @ kept_vectors.T
    )

    condition_number = (
        kept_values.max() / kept_values.min()
    )

    return inverse, rank, condition_number


def build_veff_data(raw, R, cfgs):
    arrays = {}
    sums = {}

    for T in CORRELATOR_T_VALUES:
        values = (
            raw[
                (raw["R"] == R) &
                (raw["T"] == T)
            ]
            .sort_values("cfg")["Re"]
            .to_numpy(dtype=float)
        )

        if len(values) != len(cfgs):
            raise RuntimeError(
                f"Expected {len(cfgs)} samples at "
                f"R={R}, T={T}; found {len(values)}"
            )

        arrays[T] = values
        sums[T] = values.sum()

    central = []
    jackknife_samples = []

    for T in VEFF_T_VALUES:
        mean_T = arrays[T].mean()
        mean_Tp1 = arrays[T + 1].mean()

        if mean_T <= 0.0 or mean_Tp1 <= 0.0:
            raise RuntimeError(
                f"Non-positive central mean at R={R}, T={T}"
            )

        central.append(
            np.log(mean_T / mean_Tp1)
        )

    for omit in range(len(cfgs)):
        row = []

        for T in VEFF_T_VALUES:
            mean_T = (
                sums[T] - arrays[T][omit]
            ) / (len(cfgs) - 1)

            mean_Tp1 = (
                sums[T + 1] - arrays[T + 1][omit]
            ) / (len(cfgs) - 1)

            if mean_T <= 0.0 or mean_Tp1 <= 0.0:
                raise RuntimeError(
                    "Non-positive leave-one-out mean at "
                    f"R={R}, T={T}, omit={omit}"
                )

            row.append(
                np.log(mean_T / mean_Tp1)
            )

        jackknife_samples.append(row)

    central = np.asarray(central, dtype=float)
    jackknife_samples = np.asarray(
        jackknife_samples,
        dtype=float,
    )

    centered = (
        jackknife_samples
        - jackknife_samples.mean(axis=0)
    )

    covariance = (
        (len(cfgs) - 1) / len(cfgs)
        * centered.T
        @ centered
    )

    return central, jackknife_samples, covariance


def correlated_constant_fit(
    central,
    jackknife_samples,
    covariance,
    tmin,
    tmax,
):
    indices = np.arange(tmin - 1, tmax)

    values = central[indices]
    samples = jackknife_samples[:, indices]
    covariance_window = covariance[
        np.ix_(indices, indices)
    ]

    inverse, rank, condition_number = (
        covariance_inverse(covariance_window)
    )

    ones = np.ones(len(indices))

    denominator = ones @ inverse @ ones

    if not np.isfinite(denominator) or denominator <= 0.0:
        raise RuntimeError(
            f"Invalid GLS denominator for T={tmin}..{tmax}"
        )

    weights = inverse @ ones / denominator

    estimate = weights @ values

    fit_samples = samples @ weights
    fit_error = jackknife_error(fit_samples)

    analytic_error = np.sqrt(1.0 / denominator)

    residual = values - estimate

    chi2_value = residual @ inverse @ residual
    dof = rank - 1

    p_value = (
        chi2_distribution.sf(chi2_value, dof)
        if dof > 0
        else np.nan
    )

    return {
        "Tmin": tmin,
        "Tmax": tmax,
        "Npoints": len(indices),
        "plateau_value": estimate,
        "err_plateau_jackknife": fit_error,
        "err_plateau_GLS": analytic_error,
        "chi2": chi2_value,
        "dof": dof,
        "chi2_per_dof": (
            chi2_value / dof
            if dof > 0
            else np.nan
        ),
        "p_value": p_value,
        "covariance_rank": rank,
        "covariance_condition_number": condition_number,
        "max_weight_magnitude": np.max(
            np.abs(weights)
        ),
    }


def main():
    raw = pd.read_csv(RAW_FILE)
    raw = raw.sort_values(
        ["cfg", "T", "R"]
    ).reset_index(drop=True)

    cfgs = validate_raw(raw)

    pyerrors = pd.read_csv(PYERRORS_FILE)

    required_pyerrors = {
        "T",
        "R",
        "Veff",
        "err_pyerrors",
        "status",
    }

    missing = required_pyerrors - set(pyerrors.columns)

    if missing:
        raise RuntimeError(
            f"Missing pyerrors columns: {sorted(missing)}"
        )

    all_window_rows = []
    selected_rows = []

    central_by_R = {}
    jackknife_by_R = {}
    covariance_by_R = {}

    for R in R_VALUES:
        central, samples, covariance = (
            build_veff_data(raw, R, cfgs)
        )

        central_by_R[R] = central
        jackknife_by_R[R] = samples
        covariance_by_R[R] = covariance

        for tmin in range(1, 6):
            for tmax in range(tmin + 2, 8):
                fit = correlated_constant_fit(
                    central,
                    samples,
                    covariance,
                    tmin,
                    tmax,
                )

                fit["R"] = R
                fit["selected_common_window"] = (
                    tmin == SELECTED_TMIN and
                    tmax == SELECTED_TMAX
                )

                all_window_rows.append(fit)

        selected = correlated_constant_fit(
            central,
            samples,
            covariance,
            SELECTED_TMIN,
            SELECTED_TMAX,
        )

        selected["R"] = R
        selected["fit_window"] = (
            f"{SELECTED_TMIN}-{SELECTED_TMAX}"
        )

        selected_rows.append(selected)

    windows = pd.DataFrame(all_window_rows)
    selected = pd.DataFrame(selected_rows)

    ordered_window_columns = [
        "R",
        "Tmin",
        "Tmax",
        "Npoints",
        "plateau_value",
        "err_plateau_jackknife",
        "err_plateau_GLS",
        "chi2",
        "dof",
        "chi2_per_dof",
        "p_value",
        "covariance_rank",
        "covariance_condition_number",
        "max_weight_magnitude",
        "selected_common_window",
    ]

    windows = windows[
        ordered_window_columns
    ].sort_values(
        ["R", "Tmin", "Tmax"]
    )

    ordered_selected_columns = [
        "R",
        "fit_window",
        "Tmin",
        "Tmax",
        "Npoints",
        "plateau_value",
        "err_plateau_jackknife",
        "err_plateau_GLS",
        "chi2",
        "dof",
        "chi2_per_dof",
        "p_value",
        "covariance_rank",
        "covariance_condition_number",
        "max_weight_magnitude",
    ]

    selected = selected[
        ordered_selected_columns
    ].sort_values("R")

    windows.to_csv(
        WINDOWS_OUT,
        index=False,
        lineterminator="\n",
    )

    selected.to_csv(
        SELECTED_OUT,
        index=False,
        lineterminator="\n",
    )

    # ---------------------------------------------------------
    # Per-R effective-potential plateau panels
    # ---------------------------------------------------------
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(10.0, 7.5),
        sharex=True,
    )

    for R, axis in zip(R_VALUES, axes.flat):
        group = (
            pyerrors[
                (pyerrors["R"] == R) &
                (pyerrors["status"] == "ok")
            ]
            .sort_values("T")
        )

        fit = selected[
            selected["R"] == R
        ].iloc[0]

        value = fit["plateau_value"]
        error = fit["err_plateau_jackknife"]

        axis.errorbar(
            group["T"],
            group["Veff"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=4,
            label="Gamma method",
        )

        axis.axvspan(
            SELECTED_TMIN - 0.5,
            SELECTED_TMAX + 0.5,
            alpha=0.08,
            label="Selected fit window",
        )

        axis.axhspan(
            value - error,
            value + error,
            alpha=0.20,
            label="Plateau ± 1σ",
        )

        axis.axhline(
            value,
            linestyle="--",
            linewidth=1.3,
            label="Plateau central value",
        )

        axis.set_title(
            f"R/a = {R}: "
            f"{value:.5f} ± {error:.5f}, "
            f"p = {fit['p_value']:.3f}"
        )

        axis.set_ylabel("a Veff(R,T)")
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=7)

    for axis in axes[-1, :]:
        axis.set_xlabel("T/a")

    for axis in axes.flat:
        axis.set_xticks(VEFF_T_VALUES)

    figure.suptitle(
        "Laplace effective potential and correlated "
        "T/a = 3–5 plateaux"
    )

    figure.tight_layout()
    figure.savefig(
        PLATEAU_FIGURE,
        dpi=200,
    )
    plt.close(figure)

    # ---------------------------------------------------------
    # Extracted static potential
    # ---------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    plt.errorbar(
        selected["R"],
        selected["plateau_value"],
        yerr=selected["err_plateau_jackknife"],
        marker="o",
        capsize=5,
        linestyle="-",
    )

    plt.xlabel("R/a")
    plt.ylabel("a V(R)")
    plt.title(
        "Laplace static potential from correlated "
        "T/a = 3–5 plateaux"
    )
    plt.xticks(R_VALUES)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(
        POTENTIAL_FIGURE,
        dpi=200,
    )
    plt.close()

    print("Selected common plateau window:")
    print(
        selected[[
            "R",
            "fit_window",
            "plateau_value",
            "err_plateau_jackknife",
            "chi2",
            "dof",
            "chi2_per_dof",
            "p_value",
        ]].to_string(index=False)
    )

    print()
    print("Key stability windows:")

    key_windows = windows[
        windows.apply(
            lambda row: (
                (row["Tmin"], row["Tmax"])
                in {
                    (2, 4),
                    (3, 5),
                    (3, 6),
                    (4, 6),
                    (4, 7),
                    (5, 7),
                }
            ),
            axis=1,
        )
    ]

    print(
        key_windows[[
            "R",
            "Tmin",
            "Tmax",
            "plateau_value",
            "err_plateau_jackknife",
            "chi2_per_dof",
            "p_value",
        ]].to_string(index=False)
    )

    print()
    print("Wrote:")
    print(WINDOWS_OUT)
    print(SELECTED_OUT)
    print(PLATEAU_FIGURE)
    print(POTENTIAL_FIGURE)


if __name__ == "__main__":
    main()
