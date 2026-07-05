#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_distribution


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"

RAW_INPUT = (
    RESULTS /
    "spatial_avg_RT_R0-6_T1-8_t0avg48_n1-100.csv"
)

PYERRORS_INPUT = (
    RESULTS /
    "laplace_Veff_R0-6_T1-7_t0avg48_n1-100_pyerrors.csv"
)

WINDOW_OUTPUT = (
    RESULTS /
    "laplace_Veff_R1-6_T1-7_t0avg48_correlated_plateau_windows.csv"
)

POTENTIAL_OUTPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg48_correlated.csv"
)

SAMPLES_OUTPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg48_jackknife_samples.csv"
)

COVARIANCE_OUTPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg48_covariance.csv"
)

CORRELATION_OUTPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg48_correlation.csv"
)

PLATEAU_FIGURE = (
    FIGURES /
    "laplace_Veff_correlated_plateau_T3-5_t0avg48_R1-6.png"
)

POTENTIAL_FIGURE = (
    FIGURES /
    "laplace_static_potential_correlated_T3-5_t0avg48_R1-6.png"
)

NCONFIG = 100
R_VALUES = np.arange(1, 7, dtype=int)
T_VALUES = np.arange(1, 8, dtype=int)
SELECTED_WINDOW = (3, 5)

EXPECTED_NSRC = 48 * 24**3
RCOND = 1.0e-12

KEY_WINDOWS = [
    (2, 4),
    (3, 5),
    (3, 6),
    (4, 6),
    (4, 7),
    (5, 7),
]


def jackknife_covariance(samples):
    samples = np.asarray(samples, dtype=float)

    if samples.ndim == 1:
        samples = samples[:, None]

    center = samples.mean(axis=0)
    deviations = samples - center

    return (
        (NCONFIG - 1.0) / NCONFIG
        * deviations.T
        @ deviations
    )


def symmetric_pseudoinverse(matrix):
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)

    scale = np.max(np.abs(eigenvalues))

    if not np.isfinite(scale) or scale <= 0.0:
        raise RuntimeError(
            "Covariance matrix has no positive scale"
        )

    keep = eigenvalues > RCOND * scale

    if not np.any(keep):
        raise RuntimeError(
            "Covariance pseudoinverse retained no modes"
        )

    inverse = (
        eigenvectors[:, keep]
        @ np.diag(1.0 / eigenvalues[keep])
        @ eigenvectors[:, keep].T
    )

    return (
        inverse,
        int(np.count_nonzero(keep)),
        eigenvalues,
    )


def load_inputs():
    raw = pd.read_csv(RAW_INPUT)
    pyerrors = pd.read_csv(PYERRORS_INPUT)

    expected_raw_columns = [
        "cfg",
        "T",
        "R",
        "Nsrc",
        "Re",
        "Im",
    ]

    if raw.columns.tolist() != expected_raw_columns:
        raise RuntimeError(
            f"Unexpected raw columns: {raw.columns.tolist()}"
        )

    if len(raw) != 5600:
        raise RuntimeError(
            f"Expected 5600 raw rows, found {len(raw)}"
        )

    if set(raw["cfg"]) != set(range(1, NCONFIG + 1)):
        raise RuntimeError("Configuration grid mismatch")

    if set(raw["T"]) != set(range(1, 9)):
        raise RuntimeError("Temporal grid mismatch")

    if set(raw["R"]) != set(range(0, 7)):
        raise RuntimeError("Separation grid mismatch")

    if not (raw["Nsrc"] == EXPECTED_NSRC).all():
        raise RuntimeError("Unexpected Nsrc values")

    if raw.duplicated(["cfg", "T", "R"]).any():
        raise RuntimeError("Duplicate cfg,T,R rows")

    if not np.isfinite(
        raw[["Re", "Im"]].to_numpy(dtype=float)
    ).all():
        raise RuntimeError("Non-finite raw correlators")

    if len(pyerrors) != 49:
        raise RuntimeError(
            f"Expected 49 pyerrors rows, "
            f"found {len(pyerrors)}"
        )

    return raw, pyerrors


def build_veff_jackknife(raw):
    analysis = {}

    for R in R_VALUES:
        matrix = (
            raw.loc[raw["R"] == R]
            .pivot(
                index="cfg",
                columns="T",
                values="Re",
            )
            .sort_index()
            .reindex(columns=range(1, 9))
            .to_numpy(dtype=float)
        )

        if matrix.shape != (NCONFIG, 8):
            raise RuntimeError(
                f"R={R}: unexpected matrix shape "
                f"{matrix.shape}"
            )

        sums = matrix.sum(axis=0)
        central_means = sums / NCONFIG

        leave_one_out_means = (
            sums[None, :] - matrix
        ) / (NCONFIG - 1)

        central_veff = np.full(7, np.nan)
        jackknife_veff = np.full(
            (NCONFIG, 7),
            np.nan,
        )

        status = []

        for index, T in enumerate(T_VALUES):
            central_valid = (
                central_means[index] > 0.0
                and central_means[index + 1] > 0.0
            )

            jackknife_valid = (
                np.all(
                    leave_one_out_means[:, index]
                    > 0.0
                )
                and np.all(
                    leave_one_out_means[:, index + 1]
                    > 0.0
                )
            )

            if not central_valid:
                status.append(
                    "nonpositive_central_mean"
                )
                continue

            if not jackknife_valid:
                status.append(
                    "nonpositive_jackknife_mean"
                )
                continue

            central_veff[index] = np.log(
                central_means[index]
                / central_means[index + 1]
            )

            jackknife_veff[:, index] = np.log(
                leave_one_out_means[:, index]
                / leave_one_out_means[:, index + 1]
            )

            status.append("ok")

        analysis[int(R)] = {
            "central": central_veff,
            "samples": jackknife_veff,
            "status": status,
        }

    return analysis


def fit_window(analysis, R, Tmin, Tmax):
    indices = np.arange(
        Tmin - 1,
        Tmax,
        dtype=int,
    )

    statuses = [
        analysis[R]["status"][index]
        for index in indices
    ]

    base = {
        "R": R,
        "Tmin": Tmin,
        "Tmax": Tmax,
        "fit_window": f"{Tmin}-{Tmax}",
        "n_points": len(indices),
    }

    if any(status != "ok" for status in statuses):
        return {
            **base,
            "plateau_value": np.nan,
            "err_plateau_jackknife": np.nan,
            "gls_error": np.nan,
            "chi2": np.nan,
            "dof": 0,
            "chi2_per_dof": np.nan,
            "p_value": np.nan,
            "covariance_rank": 0,
            "status": "invalid_veff_in_window",
            "_samples": None,
            "_weights": None,
        }

    central = analysis[R]["central"][indices]
    samples = analysis[R]["samples"][:, indices]

    covariance = jackknife_covariance(samples)

    inverse, rank, _ = symmetric_pseudoinverse(
        covariance
    )

    ones = np.ones(len(indices), dtype=float)
    denominator = float(ones @ inverse @ ones)

    if denominator <= 0.0:
        return {
            **base,
            "plateau_value": np.nan,
            "err_plateau_jackknife": np.nan,
            "gls_error": np.nan,
            "chi2": np.nan,
            "dof": 0,
            "chi2_per_dof": np.nan,
            "p_value": np.nan,
            "covariance_rank": rank,
            "status": "invalid_gls_denominator",
            "_samples": None,
            "_weights": None,
        }

    weights = inverse @ ones / denominator

    plateau_value = float(weights @ central)
    plateau_samples = samples @ weights

    plateau_covariance = jackknife_covariance(
        plateau_samples
    )

    plateau_error = float(
        np.sqrt(plateau_covariance[0, 0])
    )

    gls_error = float(
        np.sqrt(1.0 / denominator)
    )

    residual = central - plateau_value

    chi2 = float(
        residual @ inverse @ residual
    )

    dof = max(rank - 1, 0)

    chi2_per_dof = (
        chi2 / dof
        if dof > 0
        else np.nan
    )

    p_value = (
        float(
            chi2_distribution.sf(
                chi2,
                dof,
            )
        )
        if dof > 0
        else np.nan
    )

    return {
        **base,
        "plateau_value": plateau_value,
        "err_plateau_jackknife": plateau_error,
        "gls_error": gls_error,
        "chi2": chi2,
        "dof": dof,
        "chi2_per_dof": chi2_per_dof,
        "p_value": p_value,
        "covariance_rank": rank,
        "status": "ok",
        "_samples": plateau_samples,
        "_weights": weights,
    }


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)

    raw, pyerrors = load_inputs()
    analysis = build_veff_jackknife(raw)

    window_rows = []
    selected_rows = []
    selected_samples = []

    for R in R_VALUES:
        for Tmin in range(1, 6):
            for Tmax in range(Tmin + 2, 8):
                result = fit_window(
                    analysis,
                    int(R),
                    Tmin,
                    Tmax,
                )

                public_result = {
                    key: value
                    for key, value in result.items()
                    if not key.startswith("_")
                }

                window_rows.append(public_result)

                if (
                    Tmin,
                    Tmax,
                ) == SELECTED_WINDOW:
                    selected_rows.append(
                        public_result
                    )

                    if result["_samples"] is None:
                        raise RuntimeError(
                            f"Selected window invalid "
                            f"for R={R}"
                        )

                    selected_samples.append(
                        result["_samples"]
                    )

    windows = pd.DataFrame(window_rows)
    selected = pd.DataFrame(selected_rows)

    windows.to_csv(
        WINDOW_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    selected.to_csv(
        POTENTIAL_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    sample_matrix = np.column_stack(
        selected_samples
    )

    sample_frame = pd.DataFrame(
        sample_matrix,
        columns=[
            f"R{R}"
            for R in R_VALUES
        ],
    )

    sample_frame.insert(
        0,
        "omitted_cfg",
        np.arange(1, NCONFIG + 1),
    )

    sample_frame.to_csv(
        SAMPLES_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    potential_covariance = jackknife_covariance(
        sample_matrix
    )

    diagonal = np.sqrt(
        np.diag(potential_covariance)
    )

    potential_correlation = (
        potential_covariance
        / np.outer(diagonal, diagonal)
    )

    labels = [
        f"R{R}"
        for R in R_VALUES
    ]

    pd.DataFrame(
        potential_covariance,
        index=labels,
        columns=labels,
    ).to_csv(
        COVARIANCE_OUTPUT,
        lineterminator="\n",
    )

    pd.DataFrame(
        potential_correlation,
        index=labels,
        columns=labels,
    ).to_csv(
        CORRELATION_OUTPUT,
        lineterminator="\n",
    )

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(14.5, 8.5),
        sharex=True,
    )

    for axis, R in zip(
        axes.flat,
        R_VALUES,
    ):
        subset = (
            pyerrors.loc[pyerrors["R"] == R]
            .sort_values("T")
        )

        valid = (
            subset["status"].eq("ok")
            & np.isfinite(subset["Veff"])
            & np.isfinite(
                subset["err_pyerrors"]
            )
        )

        axis.errorbar(
            subset.loc[valid, "T"],
            subset.loc[valid, "Veff"],
            yerr=subset.loc[
                valid,
                "err_pyerrors",
            ],
            fmt="o",
            capsize=3,
            label="Gamma method",
        )

        result = selected.loc[
            selected["R"] == R
        ].iloc[0]

        value = result["plateau_value"]
        error = result[
            "err_plateau_jackknife"
        ]

        axis.axvspan(
            SELECTED_WINDOW[0] - 0.5,
            SELECTED_WINDOW[1] + 0.5,
            alpha=0.10,
            label="Selected fit window",
        )

        axis.axhspan(
            value - error,
            value + error,
            alpha=0.18,
            label=r"Plateau $\pm1\sigma$",
        )

        axis.axhline(
            value,
            linestyle="--",
            linewidth=1.3,
            label="Plateau central value",
        )

        axis.set_title(
            rf"$R/a={R}$, "
            rf"$p={result['p_value']:.3f}$"
        )

        axis.set_xticks(T_VALUES)
        axis.grid(alpha=0.3)

    for axis in axes[-1, :]:
        axis.set_xlabel(r"$T/a$")

    for axis in axes[:, 0]:
        axis.set_ylabel(r"$aV_{\mathrm{eff}}(R,T)$")

    handles, labels_legend = axes[0, 0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels_legend,
        loc="upper center",
        ncol=4,
        bbox_to_anchor=(0.5, 0.995),
    )

    figure.suptitle(
        "Laplace effective potential and correlated "
        r"$T/a=3$-$5$ plateaux",
        y=1.035,
    )

    figure.tight_layout()

    figure.savefig(
        PLATEAU_FIGURE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(figure)

    figure, axis = plt.subplots(
        figsize=(7.2, 5.0)
    )

    axis.errorbar(
        selected["R"],
        selected["plateau_value"],
        yerr=selected[
            "err_plateau_jackknife"
        ],
        fmt="o",
        capsize=4,
        label=r"Correlated $T/a=3$-$5$ plateau",
    )

    axis.plot(
        selected["R"],
        selected["plateau_value"],
        linestyle="--",
        alpha=0.7,
        label="Guide to the eye",
    )

    axis.set_xlabel(r"$R/a$")
    axis.set_ylabel(r"$aV(R)$")
    axis.set_xticks(R_VALUES)
    axis.grid(alpha=0.3)
    axis.legend()

    axis.set_title(
        "Laplace static potential from "
        r"correlated $T/a=3$-$5$ plateaux"
    )

    figure.tight_layout()

    figure.savefig(
        POTENTIAL_FIGURE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(figure)

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
            "status",
        ]].to_string(
            index=False,
            float_format=lambda value:
                f"{value:.8g}",
        )
    )

    key = windows[
        windows.apply(
            lambda row: (
                int(row["Tmin"]),
                int(row["Tmax"]),
            ) in KEY_WINDOWS,
            axis=1,
        )
    ]

    print()
    print("Key stability windows:")
    print(
        key[[
            "R",
            "Tmin",
            "Tmax",
            "plateau_value",
            "err_plateau_jackknife",
            "chi2_per_dof",
            "p_value",
            "status",
        ]].to_string(
            index=False,
            float_format=lambda value:
                f"{value:.8g}",
        )
    )

    print()
    print("Potential covariance matrix:")
    print(
        np.array2string(
            potential_covariance,
            precision=10,
            suppress_small=False,
        )
    )

    print()
    print("Wrote:")
    for path in [
        WINDOW_OUTPUT,
        POTENTIAL_OUTPUT,
        SAMPLES_OUTPUT,
        COVARIANCE_OUTPUT,
        CORRELATION_OUTPUT,
        PLATEAU_FIGURE,
        POTENTIAL_FIGURE,
    ]:
        print(path)


if __name__ == "__main__":
    main()
