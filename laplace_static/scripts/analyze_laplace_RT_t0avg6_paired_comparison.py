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
FIGDIR = RESULTS / "figures"

FIXED_INPUT = (
    RESULTS /
    "spatial_avg_RT_R0-4_T1-8_n1-100.csv"
)

AVG6_INPUT = (
    RESULTS /
    "spatial_avg_RT_R0-4_T1-8_t0avg6_n1-100.csv"
)

FIXED_SELECTED = (
    RESULTS /
    "laplace_static_potential_R1-4_T3-5_correlated.csv"
)

AVG6_SELECTED = (
    RESULTS /
    "laplace_static_potential_R1-4_T3-5_t0avg6_correlated.csv"
)

OUTPUT_CSV = (
    RESULTS /
    "laplace_static_potential_fixed_vs_t0avg6_T3-5_paired.csv"
)

OUTPUT_SUMMARY = (
    RESULTS /
    "laplace_static_potential_fixed_vs_t0avg6_T3-5_paired_summary.txt"
)

OUTPUT_FIGURE = (
    FIGDIR /
    "laplace_static_potential_fixed_vs_t0avg6_T3-5_paired.png"
)

NCONFIG = 100
R_VALUES = np.arange(1, 5, dtype=int)
T_WINDOW = np.array([3, 4, 5], dtype=int)

COVARIANCE_RCOND = 1.0e-12


def load_raw(path, expected_nsrc):
    frame = pd.read_csv(path)

    expected_columns = [
        "cfg",
        "T",
        "R",
        "Nsrc",
        "Re",
        "Im",
    ]

    if frame.columns.tolist() != expected_columns:
        raise RuntimeError(
            f"{path}: unexpected columns "
            f"{frame.columns.tolist()}"
        )

    if len(frame) != NCONFIG * 8 * 5:
        raise RuntimeError(
            f"{path}: expected 4000 rows, "
            f"found {len(frame)}"
        )

    if set(frame["cfg"]) != set(range(1, NCONFIG + 1)):
        raise RuntimeError(
            f"{path}: configuration grid mismatch"
        )

    if set(frame["T"]) != set(range(1, 9)):
        raise RuntimeError(
            f"{path}: temporal grid mismatch"
        )

    if set(frame["R"]) != set(range(0, 5)):
        raise RuntimeError(
            f"{path}: separation grid mismatch"
        )

    if not (frame["Nsrc"] == expected_nsrc).all():
        raise RuntimeError(
            f"{path}: unexpected Nsrc values"
        )

    if frame.duplicated(["cfg", "T", "R"]).any():
        raise RuntimeError(
            f"{path}: duplicate cfg,T,R rows"
        )

    values = frame[["Re", "Im"]].to_numpy(
        dtype=float
    )

    if not np.isfinite(values).all():
        raise RuntimeError(
            f"{path}: non-finite correlator values"
        )

    return frame


def jackknife_covariance(samples):
    sample_mean = samples.mean(axis=0)
    centered = samples - sample_mean

    return (
        (NCONFIG - 1.0) / NCONFIG
        * centered.T
        @ centered
    )


def symmetric_pseudoinverse(matrix):
    eigenvalues, eigenvectors = np.linalg.eigh(
        matrix
    )

    scale = np.max(
        np.abs(eigenvalues)
    )

    if not np.isfinite(scale) or scale <= 0.0:
        raise RuntimeError(
            "Covariance matrix has no positive scale"
        )

    keep = (
        eigenvalues
        > COVARIANCE_RCOND * scale
    )

    if not np.any(keep):
        raise RuntimeError(
            "Covariance pseudoinverse retained no modes"
        )

    inverse = (
        eigenvectors[:, keep]
        @ np.diag(
            1.0 / eigenvalues[keep]
        )
        @ eigenvectors[:, keep].T
    )

    return (
        inverse,
        int(np.count_nonzero(keep)),
        eigenvalues,
    )


def analyse_dataset(frame):
    analysis = {}

    for R in R_VALUES:
        matrix = (
            frame.loc[frame["R"] == R]
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
                f"R={R}: unexpected correlator "
                f"matrix shape {matrix.shape}"
            )

        sums = matrix.sum(axis=0)

        central_means = (
            sums / NCONFIG
        )

        leave_one_out_means = (
            sums[None, :] - matrix
        ) / (NCONFIG - 1)

        if (
            np.any(central_means <= 0.0)
            or np.any(
                leave_one_out_means <= 0.0
            )
        ):
            raise RuntimeError(
                f"R={R}: nonpositive ensemble or "
                "leave-one-out mean"
            )

        central_veff = np.log(
            central_means[:-1]
            / central_means[1:]
        )

        jackknife_veff = np.log(
            leave_one_out_means[:, :-1]
            / leave_one_out_means[:, 1:]
        )

        window_indices = (
            T_WINDOW - 1
        )

        central_window = (
            central_veff[window_indices]
        )

        jackknife_window = (
            jackknife_veff[
                :,
                window_indices,
            ]
        )

        covariance = jackknife_covariance(
            jackknife_window
        )

        (
            inverse,
            rank,
            eigenvalues,
        ) = symmetric_pseudoinverse(
            covariance
        )

        ones = np.ones(
            len(T_WINDOW),
            dtype=float,
        )

        weights = (
            inverse @ ones
            / (ones @ inverse @ ones)
        )

        plateau_value = float(
            weights @ central_window
        )

        plateau_samples = (
            jackknife_window @ weights
        )

        plateau_error = float(
            np.sqrt(
                jackknife_covariance(
                    plateau_samples[:, None]
                )[0, 0]
            )
        )

        gls_error = float(
            np.sqrt(
                1.0
                / (ones @ inverse @ ones)
            )
        )

        analysis[int(R)] = {
            "value": plateau_value,
            "error": plateau_error,
            "gls_error": gls_error,
            "samples": plateau_samples,
            "weights": weights,
            "covariance_rank": rank,
            "covariance_eigenvalues":
                eigenvalues,
        }

    return analysis


def validate_against_selected(
    analysis,
    path,
):
    selected = (
        pd.read_csv(path)
        .set_index("R")
    )

    for R in R_VALUES:
        row = selected.loc[int(R)]

        np.testing.assert_allclose(
            analysis[int(R)]["value"],
            row["plateau_value"],
            rtol=5.0e-11,
            atol=5.0e-13,
        )

        np.testing.assert_allclose(
            analysis[int(R)]["error"],
            row["err_plateau_jackknife"],
            rtol=5.0e-11,
            atol=5.0e-13,
        )

    return selected


def main():
    FIGDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    fixed_frame = load_raw(
        FIXED_INPUT,
        24 ** 3,
    )

    avg6_frame = load_raw(
        AVG6_INPUT,
        6 * 24 ** 3,
    )

    fixed = analyse_dataset(
        fixed_frame
    )

    avg6 = analyse_dataset(
        avg6_frame
    )

    fixed_selected = validate_against_selected(
        fixed,
        FIXED_SELECTED,
    )

    avg6_selected = validate_against_selected(
        avg6,
        AVG6_SELECTED,
    )

    delta_values = []
    delta_samples = []
    records = []

    for R in R_VALUES:
        fixed_result = fixed[int(R)]
        avg6_result = avg6[int(R)]

        paired_samples = (
            avg6_result["samples"]
            - fixed_result["samples"]
        )

        delta_value = (
            avg6_result["value"]
            - fixed_result["value"]
        )

        delta_error = float(
            np.sqrt(
                jackknife_covariance(
                    paired_samples[:, None]
                )[0, 0]
            )
        )

        paired_covariance = (
            (NCONFIG - 1.0)
            / NCONFIG
            * np.sum(
                (
                    fixed_result["samples"]
                    - fixed_result[
                        "samples"
                    ].mean()
                )
                * (
                    avg6_result["samples"]
                    - avg6_result[
                        "samples"
                    ].mean()
                )
            )
        )

        correlation = (
            paired_covariance
            / (
                fixed_result["error"]
                * avg6_result["error"]
            )
        )

        records.append({
            "R": int(R),
            "fit_window": "3-5",
            "V_fixed":
                fixed_result["value"],
            "err_fixed":
                fixed_result["error"],
            "p_fixed": float(
                fixed_selected.loc[
                    int(R),
                    "p_value",
                ]
            ),
            "V_t0avg6":
                avg6_result["value"],
            "err_t0avg6":
                avg6_result["error"],
            "p_t0avg6": float(
                avg6_selected.loc[
                    int(R),
                    "p_value",
                ]
            ),
            "delta_V":
                delta_value,
            "err_delta_paired":
                delta_error,
            "z_delta":
                delta_value / delta_error,
            "error_reduction_factor":
                (
                    fixed_result["error"]
                    / avg6_result["error"]
                ),
            "corr_fixed_t0avg6":
                correlation,
        })

        delta_values.append(
            delta_value
        )

        delta_samples.append(
            paired_samples
        )

    comparison = pd.DataFrame.from_records(
        records
    )

    delta_values = np.asarray(
        delta_values,
        dtype=float,
    )

    delta_samples = np.column_stack(
        delta_samples
    )

    delta_covariance = jackknife_covariance(
        delta_samples
    )

    (
        delta_inverse,
        delta_rank,
        delta_eigenvalues,
    ) = symmetric_pseudoinverse(
        delta_covariance
    )

    global_chi2 = float(
        delta_values
        @ delta_inverse
        @ delta_values
    )

    global_p_value = float(
        chi2_distribution.sf(
            global_chi2,
            delta_rank,
        )
    )

    comparison.to_csv(
        OUTPUT_CSV,
        index=False,
        lineterminator="\n",
    )

    summary_lines = [
        (
            "Paired fixed-source versus six-source "
            "temporal-average comparison"
        ),
        "Fit window: T/a = 3-5",
        f"Configurations: {NCONFIG}",
        f"Covariance rank: {delta_rank}",
        f"Global chi2: {global_chi2:.12g}",
        f"Global dof: {delta_rank}",
        (
            "Global p-value: "
            f"{global_p_value:.12g}"
        ),
        "",
        "Paired delta covariance matrix:",
        np.array2string(
            delta_covariance,
            precision=12,
            suppress_small=False,
        ),
        "",
        (
            "Paired delta covariance "
            "eigenvalues:"
        ),
        np.array2string(
            delta_eigenvalues,
            precision=12,
            suppress_small=False,
        ),
        "",
    ]

    OUTPUT_SUMMARY.write_text(
        "\n".join(summary_lines)
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.7),
        gridspec_kw={
            "width_ratios": [
                1.35,
                1.0,
                1.0,
            ]
        },
    )

    x = comparison["R"].to_numpy(
        dtype=float
    )

    axes[0].errorbar(
        x - 0.045,
        comparison["V_fixed"],
        yerr=comparison["err_fixed"],
        fmt="o-",
        capsize=4,
        label=r"Fixed $t_0/a=0$",
    )

    axes[0].errorbar(
        x + 0.045,
        comparison["V_t0avg6"],
        yerr=comparison["err_t0avg6"],
        fmt="s-",
        capsize=4,
        label="Six-source average",
    )

    axes[0].set_xlabel(r"$R/a$")
    axes[0].set_ylabel(r"$a\,V(R)$")
    axes[0].set_xticks(R_VALUES)
    axes[0].set_title("Static potential")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].errorbar(
        x,
        comparison["delta_V"],
        yerr=comparison[
            "err_delta_paired"
        ],
        fmt="o",
        capsize=4,
    )

    axes[1].axhline(
        0.0,
        linestyle="--",
        linewidth=1.2,
    )

    axes[1].set_xlabel(r"$R/a$")
    axes[1].set_ylabel(
        r"$a\,[V_{\mathrm{avg6}}"
        r"-V_{\mathrm{fixed}}]$"
    )

    axes[1].set_xticks(R_VALUES)

    axes[1].set_title(
        "Paired differences\n"
        rf"$\chi^2={global_chi2:.2f}$, "
        rf"$\mathrm{{dof}}={delta_rank}$, "
        rf"$p={global_p_value:.3f}$"
    )

    axes[1].grid(alpha=0.3)

    axes[2].bar(
        x,
        comparison[
            "error_reduction_factor"
        ],
    )

    axes[2].axhline(
        1.0,
        linestyle="--",
        linewidth=1.2,
    )

    axes[2].set_xlabel(r"$R/a$")

    axes[2].set_ylabel(
        r"$\sigma_{\mathrm{fixed}}"
        r"/\sigma_{\mathrm{avg6}}$"
    )

    axes[2].set_xticks(R_VALUES)

    axes[2].set_title(
        "Plateau error reduction"
    )

    axes[2].grid(
        axis="y",
        alpha=0.3,
    )

    figure.suptitle(
        r"Laplace static potential: paired "
        r"$T/a=3$-$5$ comparison"
    )

    figure.tight_layout()

    figure.savefig(
        OUTPUT_FIGURE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(figure)

    print(
        "Paired T/a = 3-5 plateau comparison:"
    )

    print(
        comparison.to_string(
            index=False,
            float_format=lambda value:
                f"{value:.8g}",
        )
    )

    print()

    print(
        "Global correlated consistency test: "
        f"chi2={global_chi2:.8f}, "
        f"dof={delta_rank}, "
        f"p={global_p_value:.8f}"
    )

    print()
    print("Wrote:")
    print(OUTPUT_CSV)
    print(OUTPUT_SUMMARY)
    print(OUTPUT_FIGURE)


if __name__ == "__main__":
    main()
