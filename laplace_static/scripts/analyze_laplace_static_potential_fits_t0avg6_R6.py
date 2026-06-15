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

POTENTIAL_INPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg6_correlated.csv"
)

SAMPLES_INPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg6_jackknife_samples.csv"
)

COVARIANCE_INPUT = (
    RESULTS /
    "laplace_static_potential_R1-6_T3-5_t0avg6_covariance.csv"
)

SUMMARY_OUTPUT = (
    RESULTS /
    "laplace_static_potential_fit_stability_t0avg6_R6.csv"
)

PARAMETER_SAMPLES_OUTPUT = (
    RESULTS /
    "laplace_static_potential_fit_parameters_t0avg6_R6_jackknife.csv"
)

MAIN_CURVE_OUTPUT = (
    RESULTS /
    "laplace_static_potential_cornell_R1-6_t0avg6_curve.csv"
)

MAIN_FIGURE = (
    FIGURES /
    "laplace_static_potential_cornell_R1-6_t0avg6.png"
)

STABILITY_FIGURE = (
    FIGURES /
    "laplace_static_potential_fit_stability_t0avg6_R6.png"
)

NCONFIG = 100
ALL_R = np.arange(1, 7, dtype=int)
RCOND = 1.0e-12

FIT_SPECS = [
    {
        "model": "cornell",
        "fit_label": "R1-6",
        "Rmin": 1,
        "Rmax": 6,
    },
    {
        "model": "cornell",
        "fit_label": "R1-5",
        "Rmin": 1,
        "Rmax": 5,
    },
    {
        "model": "cornell",
        "fit_label": "R1-4",
        "Rmin": 1,
        "Rmax": 4,
    },
    {
        "model": "cornell",
        "fit_label": "R2-6",
        "Rmin": 2,
        "Rmax": 6,
    },
    {
        "model": "cornell",
        "fit_label": "R2-5",
        "Rmin": 2,
        "Rmax": 5,
    },
    {
        "model": "cornell",
        "fit_label": "R3-6",
        "Rmin": 3,
        "Rmax": 6,
    },
    {
        "model": "linear",
        "fit_label": "R3-6",
        "Rmin": 3,
        "Rmax": 6,
    },
    {
        "model": "linear",
        "fit_label": "R3-5",
        "Rmin": 3,
        "Rmax": 5,
    },
    {
        "model": "linear",
        "fit_label": "R4-6",
        "Rmin": 4,
        "Rmax": 6,
    },
]


def jackknife_covariance(samples):
    samples = np.asarray(samples, dtype=float)

    if samples.ndim == 1:
        samples = samples[:, None]

    if samples.shape[0] != NCONFIG:
        raise RuntimeError(
            f"Expected {NCONFIG} jackknife samples, "
            f"found {samples.shape[0]}"
        )

    deviations = samples - samples.mean(axis=0)

    return (
        (NCONFIG - 1.0) / NCONFIG
        * deviations.T
        @ deviations
    )


def jackknife_error(samples):
    samples = np.asarray(samples, dtype=float)
    deviations = samples - samples.mean(axis=0)

    return np.sqrt(
        (NCONFIG - 1.0) / NCONFIG
        * np.sum(deviations**2, axis=0)
    )


def symmetric_pseudoinverse(matrix):
    matrix = np.asarray(matrix, dtype=float)

    if matrix.shape[0] != matrix.shape[1]:
        raise RuntimeError(
            f"Expected square matrix, got {matrix.shape}"
        )

    asymmetry = float(
        np.max(
            np.abs(matrix - matrix.T)
        )
    )

    matrix_scale = float(
        np.max(np.abs(matrix))
    )

    symmetry_tolerance = (
        1.0e-12
        * max(
            matrix_scale,
            np.finfo(float).tiny,
        )
    )

    if (
        not np.isfinite(asymmetry)
        or asymmetry > symmetry_tolerance
    ):
        raise RuntimeError(
            "Matrix is not symmetric within "
            "floating-point tolerance: "
            f"max asymmetry={asymmetry:.6e}, "
            f"tolerance={symmetry_tolerance:.6e}"
        )

    # Remove harmless floating-point antisymmetric roundoff.
    matrix = 0.5 * (
        matrix + matrix.T
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        matrix
    )

    scale = np.max(np.abs(eigenvalues))

    if not np.isfinite(scale) or scale <= 0.0:
        raise RuntimeError(
            "Matrix has no positive scale"
        )

    keep = eigenvalues > RCOND * scale

    if not np.any(keep):
        raise RuntimeError(
            "Pseudoinverse retained no eigenmodes"
        )

    inverse = (
        eigenvectors[:, keep]
        @ np.diag(1.0 / eigenvalues[keep])
        @ eigenvectors[:, keep].T
    )

    retained = eigenvalues[keep]

    condition = float(
        retained.max() / retained.min()
    )

    return (
        inverse,
        int(np.count_nonzero(keep)),
        eigenvalues,
        condition,
    )


def design_matrix(R, model):
    R = np.asarray(R, dtype=float)

    if model == "cornell":
        return np.column_stack([
            np.ones_like(R),
            R,
            -1.0 / R,
        ])

    if model == "linear":
        return np.column_stack([
            np.ones_like(R),
            R,
        ])

    raise RuntimeError(
        f"Unknown model: {model}"
    )


def parameter_names(model):
    if model == "cornell":
        return [
            "V0",
            "sigma",
            "alpha",
        ]

    if model == "linear":
        return [
            "V0",
            "sigma",
        ]

    raise RuntimeError(
        f"Unknown model: {model}"
    )


def load_inputs():
    potential = pd.read_csv(
        POTENTIAL_INPUT
    )

    samples = pd.read_csv(
        SAMPLES_INPUT
    )

    covariance_frame = pd.read_csv(
        COVARIANCE_INPUT,
        index_col=0,
    )

    if len(potential) != 6:
        raise RuntimeError(
            f"Expected six potential rows, "
            f"found {len(potential)}"
        )

    if potential["R"].tolist() != ALL_R.tolist():
        raise RuntimeError(
            "Potential R ordering mismatch"
        )

    if not (
        potential["status"] == "ok"
    ).all():
        raise RuntimeError(
            "At least one selected plateau is invalid"
        )

    expected_sample_columns = [
        "omitted_cfg",
        "R1",
        "R2",
        "R3",
        "R4",
        "R5",
        "R6",
    ]

    if samples.columns.tolist() != expected_sample_columns:
        raise RuntimeError(
            f"Unexpected sample columns: "
            f"{samples.columns.tolist()}"
        )

    if samples.shape != (NCONFIG, 7):
        raise RuntimeError(
            f"Unexpected samples shape: "
            f"{samples.shape}"
        )

    if samples["omitted_cfg"].tolist() != list(
        range(1, NCONFIG + 1)
    ):
        raise RuntimeError(
            "Jackknife omitted-configuration "
            "ordering mismatch"
        )

    expected_covariance_labels = [
        "R1",
        "R2",
        "R3",
        "R4",
        "R5",
        "R6",
    ]

    if (
        covariance_frame.index.tolist()
        != expected_covariance_labels
    ):
        raise RuntimeError(
            "Covariance row labels mismatch"
        )

    if (
        covariance_frame.columns.tolist()
        != expected_covariance_labels
    ):
        raise RuntimeError(
            "Covariance column labels mismatch"
        )

    covariance = covariance_frame.to_numpy(
        dtype=float
    )

    sample_matrix = samples[
        expected_covariance_labels
    ].to_numpy(dtype=float)

    central = potential[
        "plateau_value"
    ].to_numpy(dtype=float)

    reconstructed_covariance = (
        jackknife_covariance(sample_matrix)
    )

    if not np.allclose(
        covariance,
        reconstructed_covariance,
        rtol=1.0e-10,
        atol=1.0e-14,
    ):
        maximum = np.max(
            np.abs(
                covariance
                - reconstructed_covariance
            )
        )

        raise RuntimeError(
            "Stored covariance does not match "
            "jackknife samples; "
            f"maximum difference={maximum}"
        )

    return (
        potential,
        central,
        sample_matrix,
        covariance,
    )


def fit_model(
    central,
    sample_matrix,
    covariance,
    model,
    Rmin,
    Rmax,
):
    indices = np.where(
        (ALL_R >= Rmin)
        & (ALL_R <= Rmax)
    )[0]

    R = ALL_R[indices].astype(float)
    values = central[indices]
    samples = sample_matrix[:, indices]
    selected_covariance = covariance[
        np.ix_(indices, indices)
    ]

    X = design_matrix(
        R,
        model,
    )

    names = parameter_names(model)
    n_parameters = len(names)
    n_points = len(R)

    if n_points < n_parameters:
        raise RuntimeError(
            f"{model} R={Rmin}..{Rmax}: "
            "fewer points than parameters"
        )

    (
        covariance_inverse,
        covariance_rank,
        covariance_eigenvalues,
        covariance_condition,
    ) = symmetric_pseudoinverse(
        selected_covariance
    )

    normal_matrix = (
        X.T
        @ covariance_inverse
        @ X
    )

    (
        normal_inverse,
        normal_rank,
        normal_eigenvalues,
        normal_condition,
    ) = symmetric_pseudoinverse(
        normal_matrix
    )

    if normal_rank != n_parameters:
        raise RuntimeError(
            f"{model} R={Rmin}..{Rmax}: "
            f"normal-matrix rank {normal_rank}, "
            f"expected {n_parameters}"
        )

    fit_operator = (
        normal_inverse
        @ X.T
        @ covariance_inverse
    )

    parameters = (
        fit_operator @ values
    )

    parameter_samples = (
        samples @ fit_operator.T
    )

    parameter_covariance = (
        jackknife_covariance(
            parameter_samples
        )
    )

    parameter_errors = np.sqrt(
        np.diag(parameter_covariance)
    )

    if not np.allclose(
        parameter_covariance,
        normal_inverse,
        rtol=1.0e-8,
        atol=1.0e-12,
    ):
        maximum = np.max(
            np.abs(
                parameter_covariance
                - normal_inverse
            )
        )

        raise RuntimeError(
            f"{model} R={Rmin}..{Rmax}: "
            "analytic and jackknife parameter "
            "covariances disagree; "
            f"maximum difference={maximum}"
        )

    fitted_values = X @ parameters
    residuals = values - fitted_values

    chi2 = float(
        residuals
        @ covariance_inverse
        @ residuals
    )

    dof = covariance_rank - n_parameters

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

    parameter_map = dict(
        zip(names, parameters)
    )

    error_map = dict(
        zip(names, parameter_errors)
    )

    sample_map = {
        name: parameter_samples[:, index]
        for index, name in enumerate(names)
    }

    sigma = parameter_map["sigma"]
    sigma_samples = sample_map["sigma"]

    if (
        sigma > 0.0
        and np.all(sigma_samples > 0.0)
    ):
        a_sqrt_sigma = float(
            np.sqrt(sigma)
        )

        a_sqrt_sigma_samples = np.sqrt(
            sigma_samples
        )

        err_a_sqrt_sigma = float(
            jackknife_error(
                a_sqrt_sigma_samples
            )
        )
    else:
        a_sqrt_sigma = np.nan
        a_sqrt_sigma_samples = np.full(
            NCONFIG,
            np.nan,
        )
        err_a_sqrt_sigma = np.nan

    result = {
        "model": model,
        "Rmin": int(Rmin),
        "Rmax": int(Rmax),
        "n_points": int(n_points),
        "n_parameters": int(n_parameters),
        "covariance_rank": int(
            covariance_rank
        ),
        "covariance_condition": (
            covariance_condition
        ),
        "normal_condition": normal_condition,
        "chi2": chi2,
        "dof": int(dof),
        "chi2_per_dof": chi2_per_dof,
        "p_value": p_value,
        "V0": float(
            parameter_map["V0"]
        ),
        "err_V0": float(
            error_map["V0"]
        ),
        "sigma": float(sigma),
        "err_sigma": float(
            error_map["sigma"]
        ),
        "alpha": (
            float(parameter_map["alpha"])
            if model == "cornell"
            else np.nan
        ),
        "err_alpha": (
            float(error_map["alpha"])
            if model == "cornell"
            else np.nan
        ),
        "a_sqrt_sigma": a_sqrt_sigma,
        "err_a_sqrt_sigma": (
            err_a_sqrt_sigma
        ),
        "status": "ok",
    }

    sample_result = {
        "V0": sample_map["V0"],
        "sigma": sample_map["sigma"],
        "alpha": (
            sample_map["alpha"]
            if model == "cornell"
            else np.full(
                NCONFIG,
                np.nan,
            )
        ),
        "a_sqrt_sigma": (
            a_sqrt_sigma_samples
        ),
    }

    diagnostics = {
        "indices": indices,
        "R": R,
        "values": values,
        "covariance": selected_covariance,
        "covariance_eigenvalues":
            covariance_eigenvalues,
        "normal_eigenvalues":
            normal_eigenvalues,
        "parameters": parameters,
        "parameter_samples":
            parameter_samples,
        "fitted_values": fitted_values,
        "residuals": residuals,
    }

    return (
        result,
        sample_result,
        diagnostics,
    )


def main():
    FIGURES.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        potential,
        central,
        sample_matrix,
        covariance,
    ) = load_inputs()

    summary_records = []
    sample_records = []
    diagnostics_by_key = {}

    for specification in FIT_SPECS:
        (
            result,
            sample_result,
            diagnostics,
        ) = fit_model(
            central=central,
            sample_matrix=sample_matrix,
            covariance=covariance,
            model=specification["model"],
            Rmin=specification["Rmin"],
            Rmax=specification["Rmax"],
        )

        result["fit_label"] = (
            specification["fit_label"]
        )

        summary_records.append(result)

        key = (
            specification["model"],
            specification["fit_label"],
        )

        diagnostics_by_key[key] = diagnostics

        for omitted_cfg in range(
            1,
            NCONFIG + 1,
        ):
            index = omitted_cfg - 1

            sample_records.append({
                "model":
                    specification["model"],
                "fit_label":
                    specification["fit_label"],
                "Rmin":
                    specification["Rmin"],
                "Rmax":
                    specification["Rmax"],
                "omitted_cfg": omitted_cfg,
                "V0":
                    sample_result["V0"][index],
                "sigma":
                    sample_result["sigma"][index],
                "alpha":
                    sample_result["alpha"][index],
                "a_sqrt_sigma":
                    sample_result[
                        "a_sqrt_sigma"
                    ][index],
            })

    summary = pd.DataFrame(
        summary_records
    )

    summary = summary[[
        "model",
        "fit_label",
        "Rmin",
        "Rmax",
        "n_points",
        "n_parameters",
        "covariance_rank",
        "covariance_condition",
        "normal_condition",
        "chi2",
        "dof",
        "chi2_per_dof",
        "p_value",
        "V0",
        "err_V0",
        "sigma",
        "err_sigma",
        "alpha",
        "err_alpha",
        "a_sqrt_sigma",
        "err_a_sqrt_sigma",
        "status",
    ]]

    summary.to_csv(
        SUMMARY_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    parameter_samples = pd.DataFrame(
        sample_records
    )

    parameter_samples.to_csv(
        PARAMETER_SAMPLES_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    main_key = (
        "cornell",
        "R1-6",
    )

    main_summary = summary[
        (summary["model"] == main_key[0])
        & (
            summary["fit_label"]
            == main_key[1]
        )
    ].iloc[0]

    main_samples = parameter_samples[
        (
            parameter_samples["model"]
            == main_key[0]
        )
        & (
            parameter_samples["fit_label"]
            == main_key[1]
        )
    ].sort_values("omitted_cfg")

    grid = np.linspace(
        1.0,
        6.0,
        401,
    )

    grid_design = design_matrix(
        grid,
        "cornell",
    )

    main_parameters = np.array([
        main_summary["V0"],
        main_summary["sigma"],
        main_summary["alpha"],
    ])

    curve = (
        grid_design @ main_parameters
    )

    main_sample_matrix = main_samples[[
        "V0",
        "sigma",
        "alpha",
    ]].to_numpy(dtype=float)

    curve_samples = (
        main_sample_matrix
        @ grid_design.T
    )

    curve_error = jackknife_error(
        curve_samples
    )

    curve_frame = pd.DataFrame({
        "R_over_a": grid,
        "aV": curve,
        "err_aV_jackknife": curve_error,
        "lower_1sigma":
            curve - curve_error,
        "upper_1sigma":
            curve + curve_error,
    })

    curve_frame.to_csv(
        MAIN_CURVE_OUTPUT,
        index=False,
        lineterminator="\n",
    )

    figure, axis = plt.subplots(
        figsize=(7.4, 5.2)
    )

    data_error = np.sqrt(
        np.diag(covariance)
    )

    axis.errorbar(
        ALL_R,
        central,
        yerr=data_error,
        fmt="o",
        capsize=4,
        label=r"Correlated $T/a=3$-$5$ plateaux",
    )

    axis.plot(
        grid,
        curve,
        label="Correlated Cornell fit",
    )

    axis.fill_between(
        grid,
        curve - curve_error,
        curve + curve_error,
        alpha=0.20,
        label=r"Jackknife $1\sigma$ band",
    )

    axis.set_xlabel(r"$R/a$")
    axis.set_ylabel(r"$aV(R)$")
    axis.set_xticks(ALL_R)
    axis.grid(alpha=0.3)
    axis.legend()

    axis.set_title(
        "Laplace static potential: "
        "correlated Cornell fit\n"
        rf"$\chi^2/\mathrm{{dof}}="
        rf"{main_summary['chi2']:.3f}/"
        rf"{int(main_summary['dof'])}$, "
        rf"$p={main_summary['p_value']:.3f}$"
    )

    figure.tight_layout()

    figure.savefig(
        MAIN_FIGURE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(figure)

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.8),
    )

    labels = [
        f"{row.model}\n{row.fit_label}"
        for row in summary.itertuples()
    ]

    positions = np.arange(
        len(summary)
    )

    axes[0].errorbar(
        positions,
        summary["V0"],
        yerr=summary["err_V0"],
        fmt="o",
        capsize=3,
    )

    axes[0].set_ylabel(r"$V_0$")
    axes[0].set_title("Additive constant")

    axes[1].errorbar(
        positions,
        summary["sigma"],
        yerr=summary["err_sigma"],
        fmt="o",
        capsize=3,
    )

    axes[1].set_ylabel(r"$a^2\sigma$")
    axes[1].set_title("String-tension parameter")

    cornell_mask = (
        summary["model"] == "cornell"
    )

    axes[2].errorbar(
        positions[cornell_mask],
        summary.loc[
            cornell_mask,
            "alpha",
        ],
        yerr=summary.loc[
            cornell_mask,
            "err_alpha",
        ],
        fmt="o",
        capsize=3,
    )

    axes[2].set_ylabel(r"$\alpha$")
    axes[2].set_title("Coulomb coefficient")

    for axis in axes:
        axis.set_xticks(positions)
        axis.set_xticklabels(
            labels,
            rotation=55,
            ha="right",
        )
        axis.grid(alpha=0.3)

        axis.axvline(
            5.5,
            linestyle="--",
            linewidth=1.0,
        )

    figure.suptitle(
        "Correlated static-potential "
        "fit stability"
    )

    figure.tight_layout()

    figure.savefig(
        STABILITY_FIGURE,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(figure)

    print("Correlated potential-fit summary:")
    print(
        summary.to_string(
            index=False,
            float_format=lambda value:
                f"{value:.8g}",
        )
    )

    print()
    print("Main Cornell R/a = 1..6 fit:")
    print(
        main_summary.to_string()
    )

    print()
    print("Wrote:")
    for path in [
        SUMMARY_OUTPUT,
        PARAMETER_SAMPLES_OUTPUT,
        MAIN_CURVE_OUTPUT,
        MAIN_FIGURE,
        STABILITY_FIGURE,
    ]:
        print(path)


if __name__ == "__main__":
    main()
