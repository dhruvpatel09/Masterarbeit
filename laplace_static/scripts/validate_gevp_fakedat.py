#!/usr/bin/env python3
"""End-to-end validation of the production GEVP solver with supplied fake data.

The MATLAB input contains a 4x4 cell array named ``corr``.  Every cell holds
``(measurement, time)`` data.  This driver deliberately reuses the numerical
core from ``analyze_laplace_gaussian_gevp.py`` while comparing its generalized
eigenvalues against ``scipy.linalg.eigh``.

Statistical errors use delete-one-*block* jackknife samples.  This is important
because the supplied measurement histories are autocorrelated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.io import loadmat
from scipy.linalg import eigh
from scipy.stats import chi2


EXPECTED_SHA256 = (
    "0ad8a271954fa47559e7d307a2bb40983a5d97a63fee4829bcd34f27ddbff6be"
)
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    REPOSITORY_ROOT
    / "laplace_static/validation/fake_gevp/data/fakedat.mat"
)
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT
    / "laplace_static/validation/fake_gevp/results"
)
SOLVER_TOLERANCE = 5.0e-11
ORTHONORMALITY_TOLERANCE = 5.0e-10


def parse_integer_list(text: str) -> tuple[int, ...]:
    values = tuple(int(item.strip()) for item in text.split(",") if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected at least one integer")
    if any(value < 0 for value in values):
        raise argparse.ArgumentTypeError("values must be non-negative")
    if len(values) != len(set(values)):
        raise argparse.ArgumentTypeError("duplicate values are not allowed")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the production GEVP solver with fakedat.mat."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-sha256", default=EXPECTED_SHA256)
    parser.add_argument("--skip-checksum", action="store_true")
    parser.add_argument("--rcond", type=float, default=1.0e-10)
    parser.add_argument(
        "--reference-times",
        type=parse_integer_list,
        default=parse_integer_list("0,1,2,3,5,8,10"),
    )
    parser.add_argument("--primary-reference-time", type=int, default=5)
    parser.add_argument("--primary-block-size", type=int, default=10)
    parser.add_argument(
        "--block-sizes",
        type=parse_integer_list,
        default=parse_integer_list("1,2,5,10,20"),
    )
    parser.add_argument("--analysis-tmax", type=int, default=40)
    parser.add_argument("--plot-tmax", type=int, default=30)
    parser.add_argument("--fit-tmax", type=int, default=30)
    parser.add_argument("--minimum-fit-points", type=int, default=4)
    parser.add_argument("--maximum-fit-points", type=int, default=20)
    parser.add_argument("--acf-max-lag", type=int, default=20)
    parser.add_argument("--tau-window", type=int, default=10)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def production_module() -> Any:
    script_dir = Path(__file__).resolve().parent
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    import analyze_laplace_gaussian_gevp as production

    return production


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv_file(
    path: Path,
    fieldnames: Iterable[str],
    rows: Iterable[dict[str, Any]],
) -> None:
    fieldnames = list(fieldnames)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_fake_data(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    payload = loadmat(path, struct_as_record=False, squeeze_me=False)
    visible_keys = sorted(key for key in payload if not key.startswith("__"))
    if visible_keys != ["corr"]:
        raise ValueError(
            f"expected exactly the MATLAB variable 'corr'; found {visible_keys}"
        )

    corr = payload["corr"]
    if corr.shape[0] != corr.shape[1]:
        raise ValueError(f"corr must be square; found cell shape {corr.shape}")
    n_operators = int(corr.shape[0])
    if n_operators != 4:
        raise ValueError(f"expected a 4x4 cell array; found {corr.shape}")

    first = np.asarray(corr[0, 0])
    if first.ndim != 2:
        raise ValueError(f"corr{{1,1}} must be two-dimensional; found {first.shape}")
    n_measurements, n_times = map(int, first.shape)
    if (n_measurements, n_times) != (1000, 101):
        raise ValueError(
            "expected every cell to have shape (1000, 101); "
            f"first cell has {first.shape}"
        )

    matrices = np.empty(
        (n_measurements, n_times, n_operators, n_operators),
        dtype=np.complex128,
    )
    element_dtypes: set[str] = set()
    for source in range(n_operators):
        for sink in range(n_operators):
            values = np.asarray(corr[source, sink])
            if values.shape != (n_measurements, n_times):
                raise ValueError(
                    f"corr{{{source + 1},{sink + 1}}} has shape {values.shape}; "
                    f"expected {(n_measurements, n_times)}"
                )
            if not np.all(np.isfinite(values)):
                raise ValueError(
                    f"corr{{{source + 1},{sink + 1}}} contains non-finite values"
                )
            element_dtypes.add(str(values.dtype))
            matrices[:, :, source, sink] = values

    metadata = {
        "matlab_variables": visible_keys,
        "cell_shape": list(corr.shape),
        "cell_element_shape": [n_measurements, n_times],
        "n_measurements": n_measurements,
        "n_times": n_times,
        "n_operators": n_operators,
        "element_dtypes": sorted(element_dtypes),
        "maximum_input_imaginary_part": float(np.max(np.abs(matrices.imag))),
        "minimum_input_value": float(np.min(matrices.real)),
        "maximum_input_value": float(np.max(matrices.real)),
    }
    return matrices, metadata


def autocorrelation_summary(
    raw_matrices: np.ndarray,
    maximum_lag: int,
    tau_window: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n_measurements, n_times, n_operators, _ = raw_matrices.shape
    chains = raw_matrices.real.transpose(2, 3, 1, 0).reshape(
        n_operators * n_operators * n_times,
        n_measurements,
    )
    centered = chains - chains.mean(axis=1, keepdims=True)
    gamma_zero = np.mean(centered * centered, axis=1)
    if np.any(gamma_zero <= 0.0):
        raise ValueError("at least one primary input stream has zero variance")

    maximum_lag = min(maximum_lag, n_measurements - 1)
    rho = np.empty((chains.shape[0], maximum_lag + 1), dtype=float)
    rho[:, 0] = 1.0
    for lag in range(1, maximum_lag + 1):
        gamma = np.mean(centered[:, :-lag] * centered[:, lag:], axis=1)
        rho[:, lag] = gamma / gamma_zero

    rows: list[dict[str, Any]] = []
    for lag in range(maximum_lag + 1):
        values = rho[:, lag]
        rows.append(
            {
                "lag": lag,
                "rho_q05": float(np.quantile(values, 0.05)),
                "rho_median": float(np.median(values)),
                "rho_q95": float(np.quantile(values, 0.95)),
            }
        )

    effective_window = min(tau_window, maximum_lag)
    tau_values = 0.5 + np.sum(rho[:, 1 : effective_window + 1], axis=1)
    positive_tau = tau_values[tau_values > 0.0]
    if not positive_tau.size:
        raise ValueError("fixed-window integrated autocorrelation estimates are non-positive")
    summary = {
        "n_primary_streams": int(chains.shape[0]),
        "tau_window": effective_window,
        "tau_int_q05": float(np.quantile(positive_tau, 0.05)),
        "tau_int_median": float(np.median(positive_tau)),
        "tau_int_q95": float(np.quantile(positive_tau, 0.95)),
        "effective_sample_size_from_median_tau": float(
            n_measurements / (2.0 * np.median(positive_tau))
        ),
    }
    return rows, summary


def block_ensemble(matrices: np.ndarray, block_size: int) -> np.ndarray:
    if block_size <= 0:
        raise ValueError("block size must be positive")
    n_measurements = matrices.shape[0]
    if n_measurements % block_size != 0:
        raise ValueError(
            f"{n_measurements} measurements are not divisible by block size {block_size}"
        )
    n_blocks = n_measurements // block_size
    return matrices.reshape(
        n_blocks,
        block_size,
        *matrices.shape[1:],
    ).mean(axis=1)


def scipy_generalized_eigenvalues(
    correlator: np.ndarray,
    metric: np.ndarray,
) -> np.ndarray:
    values = eigh(
        correlator,
        metric,
        eigvals_only=True,
        check_finite=True,
    )
    return np.real(values[::-1])


def solve_principal_correlators(
    matrices: np.ndarray,
    block_size: int,
    reference_time: int,
    analysis_tmax: int,
    rcond: float,
    production: Any,
    compare_jackknife_reference: bool,
) -> dict[str, Any]:
    blocked = block_ensemble(matrices, block_size)
    central, jackknife = production.central_and_jackknife_matrices(blocked)
    n_blocks = blocked.shape[0]
    n_operators = blocked.shape[-1]
    final_time = min(analysis_tmax, matrices.shape[1] - 1)
    if not 0 <= reference_time < final_time:
        raise ValueError(
            f"reference time {reference_time} must lie in [0, {final_time - 1}]"
        )

    central_lambda = np.full((final_time + 1, n_operators), np.nan)
    jackknife_lambda = np.full(
        (n_blocks, final_time + 1, n_operators),
        np.nan,
    )
    central_lambda[reference_time] = 1.0
    jackknife_lambda[:, reference_time] = 1.0

    metric = central[reference_time]
    metric_diagnostics = production.metric_diagnostics(metric)
    central_metric_status = production.metric_status(metric_diagnostics, rcond)

    maximum_scaled_solver_difference = 0.0
    maximum_absolute_solver_difference = 0.0
    maximum_solver_residual = 0.0
    maximum_metric_orthonormality_residual = 0.0
    central_solve_failures = 0
    jackknife_solve_failures = 0
    jackknife_metric_valid = 0
    identity = np.eye(n_operators)

    if central_metric_status == "ok":
        for time in range(reference_time + 1, final_time + 1):
            try:
                solution = production.solve_gevp(central[time], metric, rcond)
                reference = scipy_generalized_eigenvalues(central[time], metric)
            except (ValueError, np.linalg.LinAlgError):
                central_solve_failures += 1
                continue

            central_lambda[time] = solution.eigenvalues
            difference = np.abs(solution.eigenvalues - reference)
            maximum_absolute_solver_difference = max(
                maximum_absolute_solver_difference,
                float(np.max(difference)),
            )
            scaled = difference / np.maximum(1.0, np.abs(reference))
            maximum_scaled_solver_difference = max(
                maximum_scaled_solver_difference,
                float(np.max(scaled)),
            )
            maximum_solver_residual = max(
                maximum_solver_residual,
                float(solution.max_residual),
            )
            orthonormality = np.linalg.norm(
                solution.eigenvectors.conj().T
                @ metric
                @ solution.eigenvectors
                - identity
            )
            maximum_metric_orthonormality_residual = max(
                maximum_metric_orthonormality_residual,
                float(orthonormality),
            )

    for omitted in range(n_blocks):
        metric_jk = jackknife[omitted, reference_time]
        try:
            diagnostics_jk = production.metric_diagnostics(metric_jk)
            status_jk = production.metric_status(diagnostics_jk, rcond)
        except ValueError:
            status_jk = "invalid"
        if status_jk != "ok":
            continue
        jackknife_metric_valid += 1

        for time in range(reference_time + 1, final_time + 1):
            try:
                solution = production.solve_gevp(
                    jackknife[omitted, time],
                    metric_jk,
                    rcond,
                )
            except ValueError:
                jackknife_solve_failures += 1
                continue
            jackknife_lambda[omitted, time] = solution.eigenvalues

            if compare_jackknife_reference:
                try:
                    reference = scipy_generalized_eigenvalues(
                        jackknife[omitted, time],
                        metric_jk,
                    )
                except np.linalg.LinAlgError:
                    jackknife_solve_failures += 1
                    continue
                difference = np.abs(solution.eigenvalues - reference)
                maximum_absolute_solver_difference = max(
                    maximum_absolute_solver_difference,
                    float(np.max(difference)),
                )
                scaled = difference / np.maximum(1.0, np.abs(reference))
                maximum_scaled_solver_difference = max(
                    maximum_scaled_solver_difference,
                    float(np.max(scaled)),
                )
                maximum_solver_residual = max(
                    maximum_solver_residual,
                    float(solution.max_residual),
                )
                orthonormality = np.linalg.norm(
                    solution.eigenvectors.conj().T
                    @ metric_jk
                    @ solution.eigenvectors
                    - identity
                )
                maximum_metric_orthonormality_residual = max(
                    maximum_metric_orthonormality_residual,
                    float(orthonormality),
                )

    central_effective = production.effective_from_principal(central_lambda)
    jackknife_effective = np.asarray(
        [production.effective_from_principal(sample) for sample in jackknife_lambda]
    )
    valid_effective = np.sum(np.isfinite(jackknife_effective), axis=0)
    effective_error = np.full_like(central_effective, np.nan)
    for time in range(central_effective.shape[0]):
        for state in range(n_operators):
            if valid_effective[time, state] == n_blocks:
                effective_error[time, state] = float(
                    production.jackknife_error(
                        jackknife_effective[:, time, state]
                    )
                )

    return {
        "block_size": block_size,
        "n_blocks": n_blocks,
        "reference_time": reference_time,
        "central_lambda": central_lambda,
        "jackknife_lambda": jackknife_lambda,
        "central_effective": central_effective,
        "jackknife_effective": jackknife_effective,
        "valid_effective": valid_effective,
        "effective_error": effective_error,
        "metric_diagnostics": metric_diagnostics,
        "central_metric_status": central_metric_status,
        "jackknife_metric_valid": jackknife_metric_valid,
        "central_solve_failures": central_solve_failures,
        "jackknife_solve_failures": jackknife_solve_failures,
        "maximum_absolute_solver_difference": maximum_absolute_solver_difference,
        "maximum_scaled_solver_difference": maximum_scaled_solver_difference,
        "maximum_solver_residual": maximum_solver_residual,
        "maximum_metric_orthonormality_residual": (
            maximum_metric_orthonormality_residual
        ),
    }


def correlated_constant_fit_scan(
    analysis: dict[str, Any],
    fit_tmax: int,
    minimum_points: int,
    maximum_points: int,
    production: Any,
) -> list[dict[str, Any]]:
    central = analysis["central_effective"]
    jackknife = analysis["jackknife_effective"]
    reference_time = int(analysis["reference_time"])
    n_blocks = int(analysis["n_blocks"])
    n_states = central.shape[1]
    final_fit_time = min(fit_tmax, central.shape[0] - 1)
    rows: list[dict[str, Any]] = []

    for state in range(n_states):
        for time_min in range(reference_time, final_fit_time + 1):
            final_time_max = min(
                final_fit_time,
                time_min + maximum_points - 1,
            )
            for time_max in range(
                time_min + minimum_points - 1,
                final_time_max + 1,
            ):
                indices = np.arange(time_min, time_max + 1)
                values = central[indices, state]
                samples = jackknife[:, indices, state]
                if not np.all(np.isfinite(values)):
                    continue
                if not np.all(np.isfinite(samples)):
                    continue

                sample_mean = samples.mean(axis=0)
                deviations = samples - sample_mean
                covariance = (
                    (n_blocks - 1) / n_blocks
                    * deviations.T
                    @ deviations
                )
                covariance = 0.5 * (covariance + covariance.T)
                covariance_values, covariance_vectors = np.linalg.eigh(covariance)
                largest = float(covariance_values[-1])
                if not np.isfinite(largest) or largest <= 0.0:
                    continue
                retained = covariance_values > 1.0e-10 * largest
                rank = int(np.count_nonzero(retained))
                if rank < 2:
                    continue
                inverse = (
                    covariance_vectors[:, retained]
                    * (1.0 / covariance_values[retained])[None, :]
                ) @ covariance_vectors[:, retained].T
                ones = np.ones(indices.size)
                denominator = float(ones @ inverse @ ones)
                if not np.isfinite(denominator) or denominator <= 0.0:
                    continue
                weights = inverse @ ones / denominator
                estimate = float(weights @ values)
                jackknife_estimates = samples @ weights
                error = float(production.jackknife_error(jackknife_estimates))
                residual = values - estimate
                chi_square = float(residual @ inverse @ residual)
                dof = rank - 1
                p_value = float(chi2.sf(chi_square, dof))

                rows.append(
                    {
                        "state": state + 1,
                        "time_min": time_min,
                        "time_max": time_max,
                        "n_points": int(indices.size),
                        "covariance_rank": rank,
                        "energy": estimate,
                        "error": error,
                        "chi2": chi_square,
                        "dof": dof,
                        "chi2_per_dof": chi_square / dof,
                        "p_value": p_value,
                        "accepted_p_ge_0p05": int(p_value >= 0.05),
                    }
                )
    return rows


def effective_rows(analyses: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for analysis in analyses:
        central = analysis["central_effective"]
        errors = analysis["effective_error"]
        valid = analysis["valid_effective"]
        lambdas = analysis["central_lambda"]
        n_blocks = int(analysis["n_blocks"])
        reference_time = int(analysis["reference_time"])
        for time in range(reference_time, central.shape[0]):
            for state in range(central.shape[1]):
                rows.append(
                    {
                        "block_size": analysis["block_size"],
                        "n_blocks": n_blocks,
                        "reference_time": reference_time,
                        "time": time,
                        "state": state + 1,
                        "principal_t": lambdas[time, state],
                        "principal_t_plus_1": lambdas[time + 1, state],
                        "effective_energy": central[time, state],
                        "jackknife_error": errors[time, state],
                        "jackknife_valid": int(valid[time, state]),
                        "jackknife_total": n_blocks,
                        "complete": int(valid[time, state] == n_blocks),
                    }
                )
    return rows


def conditioning_rows(analyses: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for analysis in analyses:
        diagnostics = analysis["metric_diagnostics"]
        row: dict[str, Any] = {
            "block_size": analysis["block_size"],
            "n_blocks": analysis["n_blocks"],
            "reference_time": analysis["reference_time"],
            "diagonal_min": diagnostics.diagonal_min,
            "diagonal_max": diagnostics.diagonal_max,
            "condition": diagnostics.condition,
            "central_metric_status": analysis["central_metric_status"],
            "jackknife_metric_valid": analysis["jackknife_metric_valid"],
            "jackknife_total": analysis["n_blocks"],
            "central_solve_failures": analysis["central_solve_failures"],
            "jackknife_solve_failures": analysis["jackknife_solve_failures"],
        }
        for index, value in enumerate(diagnostics.normalized_eigenvalues, start=1):
            row[f"normalized_metric_eigenvalue_{index}"] = value
        rows.append(row)
    return rows


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 180,
            "font.size": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def robust_y_limits(values: np.ndarray, errors: np.ndarray) -> tuple[float, float] | None:
    finite_values = values[np.isfinite(values)]
    if not finite_values.size:
        return None
    lower, upper = np.quantile(finite_values, [0.02, 0.98])
    finite_errors = errors[np.isfinite(errors)]
    typical_error = float(np.median(finite_errors)) if finite_errors.size else 0.0
    padding = max(0.15 * float(upper - lower), 4.0 * typical_error, 0.01)
    return float(lower - padding), float(upper + padding)


def plot_autocorrelations(rows: list[dict[str, Any]], output: Path) -> None:
    lags = np.array([row["lag"] for row in rows])
    lower = np.array([row["rho_q05"] for row in rows])
    median = np.array([row["rho_median"] for row in rows])
    upper = np.array([row["rho_q95"] for row in rows])
    fig, axis = plt.subplots(figsize=(7.0, 4.2))
    axis.fill_between(lags, lower, upper, alpha=0.25, label="5--95% range")
    axis.plot(lags, median, marker="o", markersize=3, label="median")
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("measurement lag")
    axis.set_ylabel(r"primary-stream $\rho(t)$")
    axis.set_title("Autocorrelation of the supplied fake-data streams")
    axis.legend()
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def plot_effective_energies(
    analysis: dict[str, Any],
    plot_tmax: int,
    output: Path,
) -> None:
    central = analysis["central_effective"]
    errors = analysis["effective_error"]
    valid = analysis["valid_effective"]
    total = int(analysis["n_blocks"])
    reference_time = int(analysis["reference_time"])
    final_time = min(plot_tmax, central.shape[0] - 1)
    times = np.arange(reference_time, final_time + 1)

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    for state, axis in enumerate(axes.flat):
        values = central[times, state]
        state_errors = errors[times, state]
        state_valid = valid[times, state]
        complete = (
            np.isfinite(values)
            & np.isfinite(state_errors)
            & (state_valid == total)
        )
        diagnostic = np.isfinite(values) & ~complete
        if np.any(complete):
            axis.errorbar(
                times[complete],
                values[complete],
                yerr=state_errors[complete],
                fmt="o",
                markersize=4,
                capsize=2,
                label="complete block jackknife",
            )
        if np.any(diagnostic):
            axis.scatter(
                times[diagnostic],
                values[diagnostic],
                facecolors="none",
                edgecolors="tab:orange",
                label="diagnostic only",
            )
        axis.set_title(f"state {state + 1}")
        axis.set_ylabel(r"$aE_n^{\mathrm{eff}}(t)$")
        preferred = times <= min(final_time, 20)
        limits = robust_y_limits(values[preferred], state_errors[preferred])
        if limits is not None:
            axis.set_ylim(*limits)
    for axis in axes[-1]:
        axis.set_xlabel("time t")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.947),
            ncol=2,
        )
    fig.suptitle(
        "Fake-data principal effective energies\n"
        f"t_ref={reference_time}, block={analysis['block_size']}, "
        f"blocks={total}",
        y=0.993,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.895))
    fig.subplots_adjust(hspace=0.27, wspace=0.25)
    fig.savefig(output)
    plt.close(fig)


def plot_reference_time_scan(
    analyses: list[dict[str, Any]],
    primary_block_size: int,
    plot_tmax: int,
    output: Path,
) -> None:
    selected = sorted(
        (
            analysis
            for analysis in analyses
            if analysis["block_size"] == primary_block_size
        ),
        key=lambda item: item["reference_time"],
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(selected)))
    for state, axis in enumerate(axes.flat):
        displayed_values = []
        for color, analysis in zip(colors, selected):
            reference_time = int(analysis["reference_time"])
            final_time = min(
                plot_tmax,
                20,
                analysis["central_effective"].shape[0] - 1,
            )
            times = np.arange(reference_time, final_time + 1)
            values = analysis["central_effective"][times, state]
            finite = np.isfinite(values)
            displayed_values.extend(values[finite].tolist())
            axis.plot(
                times[finite],
                values[finite],
                marker="o",
                markersize=2.5,
                linewidth=1.0,
                color=color,
                label=f"t_ref={reference_time}",
            )
        axis.set_title(f"state {state + 1}")
        axis.set_ylabel(r"$aE_n^{\mathrm{eff}}(t)$")
        if displayed_values:
            limits = robust_y_limits(
                np.asarray(displayed_values),
                np.zeros(len(displayed_values)),
            )
            if limits is not None:
                axis.set_ylim(*limits)
    for axis in axes[-1]:
        axis.set_xlabel("time t")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.947),
        ncol=min(4, len(labels)),
    )
    fig.suptitle(
        f"Reference-time scan (block size {primary_block_size})",
        y=0.993,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.895))
    fig.subplots_adjust(hspace=0.27, wspace=0.25)
    fig.savefig(output)
    plt.close(fig)


def plot_block_error_ratios(
    analyses: list[dict[str, Any]],
    primary_reference_time: int,
    output: Path,
) -> None:
    selected = sorted(
        (
            analysis
            for analysis in analyses
            if analysis["reference_time"] == primary_reference_time
        ),
        key=lambda item: item["block_size"],
    )
    by_block = {int(item["block_size"]): item for item in selected}
    if 1 not in by_block:
        return
    baseline = by_block[1]["effective_error"]
    fig, axis = plt.subplots(figsize=(7.2, 4.5))
    for state in range(baseline.shape[1]):
        ratios = []
        block_sizes = []
        for block_size, analysis in sorted(by_block.items()):
            current = analysis["effective_error"]
            first_time = primary_reference_time
            final_time = min(20, current.shape[0] - 1)
            base_values = baseline[first_time : final_time + 1, state]
            current_values = current[first_time : final_time + 1, state]
            valid = (
                np.isfinite(base_values)
                & np.isfinite(current_values)
                & (base_values > 0.0)
            )
            if not np.any(valid):
                continue
            ratios.append(float(np.median(current_values[valid] / base_values[valid])))
            block_sizes.append(block_size)
        axis.plot(
            block_sizes,
            ratios,
            marker="o",
            label=f"state {state + 1}",
        )
    axis.axhline(1.0, color="black", linewidth=0.8)
    axis.set_xlabel("contiguous block size")
    axis.set_ylabel("median error / unblocked error")
    axis.set_title("Block-jackknife error stability")
    axis.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def plot_metric_conditioning(
    analyses: list[dict[str, Any]],
    primary_block_size: int,
    output: Path,
) -> None:
    selected = sorted(
        (
            analysis
            for analysis in analyses
            if analysis["block_size"] == primary_block_size
        ),
        key=lambda item: item["reference_time"],
    )
    times = [item["reference_time"] for item in selected]
    conditions = [item["metric_diagnostics"].condition for item in selected]
    fig, axis = plt.subplots(figsize=(7.0, 4.2))
    axis.semilogy(times, conditions, marker="o")
    axis.set_xlabel("reference time")
    axis.set_ylabel("unit-diagonal metric condition number")
    axis.set_title(f"GEVP metric conditioning (block size {primary_block_size})")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def longest_diagnostic_candidates(
    fit_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    states = sorted({int(row["state"]) for row in fit_rows})
    for state in states:
        accepted = [
            row
            for row in fit_rows
            if int(row["state"]) == state
            and int(row["accepted_p_ge_0p05"]) == 1
        ]
        if not accepted:
            continue
        accepted.sort(
            key=lambda row: (
                -int(row["n_points"]),
                int(row["time_min"]),
                -float(row["p_value"]),
            )
        )
        candidates.append(accepted[0])
    return candidates


def write_report(
    path: Path,
    input_path: Path,
    input_hash: str,
    metadata: dict[str, Any],
    autocorrelation: dict[str, Any],
    analyses: list[dict[str, Any]],
    fit_rows: list[dict[str, Any]],
    primary_block_size: int,
    primary_reference_time: int,
    symmetry_summary: dict[str, float],
    numerical_gates: dict[str, Any],
    output_files: list[str],
    production: Any,
) -> None:
    primary = next(
        analysis
        for analysis in analyses
        if analysis["block_size"] == primary_block_size
        and analysis["reference_time"] == primary_reference_time
    )
    candidates = longest_diagnostic_candidates(fit_rows)
    lines = [
        "FAKE-DATA GEVP VALIDATION",
        "=========================",
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat()}",
        f"Input: {input_path}",
        f"SHA-256: {input_hash}",
        f"Production solver: {Path(production.__file__).resolve()}",
        f"Python: {platform.python_version()}",
        f"NumPy: {np.__version__}",
        f"SciPy: {scipy.__version__}",
        "",
        "Input contract",
        "--------------",
        f"Cell shape: {tuple(metadata['cell_shape'])}",
        f"Element shape: {tuple(metadata['cell_element_shape'])}",
        f"Operators: {metadata['n_operators']}",
        f"Measurements: {metadata['n_measurements']}",
        f"Times: 0..{metadata['n_times'] - 1}",
        f"Element dtypes: {', '.join(metadata['element_dtypes'])}",
        (
            "Maximum input imaginary part: "
            f"{metadata['maximum_input_imaginary_part']:.6e}"
        ),
        "Python time index t corresponds to MATLAB column t+1.",
        "",
        "Hermitian symmetrization",
        "------------------------",
        (
            "Maximum sample-mean anti-Hermitian residual over the analysis range: "
            f"{symmetry_summary['maximum_raw_mean_residual']:.6e}"
        ),
        (
            "Maximum residual after symmetrization: "
            f"{symmetry_summary['maximum_symmetrized_mean_residual']:.6e}"
        ),
        "The analysis uses (C + C^dagger)/2 for every measurement and time.",
        "",
        "Autocorrelation audit",
        "---------------------",
        f"Primary streams: {autocorrelation['n_primary_streams']}",
        f"Fixed tau-int window: {autocorrelation['tau_window']}",
        (
            "tau_int 5/50/95%: "
            f"{autocorrelation['tau_int_q05']:.4f} / "
            f"{autocorrelation['tau_int_median']:.4f} / "
            f"{autocorrelation['tau_int_q95']:.4f}"
        ),
        (
            "Effective sample size from median tau_int: "
            f"{autocorrelation['effective_sample_size_from_median_tau']:.1f}"
        ),
        "Primary statistical analysis: contiguous delete-one-block jackknife.",
        f"Primary block size: {primary_block_size}",
        f"Primary number of blocks: {primary['n_blocks']}",
        "",
        "Numerical validation gates",
        "--------------------------",
    ]
    for name, value in numerical_gates.items():
        rendered = f"{value:.6e}" if isinstance(value, float) else str(value)
        lines.append(f"{name}: {rendered}")

    lines.extend(
        [
            "",
            "Primary metric",
            "--------------",
            f"Reference time: {primary_reference_time}",
            f"Status: {primary['central_metric_status']}",
            f"Condition number: {primary['metric_diagnostics'].condition:.6e}",
            (
                "Normalized eigenvalues: "
                + np.array2string(
                    primary["metric_diagnostics"].normalized_eigenvalues,
                    precision=8,
                )
            ),
            (
                "Valid jackknife metrics: "
                f"{primary['jackknife_metric_valid']}/{primary['n_blocks']}"
            ),
            "",
            "Diagnostic fit-window scan",
            "--------------------------",
            (
                "The rows below are the longest p>=0.05 windows found by the "
                "mechanical scan. They are diagnostics, not final plateau choices, "
                "and no expected energy was supplied to the scan."
            ),
        ]
    )
    if candidates:
        for row in candidates:
            lines.append(
                "state {state}: t={time_min}..{time_max}, "
                "E={energy:.8f} +/- {error:.8f}, "
                "chi2/dof={chi2_per_dof:.3f}, p={p_value:.3f}".format(**row)
            )
    else:
        lines.append("No p>=0.05 candidate window was found.")

    lines.extend(
        [
            "",
            "Interpretation rule",
            "-------------------",
            "Do not hard-code or infer the generator energies from this report.",
            (
                "Final energy levels require inspection of time-, reference-time-, "
                "block-size-, and fit-window stability, followed by confirmation "
                "against the hidden generator values."
            ),
            "",
            "Written files",
            "-------------",
        ]
    )
    lines.extend(output_files)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def self_test(production: Any) -> None:
    energies = np.array([0.10, 0.20, 0.30, 0.40])
    overlaps = np.array(
        [
            [1.0, 0.2, 0.3, 0.1],
            [0.4, 1.1, 0.2, 0.3],
            [0.2, 0.3, 0.9, 0.4],
            [0.1, 0.2, 0.4, 1.2],
        ],
        dtype=np.complex128,
    )
    correlators = np.asarray(
        [
            overlaps
            @ np.diag(np.exp(-energies * time))
            @ overlaps.conj().T
            for time in range(12)
        ]
    )
    principal = np.full((12, 4), np.nan)
    principal[0] = 1.0
    for time in range(1, 12):
        solution = production.solve_gevp(correlators[time], correlators[0], 1e-12)
        expected = np.exp(-energies * time)
        if not np.allclose(solution.eigenvalues, expected, rtol=1e-11, atol=1e-12):
            raise AssertionError("production GEVP failed the analytic spectrum test")
        reference = scipy_generalized_eigenvalues(correlators[time], correlators[0])
        if not np.allclose(solution.eigenvalues, reference, rtol=1e-12, atol=1e-12):
            raise AssertionError("production and SciPy generalized eigenvalues disagree")
        principal[time] = solution.eigenvalues
    effective = production.effective_from_principal(principal)
    if not np.allclose(effective, energies[None, :], rtol=1e-10, atol=1e-11):
        raise AssertionError("effective-energy construction failed")

    replicated = np.repeat(correlators[None, ...], 20, axis=0)
    central, jackknife = production.central_and_jackknife_matrices(replicated)
    if not np.allclose(central, correlators):
        raise AssertionError("central mean helper failed")
    if not np.allclose(jackknife, correlators[None, ...]):
        raise AssertionError("jackknife mean helper failed")
    print("PASS: fake-data validation driver self-test")


def validate_arguments(args: argparse.Namespace, n_times: int) -> None:
    if args.rcond <= 0.0:
        raise ValueError("rcond must be positive")
    if args.primary_block_size not in args.block_sizes:
        raise ValueError("primary block size must also appear in --block-sizes")
    if args.primary_reference_time not in args.reference_times:
        raise ValueError(
            "primary reference time must also appear in --reference-times"
        )
    if args.analysis_tmax >= n_times:
        raise ValueError(
            f"analysis-tmax must be smaller than the number of times ({n_times})"
        )
    if args.plot_tmax > args.analysis_tmax - 1:
        raise ValueError("plot-tmax cannot exceed analysis-tmax-1")
    if args.fit_tmax > args.analysis_tmax - 1:
        raise ValueError("fit-tmax cannot exceed analysis-tmax-1")
    if args.minimum_fit_points < 2:
        raise ValueError("minimum-fit-points must be at least two")
    if args.maximum_fit_points < args.minimum_fit_points:
        raise ValueError("maximum-fit-points must be >= minimum-fit-points")
    if args.tau_window > args.acf_max_lag:
        raise ValueError("tau-window cannot exceed acf-max-lag")


def main() -> None:
    args = parse_args()
    production = production_module()
    if args.self_test:
        self_test(production)
        return

    input_path = args.input.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    input_hash = file_sha256(input_path)
    if not args.skip_checksum and input_hash != args.expected_sha256.lower():
        raise ValueError(
            f"input SHA-256 mismatch: expected {args.expected_sha256}, found {input_hash}"
        )

    raw_matrices, metadata = load_fake_data(input_path)
    validate_arguments(args, int(metadata["n_times"]))
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    self_test(production)
    final_audit_time = min(args.analysis_tmax, raw_matrices.shape[1] - 1)
    raw_mean = raw_matrices.mean(axis=0)
    symmetrized = production.hermitian(raw_matrices)
    symmetrized_mean = symmetrized.mean(axis=0)
    raw_residuals = [
        production.antihermitian_residual(raw_mean[time])
        for time in range(final_audit_time + 1)
    ]
    symmetrized_residuals = [
        production.antihermitian_residual(symmetrized_mean[time])
        for time in range(final_audit_time + 1)
    ]
    symmetry_summary = {
        "maximum_raw_mean_residual": float(max(raw_residuals)),
        "maximum_symmetrized_mean_residual": float(max(symmetrized_residuals)),
    }

    acf_rows, acf_summary = autocorrelation_summary(
        raw_matrices,
        args.acf_max_lag,
        args.tau_window,
    )

    combinations: list[tuple[int, int]] = []
    for reference_time in args.reference_times:
        combinations.append((args.primary_block_size, reference_time))
    for block_size in args.block_sizes:
        combinations.append((block_size, args.primary_reference_time))
    combinations = list(dict.fromkeys(combinations))

    analyses: list[dict[str, Any]] = []
    for block_size, reference_time in combinations:
        print(
            f"Analysing block={block_size}, t_ref={reference_time}",
            flush=True,
        )
        analyses.append(
            solve_principal_correlators(
                symmetrized,
                block_size,
                reference_time,
                args.analysis_tmax,
                args.rcond,
                production,
                compare_jackknife_reference=(
                    block_size == args.primary_block_size
                    and reference_time == args.primary_reference_time
                ),
            )
        )

    primary = next(
        analysis
        for analysis in analyses
        if analysis["block_size"] == args.primary_block_size
        and analysis["reference_time"] == args.primary_reference_time
    )
    fit_rows = correlated_constant_fit_scan(
        primary,
        args.fit_tmax,
        args.minimum_fit_points,
        args.maximum_fit_points,
        production,
    )

    effective_output = output_dir / "fake_gevp_effective.csv"
    conditioning_output = output_dir / "fake_gevp_conditioning.csv"
    acf_output = output_dir / "fake_gevp_autocorrelation.csv"
    fit_output = output_dir / "fake_gevp_fit_scan.csv"
    metadata_output = output_dir / "fake_gevp_metadata.json"
    report_output = output_dir / "fake_gevp_validation.txt"
    effective_plot = output_dir / "fake_gevp_effective_tref5.png"
    reference_plot = output_dir / "fake_gevp_reference_time_scan.png"
    acf_plot = output_dir / "fake_gevp_autocorrelation.png"
    block_plot = output_dir / "fake_gevp_block_error_ratios.png"
    conditioning_plot = output_dir / "fake_gevp_metric_conditioning.png"

    effective_fieldnames = [
        "block_size",
        "n_blocks",
        "reference_time",
        "time",
        "state",
        "principal_t",
        "principal_t_plus_1",
        "effective_energy",
        "jackknife_error",
        "jackknife_valid",
        "jackknife_total",
        "complete",
    ]
    conditioning_fieldnames = [
        "block_size",
        "n_blocks",
        "reference_time",
        "diagonal_min",
        "diagonal_max",
        "normalized_metric_eigenvalue_1",
        "normalized_metric_eigenvalue_2",
        "normalized_metric_eigenvalue_3",
        "normalized_metric_eigenvalue_4",
        "condition",
        "central_metric_status",
        "jackknife_metric_valid",
        "jackknife_total",
        "central_solve_failures",
        "jackknife_solve_failures",
    ]
    acf_fieldnames = ["lag", "rho_q05", "rho_median", "rho_q95"]
    fit_fieldnames = [
        "state",
        "time_min",
        "time_max",
        "n_points",
        "covariance_rank",
        "energy",
        "error",
        "chi2",
        "dof",
        "chi2_per_dof",
        "p_value",
        "accepted_p_ge_0p05",
    ]

    write_csv_file(effective_output, effective_fieldnames, effective_rows(analyses))
    write_csv_file(
        conditioning_output,
        conditioning_fieldnames,
        conditioning_rows(analyses),
    )
    write_csv_file(acf_output, acf_fieldnames, acf_rows)
    write_csv_file(fit_output, fit_fieldnames, fit_rows)

    maximum_absolute_difference = max(
        analysis["maximum_absolute_solver_difference"] for analysis in analyses
    )
    maximum_scaled_difference = max(
        analysis["maximum_scaled_solver_difference"] for analysis in analyses
    )
    maximum_residual = max(
        analysis["maximum_solver_residual"] for analysis in analyses
    )
    maximum_orthonormality = max(
        analysis["maximum_metric_orthonormality_residual"]
        for analysis in analyses
    )
    all_central_metrics_ok = all(
        analysis["central_metric_status"] == "ok" for analysis in analyses
    )
    all_selected_jackknife_metrics_ok = all(
        analysis["jackknife_metric_valid"] == analysis["n_blocks"]
        for analysis in analyses
    )
    all_central_solves_ok = all(
        analysis["central_solve_failures"] == 0 for analysis in analyses
    )
    all_selected_jackknife_solves_ok = all(
        analysis["jackknife_solve_failures"] == 0 for analysis in analyses
    )
    numerical_pass = (
        all_central_metrics_ok
        and all_selected_jackknife_metrics_ok
        and all_central_solves_ok
        and all_selected_jackknife_solves_ok
        and maximum_scaled_difference <= SOLVER_TOLERANCE
        and maximum_residual <= SOLVER_TOLERANCE
        and maximum_orthonormality <= ORTHONORMALITY_TOLERANCE
    )
    numerical_gates = {
        "production_vs_scipy_max_absolute_difference": maximum_absolute_difference,
        "production_vs_scipy_max_scaled_difference": maximum_scaled_difference,
        "maximum_generalized_eigenpair_residual": maximum_residual,
        "maximum_metric_orthonormality_residual": maximum_orthonormality,
        "all_central_metrics_ok": all_central_metrics_ok,
        "all_selected_jackknife_metrics_ok": all_selected_jackknife_metrics_ok,
        "all_central_solves_ok": all_central_solves_ok,
        "all_selected_jackknife_solves_ok": all_selected_jackknife_solves_ok,
        "overall_numerical_validation": "PASS" if numerical_pass else "FAIL",
    }

    metadata_payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(input_path),
        "input_sha256": input_hash,
        "production_solver": str(Path(production.__file__).resolve()),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input": metadata,
        "symmetry": symmetry_summary,
        "autocorrelation": acf_summary,
        "parameters": {
            "rcond": args.rcond,
            "reference_times": list(args.reference_times),
            "primary_reference_time": args.primary_reference_time,
            "block_sizes": list(args.block_sizes),
            "primary_block_size": args.primary_block_size,
            "analysis_tmax": args.analysis_tmax,
            "plot_tmax": args.plot_tmax,
            "fit_tmax": args.fit_tmax,
            "minimum_fit_points": args.minimum_fit_points,
            "maximum_fit_points": args.maximum_fit_points,
        },
        "numerical_gates": numerical_gates,
    }
    metadata_output.write_text(
        json.dumps(metadata_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    configure_plot_style()
    plot_autocorrelations(acf_rows, acf_plot)
    plot_effective_energies(primary, args.plot_tmax, effective_plot)
    plot_reference_time_scan(
        analyses,
        args.primary_block_size,
        args.plot_tmax,
        reference_plot,
    )
    plot_block_error_ratios(
        analyses,
        args.primary_reference_time,
        block_plot,
    )
    plot_metric_conditioning(
        analyses,
        args.primary_block_size,
        conditioning_plot,
    )

    output_files = [
        path.name
        for path in [
            report_output,
            metadata_output,
            effective_output,
            conditioning_output,
            acf_output,
            fit_output,
            effective_plot,
            reference_plot,
            acf_plot,
            block_plot,
            conditioning_plot,
        ]
    ]
    write_report(
        report_output,
        input_path,
        input_hash,
        metadata,
        acf_summary,
        analyses,
        fit_rows,
        args.primary_block_size,
        args.primary_reference_time,
        symmetry_summary,
        numerical_gates,
        output_files,
        production,
    )

    for path in [
        report_output,
        metadata_output,
        effective_output,
        conditioning_output,
        acf_output,
        fit_output,
        effective_plot,
        reference_plot,
        acf_plot,
        block_plot,
        conditioning_plot,
    ]:
        print(f"Wrote: {path}")
    print(f"Overall numerical validation: {'PASS' if numerical_pass else 'FAIL'}")
    if not numerical_pass:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
