#!/usr/bin/env python3
"""Analyse direct Gaussian profiles and full/reduced GEVP bases.

The PROFILE_BASIS=delta run stores the complete Nv x Nv endpoint-mode
kernel.  It can therefore be projected offline onto the three profiles

    sigma = 0.01, 0.05, infinity,

including all mixed source/sink correlators.  This script performs the direct
single-profile comparison, the full 3x3 generalized-eigenvalue analysis, and
the reduced 2x2 analyses that exclude the near-duplicate sigma=0.05/constant
pair.  It does not reread gauge fields or eigenvectors and does not require a
new Slurm job.

The raw scan's T0Fixed=0 is the temporal origin of the Wilson lines.  The
GEVP reference time used below is a different quantity and is consistently
called tau_ref.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# The deterministic self-test exercises only in-memory projection and GEVP
# routines.  Allow it to run on lightweight development hosts where h5py is
# absent; a real analysis still imports and uses h5py through the historical
# delta-analysis module.
if "--self-test" in sys.argv:
    try:
        import h5py as _h5py  # noqa: F401
    except ModuleNotFoundError:
        import types

        sys.modules["h5py"] = types.ModuleType("h5py")

import inspect_laplace_gaussian_from_delta as delta_analysis


PROFILE_SIGMAS = (0.01, 0.05, np.inf)
PROFILE_LABELS = ("sigma=0.0100", "sigma=0.0500", "constant")
PROFILE_DISPLAY = (
    r"$\sigma=0.01$",
    r"$\sigma=0.05$",
    r"$\sigma=\infty\;(\rho_i=1)$",
)
PROFILE_COLORS = ("#0072b2", "#d55e00", "#15284b")
PROFILE_MARKERS = ("s", "^", "o")
N_PROFILES = len(PROFILE_SIGMAS)
REDUCED_BASES = (
    ("sigma001_sigma005", (0, 1)),
    ("sigma001_constant", (0, 2)),
)
REDUCED_BASIS_DISPLAY = {
    "sigma001_sigma005": r"$\{0.01,0.05\}$",
    "sigma001_constant": r"$\{0.01,\infty\}$",
}


@dataclass(frozen=True)
class MetricDiagnostics:
    diagonal_min: float
    diagonal_max: float
    normalized_eigenvalues: np.ndarray
    condition: float


@dataclass(frozen=True)
class GevpSolution:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    metric: MetricDiagnostics
    max_residual: float


@dataclass(frozen=True)
class ScalarResult:
    profile_index: int
    radius: int
    effective: np.ndarray
    effective_error: np.ndarray
    fit: delta_analysis.FitResult


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Project sigma={0.01,0.05,infinity} from a fixed-source-time "
            "delta scan and analyse the full 3x3 and reduced 2x2 GEVPs."
        )
    )
    parser.add_argument(
        "delta_dir",
        type=Path,
        nargs="?",
        help="directory containing cfg<jobid>_<cfg>.out files",
    )
    parser.add_argument(
        "jobid",
        nargs="?",
        help="Slurm array job ID used in the output filenames",
    )
    parser.add_argument(
        "--eigenvalue-root",
        type=Path,
        default=Path("mental/runs_Em1p4_Nv10_qcdnew_full"),
        help=(
            "root containing n<cfg>/eigenvalues/*.h5 "
            "(default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("laplace_static/results/weighted_profile_gevp"),
        help="directory for CSV, TXT, and PNG outputs",
    )
    parser.add_argument(
        "--gevp-reference-times",
        type=int,
        nargs="+",
        default=(1, 2, 3),
        metavar="TAU_REF",
        help="GEVP reference times tau_ref/a to scan (default: 1 2 3)",
    )
    parser.add_argument(
        "--plot-reference-time",
        type=int,
        default=2,
        help="tau_ref/a used in the spectrum and optimized plots",
    )
    parser.add_argument(
        "--optimization-time",
        type=int,
        default=3,
        help=(
            "tau_d/a at which fixed optimized vectors are defined "
            "(default: 3)"
        ),
    )
    parser.add_argument(
        "--gevp-rcond",
        type=float,
        default=1.0e-10,
        help=(
            "minimum allowed eigenvalue of the unit-diagonal GEVP "
            "metric relative to its maximum (default: %(default).1e)"
        ),
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run deterministic in-memory tests and exit",
    )
    parser.add_argument(
        "--self-test-artifact-dir",
        type=Path,
        help=(
            "optional directory in which --self-test renders synthetic "
            "CSV, TXT, and PNG outputs"
        ),
    )
    args = parser.parse_args()

    if not args.self_test and (args.delta_dir is None or args.jobid is None):
        parser.error("delta_dir and jobid are required unless --self-test is used")
    return args


def hermitian(matrix: np.ndarray) -> np.ndarray:
    return 0.5 * (matrix + np.swapaxes(matrix.conj(), -1, -2))


def antihermitian_residual(matrix: np.ndarray) -> float:
    denominator = float(np.linalg.norm(matrix))
    if denominator == 0.0:
        return 0.0
    return float(
        np.linalg.norm(matrix - matrix.conj().T) / denominator
    )


def jackknife_error(samples: np.ndarray) -> np.ndarray:
    """Return delete-one jackknife errors along axis zero."""
    n_cfg = samples.shape[0]
    centered = samples - samples.mean(axis=0)
    return np.sqrt(
        (n_cfg - 1) / n_cfg
        * np.sum(np.abs(centered) ** 2, axis=0)
    )


def endpoint_weights(eigenvalues: np.ndarray) -> np.ndarray:
    """Return w=rho^2 for every cfg, profile, time, and Laplace mode."""
    weights = np.ones(
        (
            eigenvalues.shape[0],
            N_PROFILES,
            eigenvalues.shape[1],
            eigenvalues.shape[2],
        ),
        dtype=float,
    )
    for profile, sigma in enumerate(PROFILE_SIGMAS):
        if np.isinf(sigma):
            continue
        scaled = eigenvalues / sigma
        weights[:, profile] = np.exp(-0.5 * scaled * scaled)
    return weights


def project_profile_matrix(
    data: np.ndarray,
    delta_matrix: np.ndarray,
    eigenvalues: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Project the delta kernel onto all 3x3 mixed profile correlators."""
    weights = endpoint_weights(eigenvalues)
    source = weights[:, :, 0, :]
    sink = weights[:, :, 1 : delta_analysis.N_T + 1, :]
    profiles = np.einsum(
        "cpi,ctrij,cqtj->ctrpq",
        source,
        delta_matrix,
        sink,
        optimize=True,
    )

    constant = profiles[..., 2, 2]
    if not np.allclose(constant, data, rtol=5.0e-12, atol=5.0e-14):
        difference = np.abs(constant - data)
        raise ValueError(
            "sigma=infinity projection does not reproduce DATA; "
            f"max_abs={difference.max():.3e}"
        )

    # Cross-check the sigma=0.05 diagonal against the historical projector.
    old_labels, old_profiles = delta_analysis.project_profiles(
        data,
        delta_matrix,
        eigenvalues,
    )
    old_index = old_labels.index("sigma=0.0500")
    if not np.allclose(
        profiles[..., 1, 1],
        old_profiles[..., old_index],
        rtol=5.0e-12,
        atol=5.0e-14,
    ):
        difference = np.abs(
            profiles[..., 1, 1] - old_profiles[..., old_index]
        )
        raise ValueError(
            "sigma=0.05 regression against historical projector failed; "
            f"max_abs={difference.max():.3e}"
        )

    return profiles, weights


def metric_diagnostics(metric: np.ndarray) -> MetricDiagnostics:
    metric = hermitian(metric)
    diagonal = np.real(np.diag(metric))
    if not np.all(np.isfinite(metric)):
        raise ValueError("GEVP metric contains non-finite values")
    if np.any(diagonal <= 0.0):
        raise ValueError(
            "GEVP metric has a non-positive diagonal entry: "
            + np.array2string(diagonal, precision=5)
        )

    scale = 1.0 / np.sqrt(diagonal)
    normalized = hermitian(
        scale[:, None] * metric * scale[None, :]
    )
    values = np.linalg.eigvalsh(normalized)
    if values[-1] <= 0.0:
        condition = np.inf
    elif values[0] <= 0.0:
        condition = np.inf
    else:
        condition = float(values[-1] / values[0])
    return MetricDiagnostics(
        diagonal_min=float(diagonal.min()),
        diagonal_max=float(diagonal.max()),
        normalized_eigenvalues=values,
        condition=condition,
    )


def metric_status(
    diagnostics: MetricDiagnostics,
    rcond: float,
) -> str:
    """Classify a metric using the same acceptance rule as solve_gevp."""
    values = diagnostics.normalized_eigenvalues
    if values[0] <= 0.0:
        return "non_positive"
    if values[0] <= rcond * values[-1]:
        return "below_rcond"
    return "ok"


def solve_gevp(
    correlator: np.ndarray,
    metric: np.ndarray,
    rcond: float,
) -> GevpSolution:
    """Solve A v=lambda B v using unit-diagonal equilibration and whitening."""
    correlator = hermitian(correlator)
    metric = hermitian(metric)
    diagnostics = metric_diagnostics(metric)
    b_values = diagnostics.normalized_eigenvalues
    if b_values[0] <= 0.0 or b_values[0] <= rcond * b_values[-1]:
        raise ValueError(
            "GEVP metric is not positive definite at the specified cutoff; "
            f"normalized_eigenvalues={np.array2string(b_values, precision=5)}, "
            f"rcond={rcond:.3e}"
        )

    diagonal = np.real(np.diag(metric))
    scale = 1.0 / np.sqrt(diagonal)
    a_scaled = hermitian(
        scale[:, None] * correlator * scale[None, :]
    )
    b_scaled = hermitian(
        scale[:, None] * metric * scale[None, :]
    )

    b_values, b_vectors = np.linalg.eigh(b_scaled)
    b_inverse_sqrt = (
        b_vectors
        * (1.0 / np.sqrt(b_values))[None, :]
    ) @ b_vectors.conj().T
    standard = hermitian(
        b_inverse_sqrt @ a_scaled @ b_inverse_sqrt
    )
    values, standard_vectors = np.linalg.eigh(standard)
    order = np.argsort(values)[::-1]
    values = np.real(values[order])
    standard_vectors = standard_vectors[:, order]

    scaled_vectors = b_inverse_sqrt @ standard_vectors
    vectors = scale[:, None] * scaled_vectors

    # Fix arbitrary phases and explicitly normalize in the B metric.
    for state in range(vectors.shape[1]):
        vector = vectors[:, state]
        norm2 = float(np.real(vector.conj() @ metric @ vector))
        if not np.isfinite(norm2) or norm2 <= 0.0:
            raise ValueError("invalid GEVP eigenvector normalization")
        vector = vector / np.sqrt(norm2)
        pivot = int(np.argmax(np.abs(vector)))
        if abs(vector[pivot]) > 0.0:
            vector = vector * np.exp(-1j * np.angle(vector[pivot]))
        vectors[:, state] = vector

    residuals = []
    norm_a = float(np.linalg.norm(correlator))
    norm_b = float(np.linalg.norm(metric))
    for state, value in enumerate(values):
        vector = vectors[:, state]
        numerator = float(
            np.linalg.norm(correlator @ vector - value * metric @ vector)
        )
        denominator = (
            norm_a * float(np.linalg.norm(vector))
            + abs(value) * norm_b * float(np.linalg.norm(vector))
        )
        residuals.append(numerator / max(denominator, 1.0e-300))

    return GevpSolution(
        eigenvalues=values,
        eigenvectors=vectors,
        metric=diagnostics,
        max_residual=max(residuals),
    )


def scalar_statistics(
    correlator: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central and jackknife effective energies for one scalar channel."""
    real = np.real(correlator)
    n_cfg = real.shape[0]
    total = real.sum(axis=0)
    mean = total / n_cfg
    leave_one_out = (total[None, :] - real) / (n_cfg - 1)

    effective = np.full(delta_analysis.N_T - 1, np.nan)
    valid = (mean[:-1] > 0.0) & (mean[1:] > 0.0)
    effective[valid] = np.log(mean[:-1][valid] / mean[1:][valid])

    jk_effective = np.full(
        (n_cfg, delta_analysis.N_T - 1), np.nan
    )
    valid_jk = (
        (leave_one_out[:, :-1] > 0.0)
        & (leave_one_out[:, 1:] > 0.0)
    )
    jk_effective[valid_jk] = np.log(
        leave_one_out[:, :-1][valid_jk]
        / leave_one_out[:, 1:][valid_jk]
    )

    errors = np.full(delta_analysis.N_T - 1, np.nan)
    for time in range(delta_analysis.N_T - 1):
        if np.all(np.isfinite(jk_effective[:, time])):
            errors[time] = float(jackknife_error(jk_effective[:, time]))
    return effective, errors, jk_effective


def analyse_scalars(
    profile_matrix: np.ndarray,
) -> tuple[list[ScalarResult], list[str]]:
    results = []
    failures = []
    for profile_index, label in enumerate(PROFILE_LABELS):
        for radius in range(1, delta_analysis.N_R):
            correlator = profile_matrix[:, :, radius, profile_index, profile_index]
            effective, error, _ = scalar_statistics(correlator)
            try:
                fit = delta_analysis.fit_profile(label, radius, correlator)
            except ValueError as fit_error:
                if label == "constant":
                    raise
                failures.append(
                    f"scalar {label}, R={radius}: {fit_error}"
                )
                fit = delta_analysis.FitResult(
                    profile=label,
                    radius=radius,
                    effective=effective,
                    fit=np.nan,
                    error=np.nan,
                    chi2_dof=np.nan,
                    p_value=np.nan,
                    early_z=np.nan,
                    jk_fit=np.full(profile_matrix.shape[0], np.nan),
                )
            results.append(
                ScalarResult(
                    profile_index=profile_index,
                    radius=radius,
                    effective=effective,
                    effective_error=error,
                    fit=fit,
                )
            )
    delta_analysis.add_constant_comparisons(
        [result.fit for result in results]
    )
    return results, failures


def central_and_jackknife_matrices(
    matrices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the full mean and all delete-one means."""
    n_cfg = matrices.shape[0]
    total = matrices.sum(axis=0)
    central = total / n_cfg
    jackknife = (total[None, ...] - matrices) / (n_cfg - 1)
    return central, jackknife


def validate_reference_times(values: Iterable[int]) -> tuple[int, ...]:
    result = tuple(sorted(set(values)))
    if not result:
        raise ValueError("at least one GEVP reference time is required")
    for value in result:
        if not (1 <= value < delta_analysis.N_T):
            raise ValueError(
                "GEVP reference times must lie in "
                f"[1,{delta_analysis.N_T - 1}], got {value}"
            )
    return result


def principal_correlators(
    matrices: np.ndarray,
    tau_ref: int,
    rcond: float,
) -> tuple[np.ndarray, np.ndarray, list[str], float]:
    """Return central and jackknife principal correlators for one R,tau_ref."""
    central, jackknife = central_and_jackknife_matrices(matrices)
    n_cfg = matrices.shape[0]
    n_operators = matrices.shape[-1]
    central_lambda = np.full(
        (delta_analysis.N_T, n_operators), np.nan
    )
    jk_lambda = np.full(
        (n_cfg, delta_analysis.N_T, n_operators), np.nan
    )
    central_lambda[tau_ref - 1] = 1.0
    jk_lambda[:, tau_ref - 1] = 1.0
    failures = []
    max_residual = 0.0

    metric = central[tau_ref - 1]
    for tau in range(tau_ref + 1, delta_analysis.N_T + 1):
        try:
            solution = solve_gevp(
                central[tau - 1], metric, rcond
            )
        except ValueError as error:
            failures.append(f"central tau={tau}: {error}")
        else:
            central_lambda[tau - 1] = solution.eigenvalues
            max_residual = max(max_residual, solution.max_residual)

    for omitted in range(n_cfg):
        metric_jk = jackknife[omitted, tau_ref - 1]
        for tau in range(tau_ref + 1, delta_analysis.N_T + 1):
            try:
                solution = solve_gevp(
                    jackknife[omitted, tau - 1], metric_jk, rcond
                )
            except ValueError as error:
                failures.append(
                    f"omit={omitted + 1} tau={tau}: {error}"
                )
            else:
                jk_lambda[omitted, tau - 1] = solution.eigenvalues

    return central_lambda, jk_lambda, failures, max_residual


def effective_from_principal(
    principal: np.ndarray,
) -> np.ndarray:
    """Compute log(lambda(t)/lambda(t+1)) along the time axis."""
    result = np.full(
        (principal.shape[0] - 1, principal.shape[1]), np.nan
    )
    valid = (
        np.isfinite(principal[:-1])
        & np.isfinite(principal[1:])
        & (principal[:-1] > 0.0)
        & (principal[1:] > 0.0)
    )
    result[valid] = np.log(
        principal[:-1][valid] / principal[1:][valid]
    )
    return result


def select_profile_basis(
    profile_matrix: np.ndarray,
    profile_indices: tuple[int, ...],
) -> np.ndarray:
    """Select the same ordered profile subset at source and sink."""
    selected = np.take(profile_matrix, profile_indices, axis=-2)
    return np.take(selected, profile_indices, axis=-1)


def analyse_gevp_basis(
    profile_matrix: np.ndarray,
    reference_times: tuple[int, ...],
    rcond: float,
    basis_name: str | None,
    basis_labels: tuple[str, ...],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Analyse one square profile basis at all radii and reference times."""
    n_operators = profile_matrix.shape[-1]
    if profile_matrix.shape[-2] != n_operators:
        raise ValueError("profile basis must be square")
    if len(basis_labels) != n_operators:
        raise ValueError("basis-label count does not match profile basis")

    raw_mean = profile_matrix.mean(axis=0)
    matrices = hermitian(profile_matrix)
    central, jackknife = central_and_jackknife_matrices(matrices)
    conditioning_rows: list[dict[str, object]] = []
    principal_rows: list[dict[str, object]] = []
    all_failures: list[str] = []
    basis_fields: dict[str, object] = {}
    failure_prefix = ""
    if basis_name is not None:
        basis_fields = {
            "basis": basis_name,
            "basis_profiles": "|".join(basis_labels),
            "basis_size": n_operators,
        }
        failure_prefix = f"basis={basis_name}, "

    for radius in range(1, delta_analysis.N_R):
        radius_matrices = matrices[:, :, radius]
        for tau_ref in reference_times:
            metric = central[tau_ref - 1, radius]
            try:
                diagnostics = metric_diagnostics(metric)
            except ValueError as error:
                diagnostics = None
                central_metric_status = f"invalid: {error}"
            else:
                central_metric_status = metric_status(diagnostics, rcond)

            jk_metric_valid = 0
            for omitted in range(matrices.shape[0]):
                try:
                    item = metric_diagnostics(
                        jackknife[omitted, tau_ref - 1, radius]
                    )
                except ValueError:
                    continue
                if (
                    item.normalized_eigenvalues[0] > 0.0
                    and item.normalized_eigenvalues[0]
                    > rcond * item.normalized_eigenvalues[-1]
                ):
                    jk_metric_valid += 1

            normalized_eigenvalues = (
                np.full(n_operators, np.nan)
                if diagnostics is None
                else diagnostics.normalized_eigenvalues
            )
            conditioning_rows.append(
                {
                    **basis_fields,
                    "radius": radius,
                    "tau_ref": tau_ref,
                    "antihermitian_residual": antihermitian_residual(
                        raw_mean[tau_ref - 1, radius]
                    ),
                    "diagonal_min": (
                        np.nan
                        if diagnostics is None
                        else diagnostics.diagonal_min
                    ),
                    "diagonal_max": (
                        np.nan
                        if diagnostics is None
                        else diagnostics.diagonal_max
                    ),
                    "metric_eval_1": normalized_eigenvalues[0],
                    "metric_eval_2": (
                        normalized_eigenvalues[1]
                        if n_operators >= 2
                        else np.nan
                    ),
                    "metric_eval_3": (
                        normalized_eigenvalues[2]
                        if n_operators >= 3
                        else np.nan
                    ),
                    "metric_condition": (
                        np.inf
                        if diagnostics is None
                        else diagnostics.condition
                    ),
                    "central_metric_status": central_metric_status,
                    "central_metric_valid": int(
                        central_metric_status == "ok"
                    ),
                    "jackknife_metric_valid": jk_metric_valid,
                    "jackknife_total": matrices.shape[0],
                }
            )

            central_lambda, jk_lambda, failures, max_residual = (
                principal_correlators(
                    radius_matrices, tau_ref, rcond
                )
            )
            all_failures.extend(
                f"{failure_prefix}R={radius}, tau_ref={tau_ref}, {failure}"
                for failure in failures
            )
            central_effective = effective_from_principal(central_lambda)
            jk_effective = np.stack(
                [
                    effective_from_principal(sample)
                    for sample in jk_lambda
                ]
            )

            for tau in range(tau_ref + 1, delta_analysis.N_T):
                for state in range(n_operators):
                    samples = jk_effective[:, tau - 1, state]
                    valid_count = int(np.count_nonzero(np.isfinite(samples)))
                    error = (
                        float(jackknife_error(samples))
                        if valid_count == matrices.shape[0]
                        else np.nan
                    )
                    principal_rows.append(
                        {
                            **basis_fields,
                            "radius": radius,
                            "tau_ref": tau_ref,
                            "effective_time": tau,
                            "state": state + 1,
                            "principal_effective_energy": central_effective[
                                tau - 1, state
                            ],
                            "jackknife_error": error,
                            "jackknife_valid": valid_count,
                            "jackknife_total": matrices.shape[0],
                            "reference_time_condition_met": int(
                                tau_ref >= (tau + 1) / 2.0
                            ),
                            "max_central_gevp_residual": max_residual,
                        }
                    )

    return conditioning_rows, principal_rows, all_failures


def analyse_gevp(
    profile_matrix: np.ndarray,
    reference_times: tuple[int, ...],
    rcond: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Analyse the full three-profile basis."""
    return analyse_gevp_basis(
        profile_matrix,
        reference_times,
        rcond,
        None,
        PROFILE_LABELS,
    )


def analyse_reduced_gevps(
    profile_matrix: np.ndarray,
    reference_times: tuple[int, ...],
    rcond: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Analyse the two reduced bases containing the narrow profile."""
    conditioning_rows: list[dict[str, object]] = []
    principal_rows: list[dict[str, object]] = []
    failures: list[str] = []
    for basis_name, profile_indices in REDUCED_BASES:
        labels = tuple(PROFILE_LABELS[index] for index in profile_indices)
        selected = select_profile_basis(profile_matrix, profile_indices)
        basis_conditioning, basis_principal, basis_failures = (
            analyse_gevp_basis(
                selected,
                reference_times,
                rcond,
                basis_name,
                labels,
            )
        )
        conditioning_rows.extend(basis_conditioning)
        principal_rows.extend(basis_principal)
        failures.extend(basis_failures)
    return conditioning_rows, principal_rows, failures


def analyse_fixed_vectors(
    profile_matrix: np.ndarray,
    tau_ref: int,
    tau_d: int,
    rcond: float,
    basis_name: str | None = None,
    basis_labels: tuple[str, ...] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Analyse V^dagger C(tau) V using V from (tau_d,tau_ref)."""
    if tau_d <= tau_ref:
        raise ValueError("optimization time must be larger than tau_ref")
    n_operators = profile_matrix.shape[-1]
    if profile_matrix.shape[-2] != n_operators:
        raise ValueError("profile basis must be square")
    if basis_labels is None:
        basis_labels = PROFILE_LABELS
    if len(basis_labels) != n_operators:
        raise ValueError("basis-label count does not match profile basis")
    basis_fields: dict[str, object] = {}
    failure_prefix = ""
    if basis_name is not None:
        basis_fields = {
            "basis": basis_name,
            "basis_profiles": "|".join(basis_labels),
            "basis_size": n_operators,
        }
        failure_prefix = f"basis={basis_name}, "
    matrices = hermitian(profile_matrix)
    n_cfg = matrices.shape[0]
    central, jackknife = central_and_jackknife_matrices(matrices)
    projected_rows: list[dict[str, object]] = []
    vector_rows: list[dict[str, object]] = []
    failures: list[str] = []

    for radius in range(1, delta_analysis.N_R):
        try:
            solution = solve_gevp(
                central[tau_d - 1, radius],
                central[tau_ref - 1, radius],
                rcond,
            )
        except ValueError as error:
            failures.append(
                f"{failure_prefix}R={radius}, central optimized vectors: "
                f"{error}"
            )
            continue

        vectors = solution.eigenvectors
        for state in range(n_operators):
            for profile in range(n_operators):
                coefficient = vectors[profile, state]
                vector_rows.append(
                    {
                        **basis_fields,
                        "radius": radius,
                        "tau_ref": tau_ref,
                        "tau_d": tau_d,
                        "state": state + 1,
                        "profile": basis_labels[profile],
                        "coefficient_real": float(np.real(coefficient)),
                        "coefficient_imag": float(np.imag(coefficient)),
                    }
                )

        projected = np.empty(
            (delta_analysis.N_T, n_operators, n_operators),
            dtype=complex,
        )
        for tau in range(1, delta_analysis.N_T + 1):
            projected[tau - 1] = (
                vectors.conj().T
                @ central[tau - 1, radius]
                @ vectors
            )

        jk_diagonal = np.full(
            (n_cfg, delta_analysis.N_T, n_operators), np.nan
        )
        for omitted in range(n_cfg):
            try:
                jk_solution = solve_gevp(
                    jackknife[omitted, tau_d - 1, radius],
                    jackknife[omitted, tau_ref - 1, radius],
                    rcond,
                )
            except ValueError as error:
                failures.append(
                    f"{failure_prefix}R={radius}, omit={omitted + 1}, "
                    f"optimized vectors: {error}"
                )
                continue
            jk_vectors = jk_solution.eigenvectors
            for tau in range(1, delta_analysis.N_T + 1):
                block = (
                    jk_vectors.conj().T
                    @ jackknife[omitted, tau - 1, radius]
                    @ jk_vectors
                )
                jk_diagonal[omitted, tau - 1] = np.real(np.diag(block))

        central_diagonal = np.real(
            np.diagonal(projected, axis1=-2, axis2=-1)
        )
        central_effective = np.full(
            (delta_analysis.N_T - 1, n_operators), np.nan
        )
        valid = (
            (central_diagonal[:-1] > 0.0)
            & (central_diagonal[1:] > 0.0)
        )
        central_effective[valid] = np.log(
            central_diagonal[:-1][valid]
            / central_diagonal[1:][valid]
        )

        jk_effective = np.full(
            (n_cfg, delta_analysis.N_T - 1, n_operators), np.nan
        )
        valid_jk = (
            np.isfinite(jk_diagonal[:, :-1])
            & np.isfinite(jk_diagonal[:, 1:])
            & (jk_diagonal[:, :-1] > 0.0)
            & (jk_diagonal[:, 1:] > 0.0)
        )
        jk_effective[valid_jk] = np.log(
            jk_diagonal[:, :-1][valid_jk]
            / jk_diagonal[:, 1:][valid_jk]
        )

        for tau in range(1, delta_analysis.N_T):
            for state in range(n_operators):
                samples = jk_effective[:, tau - 1, state]
                valid_count = int(np.count_nonzero(np.isfinite(samples)))
                error = (
                    float(jackknife_error(samples))
                    if valid_count == n_cfg
                    else np.nan
                )
                projected_rows.append(
                    {
                        **basis_fields,
                        "radius": radius,
                        "tau_ref": tau_ref,
                        "tau_d": tau_d,
                        "effective_time": tau,
                        "state": state + 1,
                        "projected_effective_energy": central_effective[
                            tau - 1, state
                        ],
                        "jackknife_error": error,
                        "jackknife_valid": valid_count,
                        "jackknife_total": n_cfg,
                    }
                )

    return projected_rows, vector_rows, failures


def analyse_reduced_fixed_vectors(
    profile_matrix: np.ndarray,
    tau_ref: int,
    tau_d: int,
    rcond: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Analyse fixed vectors for both reduced two-profile bases."""
    projected_rows: list[dict[str, object]] = []
    vector_rows: list[dict[str, object]] = []
    failures: list[str] = []
    for basis_name, profile_indices in REDUCED_BASES:
        labels = tuple(PROFILE_LABELS[index] for index in profile_indices)
        selected = select_profile_basis(profile_matrix, profile_indices)
        basis_projected, basis_vectors, basis_failures = analyse_fixed_vectors(
            selected,
            tau_ref,
            tau_d,
            rcond,
            basis_name,
            labels,
        )
        projected_rows.extend(basis_projected)
        vector_rows.extend(basis_vectors)
        failures.extend(basis_failures)
    return projected_rows, vector_rows, failures


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    if fieldnames is None:
        if not rows:
            raise ValueError(
                f"fieldnames are required to write empty CSV {path}"
            )
        fieldnames = list(rows[0])
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def scalar_rows(results: list[ScalarResult]) -> list[dict[str, object]]:
    rows = []
    for result in results:
        for time in range(1, delta_analysis.N_T):
            rows.append(
                {
                    "profile": PROFILE_LABELS[result.profile_index],
                    "sigma": (
                        "inf"
                        if np.isinf(PROFILE_SIGMAS[result.profile_index])
                        else f"{PROFILE_SIGMAS[result.profile_index]:.4f}"
                    ),
                    "radius": result.radius,
                    "effective_time": time,
                    "effective_energy": result.effective[time - 1],
                    "jackknife_error": result.effective_error[time - 1],
                    "plateau_fit_T3_5": result.fit.fit,
                    "plateau_fit_error": result.fit.error,
                    "chi2_dof": result.fit.chi2_dof,
                    "p_value": result.fit.p_value,
                    "z2": result.fit.early_z,
                    "shift_from_constant": result.fit.delta_constant,
                    "paired_shift_z": result.fit.delta_constant_z,
                }
            )
    return rows


def index_scalar_results(
    results: list[ScalarResult],
) -> dict[tuple[int, int], ScalarResult]:
    return {
        (result.profile_index, result.radius): result
        for result in results
    }


def index_rows(
    rows: list[dict[str, object]],
    keys: tuple[str, ...],
) -> dict[tuple[object, ...], dict[str, object]]:
    return {tuple(row[key] for key in keys): row for row in rows}


def row_validity_masks(
    rows: list[dict[str, object]],
    value_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return sorted series arrays and complete/incomplete validity masks."""
    selected = sorted(rows, key=lambda row: int(row["effective_time"]))
    x = np.array([row["effective_time"] for row in selected], dtype=float)
    y = np.array([row[value_key] for row in selected], dtype=float)
    error = np.array([row["jackknife_error"] for row in selected], dtype=float)
    complete = np.array(
        [
            row["jackknife_valid"] == row["jackknife_total"]
            for row in selected
        ],
        dtype=bool,
    )
    complete &= np.isfinite(y) & np.isfinite(error)
    diagnostic = np.isfinite(y) & ~complete
    return x, y, error, complete, diagnostic


def plot_scalar_comparison(
    path: Path,
    results: list[ScalarResult],
) -> None:
    indexed = index_scalar_results(results)
    figure, axes = plt.subplots(2, 3, figsize=(14.5, 8.7), sharex=True)
    times = np.arange(1, delta_analysis.N_T)
    for radius, axis in zip(range(1, delta_analysis.N_R), axes.flat):
        axis.axhline(0.0, color="#555555", linewidth=0.9)
        axis.axvspan(3, 5, color="#d9e6f2", alpha=0.45)
        for profile in range(N_PROFILES):
            item = indexed[(profile, radius)]
            if not np.isfinite(item.fit.fit):
                continue
            residual = item.effective - item.fit.fit
            complete = np.isfinite(residual) & np.isfinite(
                item.effective_error
            )
            diagnostic = np.isfinite(residual) & ~complete
            if np.any(complete):
                axis.errorbar(
                    times[complete],
                    residual[complete],
                    yerr=item.effective_error[complete],
                    color=PROFILE_COLORS[profile],
                    marker=PROFILE_MARKERS[profile],
                    markersize=4.5,
                    linewidth=1.2,
                    capsize=2.0,
                    label=PROFILE_DISPLAY[profile],
                )
            if np.any(diagnostic):
                axis.scatter(
                    times[diagnostic],
                    residual[diagnostic],
                    edgecolors=PROFILE_COLORS[profile],
                    facecolors="none",
                    marker=PROFILE_MARKERS[profile],
                    s=30,
                    alpha=0.55,
                    label=(
                        PROFILE_DISPLAY[profile]
                        if not np.any(complete)
                        else None
                    ),
                )
        axis.set_title(rf"$R/a={radius}$")
        axis.grid(alpha=0.22)
    for axis in axes[-1]:
        axis.set_xlabel(r"$\tau/a$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$V_{\rm eff}(\tau)-\widehat V_{3:5}$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=3,
        frameon=False,
    )
    figure.suptitle(
        "Direct profile comparison: three Gaussian widths",
        y=0.975,
        fontsize=15,
    )
    figure.text(
        0.5,
        0.025,
        "Fixed source time 0, x-oriented 100-configuration pilot. "
        "Each profile is centred on its own correlated tau/a=3..5 fit.",
        ha="center",
        fontsize=9.5,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.825,
        wspace=0.17,
        hspace=0.25,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_ground_reference_scan(
    path: Path,
    principal_rows: list[dict[str, object]],
    scalar_results: list[ScalarResult],
    reference_times: tuple[int, ...],
) -> None:
    scalar_index = index_scalar_results(scalar_results)
    colors = ("#0072b2", "#d55e00", "#009e73", "#cc79a7")
    markers = ("o", "s", "^", "D")
    figure, axes = plt.subplots(2, 3, figsize=(14.5, 8.7), sharex=True)

    for radius, axis in zip(range(1, delta_analysis.N_R), axes.flat):
        reference = scalar_index[(2, radius)]
        axis.axhspan(
            reference.fit.fit - reference.fit.error,
            reference.fit.fit + reference.fit.error,
            color="#15284b",
            alpha=0.12,
            label=r"$\rho_i=1$ scalar plateau",
        )
        axis.axhline(reference.fit.fit, color="#15284b", linewidth=0.9)
        for index, tau_ref in enumerate(reference_times):
            selected = [
                row
                for row in principal_rows
                if row["radius"] == radius
                and row["tau_ref"] == tau_ref
                and row["state"] == 1
            ]
            if not selected:
                continue
            x, y, error, complete, diagnostic = row_validity_masks(
                selected, "principal_effective_energy"
            )
            preferred = np.array(
                [
                    bool(row["reference_time_condition_met"])
                    for row in sorted(
                        selected,
                        key=lambda row: int(row["effective_time"]),
                    )
                ]
            )
            label = rf"$\tau_{{\rm ref}}/a={tau_ref}$"
            color = colors[index % len(colors)]
            marker = markers[index % len(markers)]
            if np.any(complete):
                axis.errorbar(
                    x[complete],
                    y[complete],
                    yerr=error[complete],
                    color=color,
                    marker=marker,
                    markersize=4.5,
                    linewidth=1.1,
                    capsize=2.0,
                    alpha=0.45,
                    label=label,
                )
            if np.any(diagnostic):
                axis.scatter(
                    x[diagnostic],
                    y[diagnostic],
                    edgecolors=color,
                    facecolors="none",
                    marker=marker,
                    s=31,
                    alpha=0.55,
                    label=label if not np.any(complete) else None,
                )
            preferred_complete = preferred & complete
            if np.any(preferred_complete):
                axis.scatter(
                    x[preferred_complete],
                    y[preferred_complete],
                    color=color,
                    marker=marker,
                    s=31,
                    zorder=4,
                )
        axis.set_title(rf"$R/a={radius}$")
        axis.grid(alpha=0.22)
    for axis in axes[-1]:
        axis.set_xlabel(r"effective time $\tau/a$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$E^{\rm eff}_1(\tau,\tau_{\rm ref})$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    figure.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=4,
        frameon=False,
    )
    figure.suptitle(
        "Three-profile GEVP ground-state reference-time scan",
        y=0.975,
        fontsize=15,
    )
    figure.text(
        0.5,
        0.025,
        "Filled markers satisfy tau_ref >= (tau+1)/2 for the adjacent-time "
        "effective energy. Hollow markers have incomplete jackknife support.",
        ha="center",
        fontsize=9.5,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.825,
        wspace=0.17,
        hspace=0.25,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_states(
    path: Path,
    rows: list[dict[str, object]],
    value_key: str,
    tau_ref: int,
    scalar_results: list[ScalarResult],
    title: str,
    footer: str,
) -> None:
    scalar_index = index_scalar_results(scalar_results)
    state_colors = ("#0072b2", "#d55e00", "#009e73")
    state_markers = ("o", "s", "^")
    figure, axes = plt.subplots(2, 3, figsize=(14.5, 8.7), sharex=True)
    for radius, axis in zip(range(1, delta_analysis.N_R), axes.flat):
        reference = scalar_index[(2, radius)]
        axis.axhspan(
            reference.fit.fit - reference.fit.error,
            reference.fit.fit + reference.fit.error,
            color="#15284b",
            alpha=0.10,
        )
        axis.axhline(reference.fit.fit, color="#15284b", linewidth=0.8)
        for state in range(1, N_PROFILES + 1):
            selected = [
                row
                for row in rows
                if row["radius"] == radius
                and row["tau_ref"] == tau_ref
                and row["state"] == state
            ]
            if not selected:
                continue
            x, y, error, complete, diagnostic = row_validity_masks(
                selected, value_key
            )
            label = rf"state $n={state}$"
            if np.any(complete):
                axis.errorbar(
                    x[complete],
                    y[complete],
                    yerr=error[complete],
                    color=state_colors[state - 1],
                    marker=state_markers[state - 1],
                    markersize=4.5,
                    linewidth=1.1,
                    capsize=2.0,
                    label=label,
                )
            if np.any(diagnostic):
                axis.scatter(
                    x[diagnostic],
                    y[diagnostic],
                    edgecolors=state_colors[state - 1],
                    facecolors="none",
                    marker=state_markers[state - 1],
                    s=31,
                    alpha=0.55,
                    label=label if not np.any(complete) else None,
                )
        if not any(
            row["radius"] == radius and row["tau_ref"] == tau_ref
            for row in rows
        ):
            axis.text(
                0.5,
                0.5,
                "GEVP unavailable",
                transform=axis.transAxes,
                ha="center",
                va="center",
                color="#777777",
            )
        axis.set_title(rf"$R/a={radius}$")
        axis.grid(alpha=0.22)
    for axis in axes[-1]:
        axis.set_xlabel(r"effective time $\tau/a$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"effective energy")
    handles = []
    labels = []
    for axis in axes.flat:
        new_handles, new_labels = axis.get_legend_handles_labels()
        handles.extend(new_handles)
        labels.extend(new_labels)
    unique = dict(zip(labels, handles))
    figure.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=3,
        frameon=False,
    )
    figure.suptitle(title, y=0.975, fontsize=15)
    figure.text(0.5, 0.025, footer, ha="center", fontsize=9.5)
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.825,
        wspace=0.17,
        hspace=0.25,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def plot_reduced_states(
    path: Path,
    rows: list[dict[str, object]],
    value_key: str,
    tau_ref: int,
    scalar_results: list[ScalarResult],
    title: str,
    footer: str,
) -> None:
    """Plot both states from each reduced basis with validity markers."""
    scalar_index = index_scalar_results(scalar_results)
    basis_colors = {
        "sigma001_sigma005": "#0072b2",
        "sigma001_constant": "#d55e00",
    }
    state_markers = {1: "o", 2: "s"}
    state_linestyles = {1: "-", 2: "--"}
    figure, axes = plt.subplots(2, 3, figsize=(14.5, 8.7), sharex=True)

    for radius, axis in zip(range(1, delta_analysis.N_R), axes.flat):
        reference = scalar_index[(2, radius)]
        axis.axhspan(
            reference.fit.fit - reference.fit.error,
            reference.fit.fit + reference.fit.error,
            color="#15284b",
            alpha=0.10,
        )
        axis.axhline(reference.fit.fit, color="#15284b", linewidth=0.8)
        finite_series = 0
        for basis_name, _ in REDUCED_BASES:
            for state in (1, 2):
                selected = [
                    row
                    for row in rows
                    if row["basis"] == basis_name
                    and row["radius"] == radius
                    and row["tau_ref"] == tau_ref
                    and row["state"] == state
                ]
                if not selected:
                    continue
                x, y, error, complete, diagnostic = row_validity_masks(
                    selected, value_key
                )
                label = (
                    f"{REDUCED_BASIS_DISPLAY[basis_name]}, "
                    rf"$n={state}$"
                )
                color = basis_colors[basis_name]
                marker = state_markers[state]
                if np.any(complete):
                    finite_series += 1
                    axis.errorbar(
                        x[complete],
                        y[complete],
                        yerr=error[complete],
                        color=color,
                        marker=marker,
                        linestyle=state_linestyles[state],
                        markersize=4.5,
                        linewidth=1.1,
                        capsize=2.0,
                        label=label,
                    )
                if np.any(diagnostic):
                    finite_series += 1
                    axis.scatter(
                        x[diagnostic],
                        y[diagnostic],
                        edgecolors=color,
                        facecolors="none",
                        marker=marker,
                        s=31,
                        alpha=0.55,
                        label=label if not np.any(complete) else None,
                    )
        if finite_series == 0:
            axis.text(
                0.5,
                0.5,
                "No finite estimates",
                transform=axis.transAxes,
                ha="center",
                va="center",
                color="#777777",
            )
        axis.set_title(rf"$R/a={radius}$")
        axis.grid(alpha=0.22)

    for axis in axes[-1]:
        axis.set_xlabel(r"effective time $\tau/a$")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"effective energy")
    handles: list[object] = []
    labels: list[str] = []
    for axis in axes.flat:
        new_handles, new_labels = axis.get_legend_handles_labels()
        handles.extend(new_handles)
        labels.extend(new_labels)
    unique = dict(zip(labels, handles))
    figure.legend(
        unique.values(),
        unique.keys(),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=4,
        frameon=False,
    )
    figure.suptitle(title, y=0.975, fontsize=15)
    figure.text(0.5, 0.025, footer, ha="center", fontsize=9.5)
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.825,
        wspace=0.17,
        hspace=0.25,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_report(
    path: Path,
    jobid: str,
    n_cfg: int,
    eigenvalues: np.ndarray,
    weights: np.ndarray,
    scalar_results: list[ScalarResult],
    conditioning_rows: list[dict[str, object]],
    principal_rows: list[dict[str, object]],
    projected_rows: list[dict[str, object]],
    reduced_conditioning_rows: list[dict[str, object]],
    reduced_principal_rows: list[dict[str, object]],
    reduced_projected_rows: list[dict[str, object]],
    failures: list[str],
    reference_times: tuple[int, ...],
    plot_reference_time: int,
    optimization_time: int,
    rcond: float,
) -> None:
    scalar_index = index_scalar_results(scalar_results)
    with path.open("w") as stream:
        print("GAUSSIAN PROFILE AND GEVP ANALYSIS", file=stream)
        print(f"Source job: {jobid}", file=stream)
        print(f"Configurations: {n_cfg}", file=stream)
        print(
            "Profiles: sigma=0.01, sigma=0.05, sigma=infinity (rho_i=1)",
            file=stream,
        )
        print(
            "Scope: fixed Wilson-line source time 0, x-oriented separation, "
            "R/a=1..6, tau/a=1..8",
            file=stream,
        )
        print(
            "GEVP reference times: " + ", ".join(map(str, reference_times)),
            file=stream,
        )
        print(f"GEVP metric rcond: {rcond:.3e}", file=stream)
        print(
            "Note: T0Fixed=0 in the raw scan and tau_ref in the GEVP are "
            "different quantities.\n",
            file=stream,
        )

        print("Eigenvalue and endpoint-weight ranges", file=stream)
        print(
            f"lambda retained: {eigenvalues.min():.9f} .. "
            f"{eigenvalues.max():.9f}",
            file=stream,
        )
        for profile, label in enumerate(PROFILE_LABELS):
            print(
                f"{label:13s} w=rho^2 range: "
                f"{weights[:, profile].min():.6e} .. "
                f"{weights[:, profile].max():.6e}",
                file=stream,
            )

        print("\nDirect scalar profiles: correlated tau/a=3..5 fits", file=stream)
        print(
            "profile        R       fit       err chi2/dof      p      z2 "
            " dEconst  zconst",
            file=stream,
        )
        for profile in range(N_PROFILES):
            for radius in range(1, delta_analysis.N_R):
                fit = scalar_index[(profile, radius)].fit
                print(
                    f"{PROFILE_LABELS[profile]:13s} {radius:d} "
                    f"{fit.fit:+9.6f} {fit.error:9.6f} "
                    f"{fit.chi2_dof:8.3f} {fit.p_value:6.3f} "
                    f"{fit.early_z:+7.2f} {fit.delta_constant:+9.6f} "
                    f"{fit.delta_constant_z:+7.2f}",
                    file=stream,
                )

        print("\nFull 3x3 GEVP metric diagnostics", file=stream)
        print(
            "R tref antiherm  b1(unitdiag) b2(unitdiag) b3(unitdiag) "
            "condition status JKmetric",
            file=stream,
        )
        for row in conditioning_rows:
            print(
                f"{row['radius']} {row['tau_ref']} "
                f"{float(row['antihermitian_residual']):8.2e} "
                f"{float(row['metric_eval_1']):12.4e} "
                f"{float(row['metric_eval_2']):12.4e} "
                f"{float(row['metric_eval_3']):12.4e} "
                f"{float(row['metric_condition']):10.3e} "
                f"{row['central_metric_status']} "
                f"{row['jackknife_metric_valid']}/{row['jackknife_total']}",
                file=stream,
            )

        print("\nFull 3x3 GEVP principal effective energies", file=stream)
        print("R tref tau state Eeff err JKvalid tref_condition", file=stream)
        for row in principal_rows:
            print(
                f"{row['radius']} {row['tau_ref']} "
                f"{row['effective_time']} {row['state']} "
                f"{float(row['principal_effective_energy']):+.8e} "
                f"{float(row['jackknife_error']):.8e} "
                f"{row['jackknife_valid']}/{row['jackknife_total']} "
                f"{row['reference_time_condition_met']}",
                file=stream,
            )

        print("\nFixed-vector optimized correlators", file=stream)
        print(
            f"Vectors: tau_ref/a={plot_reference_time}, "
            f"tau_d/a={optimization_time}",
            file=stream,
        )
        print(
            "C_opt(tau) = V^dagger C(tau) V.  Its n=2 diagonal is the "
            "candidate first-excited fixed-vector correlator. "
            "It equals the corresponding principal correlator at tau_d by "
            "construction, but need not equal lambda_2 at every tau.",
            file=stream,
        )
        print("R tau state Eeff err JKvalid", file=stream)
        for row in projected_rows:
            print(
                f"{row['radius']} {row['effective_time']} {row['state']} "
                f"{float(row['projected_effective_energy']):+.8e} "
                f"{float(row['jackknife_error']):.8e} "
                f"{row['jackknife_valid']}/{row['jackknife_total']}",
                file=stream,
            )

        print("\nReduced 2x2 GEVP metric diagnostics", file=stream)
        print(
            "basis R tref antiherm b1(unitdiag) b2(unitdiag) condition "
            "status JKmetric",
            file=stream,
        )
        for row in reduced_conditioning_rows:
            print(
                f"{row['basis']} {row['radius']} {row['tau_ref']} "
                f"{float(row['antihermitian_residual']):8.2e} "
                f"{float(row['metric_eval_1']):12.4e} "
                f"{float(row['metric_eval_2']):12.4e} "
                f"{float(row['metric_condition']):10.3e} "
                f"{row['central_metric_status']} "
                f"{row['jackknife_metric_valid']}/{row['jackknife_total']}",
                file=stream,
            )

        print(
            "\nReduced 2x2 principal effective energies in the "
            "reference-time regime",
            file=stream,
        )
        print(
            "basis R tref tau state Eeff err JKvalid",
            file=stream,
        )
        for row in reduced_principal_rows:
            if not row["reference_time_condition_met"]:
                continue
            print(
                f"{row['basis']} {row['radius']} {row['tau_ref']} "
                f"{row['effective_time']} {row['state']} "
                f"{float(row['principal_effective_energy']):+.8e} "
                f"{float(row['jackknife_error']):.8e} "
                f"{row['jackknife_valid']}/{row['jackknife_total']}",
                file=stream,
            )

        print(
            "\nReduced 2x2 candidate first-excited fixed-vector "
            "correlators (tau/a=2..5)",
            file=stream,
        )
        print("basis R tau Eeff err JKvalid", file=stream)
        for row in reduced_projected_rows:
            if row["state"] != 2 or not (2 <= row["effective_time"] <= 5):
                continue
            print(
                f"{row['basis']} {row['radius']} "
                f"{row['effective_time']} "
                f"{float(row['projected_effective_energy']):+.8e} "
                f"{float(row['jackknife_error']):.8e} "
                f"{row['jackknife_valid']}/{row['jackknife_total']}",
                file=stream,
            )

        print("\nInterpretation guardrails", file=stream)
        print(
            "- A usable GEVP requires a positive, well-conditioned "
            "C(tau_ref) in the full sample and every delete-one sample.",
            file=stream,
        )
        print(
            "- A positive lambda_n and a finite all-sample jackknife error "
            "are required before quoting E_n.",
            file=stream,
        )
        print(
            "- Stability across tau_ref, basis choice, and an effective-energy "
            "plateau is required before interpreting n=2 or n=3 as an "
            "excited static potential.",
            file=stream,
        )
        print(
            "- The condition tau_ref >= (tau+1)/2 marks the regime advocated "
            "in arXiv:0902.1265 for the adjacent-time effective energy.",
            file=stream,
        )
        print(
            "- The 192x64^3 ensemble and the full 3D clover profile are not "
            "included in this calculation.",
            file=stream,
        )

        if failures:
            print(
                f"\nDiagnostic exclusions ({len(failures)} attempted "
                "evaluations)",
                file=stream,
            )
            for failure in failures[:40]:
                print(f"- {failure}", file=stream)
            if len(failures) > 40:
                print(
                    f"- ... {len(failures) - 40} additional exclusions; "
                    "see CSV validity counts.",
                    file=stream,
                )
        else:
            print("\nDiagnostic exclusions: none", file=stream)


def synthetic_profile_ensemble(
    n_cfg: int,
) -> np.ndarray:
    """Build a positive-definite three-state ensemble for integration tests."""
    rng = np.random.default_rng(20260813)
    result = np.empty(
        (
            n_cfg,
            delta_analysis.N_T,
            delta_analysis.N_R,
            N_PROFILES,
            N_PROFILES,
        ),
        dtype=complex,
    )
    base_overlap = np.array(
        [
            [1.00, 0.48 + 0.08j, 0.18],
            [0.72 - 0.05j, 0.92, 0.43 + 0.04j],
            [0.45, 0.62 - 0.06j, 1.08],
        ],
        dtype=complex,
    )
    for cfg in range(n_cfg):
        amplitude = np.exp(rng.normal(0.0, 0.025, N_PROFILES))
        for radius in range(delta_analysis.N_R):
            energies = np.array(
                [
                    0.36 + 0.105 * radius,
                    0.76 + 0.090 * radius,
                    1.18 + 0.075 * radius,
                ]
            )
            energies += rng.normal(0.0, 0.004, N_PROFILES)
            overlap = base_overlap * (
                1.0 + rng.normal(0.0, 0.006, base_overlap.shape)
            )
            for tau in range(1, delta_analysis.N_T + 1):
                result[cfg, tau - 1, radius] = (
                    overlap
                    @ np.diag(amplitude * np.exp(-energies * tau))
                    @ overlap.conj().T
                )
    return result


def render_self_test_artifacts(
    output_dir: Path,
    profile_matrix: np.ndarray,
    eigenvalues: np.ndarray,
    weights: np.ndarray,
    scalar_results: list[ScalarResult],
    conditioning_rows: list[dict[str, object]],
    principal_rows: list[dict[str, object]],
    projected_rows: list[dict[str, object]],
    vector_rows: list[dict[str, object]],
    reduced_conditioning_rows: list[dict[str, object]],
    reduced_principal_rows: list[dict[str, object]],
    reduced_projected_rows: list[dict[str, object]],
    reduced_vector_rows: list[dict[str, object]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "selftest_scalar.csv", scalar_rows(scalar_results))
    write_csv(output_dir / "selftest_conditioning.csv", conditioning_rows)
    write_csv(output_dir / "selftest_principal.csv", principal_rows)
    write_csv(output_dir / "selftest_projected.csv", projected_rows)
    write_csv(output_dir / "selftest_vectors.csv", vector_rows)
    write_csv(
        output_dir / "selftest_reduced_conditioning.csv",
        reduced_conditioning_rows,
    )
    write_csv(
        output_dir / "selftest_reduced_principal.csv",
        reduced_principal_rows,
    )
    write_csv(
        output_dir / "selftest_reduced_projected.csv",
        reduced_projected_rows,
    )
    write_csv(
        output_dir / "selftest_reduced_vectors.csv",
        reduced_vector_rows,
    )
    write_report(
        output_dir / "selftest.txt",
        "selftest",
        profile_matrix.shape[0],
        eigenvalues,
        weights,
        scalar_results,
        conditioning_rows,
        principal_rows,
        projected_rows,
        reduced_conditioning_rows,
        reduced_principal_rows,
        reduced_projected_rows,
        [],
        (1, 2, 3),
        2,
        3,
        1.0e-10,
    )
    plot_scalar_comparison(output_dir / "selftest_scalar.png", scalar_results)
    plot_ground_reference_scan(
        output_dir / "selftest_ground_tref_scan.png",
        principal_rows,
        scalar_results,
        (1, 2, 3),
    )
    plot_states(
        output_dir / "selftest_levels_tref2.png",
        principal_rows,
        "principal_effective_energy",
        2,
        scalar_results,
        "Synthetic three-profile GEVP principal effective energies",
        "Layout and uncertainty-rendering integration test.",
    )
    plot_reduced_states(
        output_dir / "selftest_reduced_levels_tref2.png",
        reduced_principal_rows,
        "principal_effective_energy",
        2,
        scalar_results,
        "Synthetic reduced-basis GEVP principal effective energies",
        "Hollow markers indicate incomplete delete-one jackknife support.",
    )
    plot_reduced_states(
        output_dir / "selftest_reduced_projected_tref2_td3.png",
        reduced_projected_rows,
        "projected_effective_energy",
        2,
        scalar_results,
        "Synthetic reduced-basis fixed-vector correlators",
        "Hollow markers indicate incomplete delete-one jackknife support.",
    )
    plot_states(
        output_dir / "selftest_projected_tref2_td3.png",
        projected_rows,
        "projected_effective_energy",
        2,
        scalar_results,
        "Synthetic fixed-vector optimized correlators",
        "Layout and uncertainty-rendering integration test.",
    )


def run_self_test(artifact_dir: Path | None = None) -> None:
    energies = np.array([0.35, 0.72, 1.15])
    overlap = np.array(
        [
            [1.0, 0.4 + 0.1j, 0.2],
            [0.3 - 0.2j, 1.1, 0.5],
            [0.15, 0.35 + 0.05j, 0.9],
        ],
        dtype=complex,
    )

    def spectral(time: int) -> np.ndarray:
        return overlap @ np.diag(np.exp(-energies * time)) @ overlap.conj().T

    tau_ref = 2
    tau = 5
    solution = solve_gevp(spectral(tau), spectral(tau_ref), 1.0e-12)
    expected = np.exp(-energies * (tau - tau_ref))
    expected = np.sort(expected)[::-1]
    if not np.allclose(solution.eigenvalues, expected, rtol=2.0e-12, atol=1.0e-13):
        raise AssertionError(
            f"GEVP eigenvalue self-test failed: {solution.eigenvalues} vs {expected}"
        )
    metric_identity = (
        solution.eigenvectors.conj().T
        @ spectral(tau_ref)
        @ solution.eigenvectors
    )
    if not np.allclose(metric_identity, np.eye(3), atol=2.0e-12):
        raise AssertionError("GEVP metric-normalization self-test failed")

    indefinite_metric = np.array([[1.0, 1.01], [1.01, 1.0]])
    indefinite_diagnostics = metric_diagnostics(indefinite_metric)
    if metric_status(indefinite_diagnostics, 1.0e-10) != "non_positive":
        raise AssertionError("indefinite metric status self-test failed")
    near_singular_metric = np.array(
        [[1.0, 1.0 - 1.0e-12], [1.0 - 1.0e-12, 1.0]]
    )
    near_singular_diagnostics = metric_diagnostics(near_singular_metric)
    if metric_status(near_singular_diagnostics, 1.0e-10) != "below_rcond":
        raise AssertionError("near-singular metric status self-test failed")

    rng = np.random.default_rng(20260812)
    n_cfg = 2
    data_shape = (
        n_cfg,
        delta_analysis.N_T,
        delta_analysis.N_R,
    )
    delta_shape = data_shape + (
        delta_analysis.N_V,
        delta_analysis.N_V,
    )
    kernel = rng.normal(size=delta_shape) + 1j * rng.normal(size=delta_shape)
    data = kernel.sum(axis=(-2, -1))
    values = np.linspace(0.04, 0.09, delta_analysis.N_V)
    eigenvalues = np.broadcast_to(
        values,
        (n_cfg, delta_analysis.N_T + 1, delta_analysis.N_V),
    ).copy()
    profiles, _ = project_profile_matrix(data, kernel, eigenvalues)
    if not np.allclose(profiles[..., 2, 2], data):
        raise AssertionError("constant-profile projection self-test failed")

    synthetic = synthetic_profile_ensemble(16)
    scalar_results, scalar_failures = analyse_scalars(synthetic)
    conditioning, principal, failures = analyse_gevp(
        synthetic, (1, 2, 3), 1.0e-10
    )
    projected, vectors, projected_failures = analyse_fixed_vectors(
        synthetic, 2, 3, 1.0e-10
    )
    reduced_conditioning, reduced_principal, reduced_failures = (
        analyse_reduced_gevps(synthetic, (1, 2, 3), 1.0e-10)
    )
    reduced_projected, reduced_vectors, reduced_projected_failures = (
        analyse_reduced_fixed_vectors(synthetic, 2, 3, 1.0e-10)
    )
    integration_failures = (
        scalar_failures
        + failures
        + projected_failures
        + reduced_failures
        + reduced_projected_failures
    )
    if integration_failures:
        raise AssertionError(
            "synthetic integration analysis unexpectedly failed: "
            + "; ".join(integration_failures[:3])
        )
    if len(conditioning) != 6 * 3:
        raise AssertionError("unexpected conditioning-row count")
    if len(principal) != 270:
        raise AssertionError("unexpected principal-row count")
    if len(projected) != 6 * 7 * 3:
        raise AssertionError("unexpected projected-row count")
    if len(vectors) != 6 * 3 * 3:
        raise AssertionError("unexpected optimized-vector-row count")
    if len(reduced_conditioning) != 2 * 6 * 3:
        raise AssertionError("unexpected reduced-conditioning-row count")
    if len(reduced_principal) != 2 * 6 * 30:
        raise AssertionError("unexpected reduced-principal-row count")
    if len(reduced_projected) != 2 * 6 * 7 * 2:
        raise AssertionError("unexpected reduced-projected-row count")
    if len(reduced_vectors) != 2 * 6 * 2 * 2:
        raise AssertionError("unexpected reduced-vector-row count")

    if artifact_dir is not None:
        synthetic_values = np.broadcast_to(
            values,
            (
                synthetic.shape[0],
                delta_analysis.N_T + 1,
                delta_analysis.N_V,
            ),
        ).copy()
        synthetic_weights = endpoint_weights(synthetic_values)
        render_self_test_artifacts(
            artifact_dir,
            synthetic,
            synthetic_values,
            synthetic_weights,
            scalar_results,
            conditioning,
            principal,
            projected,
            vectors,
            reduced_conditioning,
            reduced_principal,
            reduced_projected,
            reduced_vectors,
        )
    print(
        "PASS: analytic GEVP, profile projection, jackknife analysis, "
        "and output integration self-tests"
    )


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_test(args.self_test_artifact_dir)
        return

    reference_times = validate_reference_times(args.gevp_reference_times)
    if args.plot_reference_time not in reference_times:
        raise ValueError(
            "--plot-reference-time must be included in "
            "--gevp-reference-times"
        )
    if not (
        args.plot_reference_time < args.optimization_time <= delta_analysis.N_T
    ):
        raise ValueError(
            "optimization time must satisfy plot_reference_time < "
            f"optimization_time <= {delta_analysis.N_T}"
        )
    if not (0.0 < args.gevp_rcond < 1.0):
        raise ValueError("--gevp-rcond must lie strictly between zero and one")

    outputs = delta_analysis.discover_outputs(args.delta_dir, args.jobid)
    data, delta_matrix = delta_analysis.read_delta_outputs(outputs)
    eigenvalues = delta_analysis.read_eigenvalues(
        args.eigenvalue_root, len(outputs)
    )
    profile_matrix, weights = project_profile_matrix(
        data, delta_matrix, eigenvalues
    )
    scalar_results, scalar_failures = analyse_scalars(profile_matrix)
    conditioning_rows, principal_rows, gevp_failures = analyse_gevp(
        profile_matrix, reference_times, args.gevp_rcond
    )
    projected_rows, vector_rows, projected_failures = analyse_fixed_vectors(
        profile_matrix,
        args.plot_reference_time,
        args.optimization_time,
        args.gevp_rcond,
    )
    reduced_conditioning_rows, reduced_principal_rows, reduced_failures = (
        analyse_reduced_gevps(
            profile_matrix,
            reference_times,
            args.gevp_rcond,
        )
    )
    reduced_projected_rows, reduced_vector_rows, reduced_projected_failures = (
        analyse_reduced_fixed_vectors(
            profile_matrix,
            args.plot_reference_time,
            args.optimization_time,
            args.gevp_rcond,
        )
    )
    failures = (
        scalar_failures
        + gevp_failures
        + projected_failures
        + reduced_failures
        + reduced_projected_failures
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"gaussian_gevp_{args.jobid}"
    paths = {
        "scalar_csv": args.output_dir / f"{stem}_scalar.csv",
        "conditioning_csv": args.output_dir / f"{stem}_conditioning.csv",
        "principal_csv": args.output_dir / f"{stem}_principal.csv",
        "projected_csv": args.output_dir / f"{stem}_projected.csv",
        "vectors_csv": args.output_dir / f"{stem}_vectors.csv",
        "reduced_conditioning_csv": args.output_dir
        / f"{stem}_reduced_conditioning.csv",
        "reduced_principal_csv": args.output_dir
        / f"{stem}_reduced_principal.csv",
        "reduced_projected_csv": args.output_dir
        / f"{stem}_reduced_projected.csv",
        "reduced_vectors_csv": args.output_dir
        / f"{stem}_reduced_vectors.csv",
        "report": args.output_dir / f"{stem}.txt",
        "scalar_plot": args.output_dir / f"{stem}_scalar.png",
        "ground_plot": args.output_dir / f"{stem}_ground_tref_scan.png",
        "levels_plot": args.output_dir
        / f"{stem}_levels_tref{args.plot_reference_time}.png",
        "projected_plot": args.output_dir
        / (
            f"{stem}_projected_tref{args.plot_reference_time}"
            f"_td{args.optimization_time}.png"
        ),
        "reduced_levels_plot": args.output_dir
        / f"{stem}_reduced_levels_tref{args.plot_reference_time}.png",
        "reduced_projected_plot": args.output_dir
        / (
            f"{stem}_reduced_projected_tref{args.plot_reference_time}"
            f"_td{args.optimization_time}.png"
        ),
    }

    write_csv(paths["scalar_csv"], scalar_rows(scalar_results))
    write_csv(paths["conditioning_csv"], conditioning_rows)
    write_csv(paths["principal_csv"], principal_rows)
    write_csv(
        paths["projected_csv"],
        projected_rows,
        [
            "radius",
            "tau_ref",
            "tau_d",
            "effective_time",
            "state",
            "projected_effective_energy",
            "jackknife_error",
            "jackknife_valid",
            "jackknife_total",
        ],
    )
    write_csv(
        paths["vectors_csv"],
        vector_rows,
        [
            "radius",
            "tau_ref",
            "tau_d",
            "state",
            "profile",
            "coefficient_real",
            "coefficient_imag",
        ],
    )
    write_csv(paths["reduced_conditioning_csv"], reduced_conditioning_rows)
    write_csv(paths["reduced_principal_csv"], reduced_principal_rows)
    write_csv(
        paths["reduced_projected_csv"],
        reduced_projected_rows,
        [
            "basis",
            "basis_profiles",
            "basis_size",
            "radius",
            "tau_ref",
            "tau_d",
            "effective_time",
            "state",
            "projected_effective_energy",
            "jackknife_error",
            "jackknife_valid",
            "jackknife_total",
        ],
    )
    write_csv(
        paths["reduced_vectors_csv"],
        reduced_vector_rows,
        [
            "basis",
            "basis_profiles",
            "basis_size",
            "radius",
            "tau_ref",
            "tau_d",
            "state",
            "profile",
            "coefficient_real",
            "coefficient_imag",
        ],
    )
    write_report(
        paths["report"],
        args.jobid,
        len(outputs),
        eigenvalues,
        weights,
        scalar_results,
        conditioning_rows,
        principal_rows,
        projected_rows,
        reduced_conditioning_rows,
        reduced_principal_rows,
        reduced_projected_rows,
        failures,
        reference_times,
        args.plot_reference_time,
        args.optimization_time,
        args.gevp_rcond,
    )
    plot_scalar_comparison(paths["scalar_plot"], scalar_results)
    plot_ground_reference_scan(
        paths["ground_plot"],
        principal_rows,
        scalar_results,
        reference_times,
    )
    plot_states(
        paths["levels_plot"],
        principal_rows,
        "principal_effective_energy",
        args.plot_reference_time,
        scalar_results,
        (
            "Three-profile GEVP principal effective energies "
            rf"($\tau_{{\rm ref}}/a={args.plot_reference_time}$)"
        ),
        (
            "States are ordered by descending principal correlator. "
            "Excited-state interpretation requires reference-time stability."
        ),
    )
    plot_reduced_states(
        paths["reduced_levels_plot"],
        reduced_principal_rows,
        "principal_effective_energy",
        args.plot_reference_time,
        scalar_results,
        (
            "Reduced-basis GEVP principal effective energies "
            rf"($\tau_{{\rm ref}}/a={args.plot_reference_time}$)"
        ),
        (
            "Solid/dashed markers denote n=1/n=2. Hollow markers have "
            "incomplete delete-one jackknife support."
        ),
    )
    plot_reduced_states(
        paths["reduced_projected_plot"],
        reduced_projected_rows,
        "projected_effective_energy",
        args.plot_reference_time,
        scalar_results,
        (
            "Reduced-basis fixed-vector optimized correlators "
            rf"($\tau_{{\rm ref}}/a={args.plot_reference_time}$, "
            rf"$\tau_d/a={args.optimization_time}$)"
        ),
        (
            "Solid/dashed markers denote n=1/n=2. Hollow markers have "
            "incomplete delete-one jackknife support."
        ),
    )
    plot_states(
        paths["projected_plot"],
        projected_rows,
        "projected_effective_energy",
        args.plot_reference_time,
        scalar_results,
        (
            "Fixed-vector optimized correlators "
            rf"($\tau_{{\rm ref}}/a={args.plot_reference_time}$, "
            rf"$\tau_d/a={args.optimization_time}$)"
        ),
        (
            r"Diagonal channels of $V^\dagger C(\tau)V$; $n=2$ is the "
            "candidate first-excited optimized correlator."
        ),
    )

    print(f"Validated configurations: {len(outputs)}")
    print("Profiles: sigma=0.01, sigma=0.05, sigma=infinity")
    print(
        "Maximum constant reconstruction residual: "
        f"{np.max(np.abs(profile_matrix[..., 2, 2] - data)):.3e}"
    )
    print(f"GEVP diagnostic exclusions recorded: {len(failures)}")
    for path in paths.values():
        print(f"Wrote: {path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, AssertionError) as error:
        raise SystemExit(f"ERROR: {error}") from error
