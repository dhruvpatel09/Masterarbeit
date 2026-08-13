#!/usr/bin/env python3
"""Plot the early-plateau comparison for Gaussian Laplace profiles.

This companion to ``inspect_laplace_gaussian_from_delta.py`` reads the same
validated PROFILE_BASIS=delta outputs and reconstructs the constant profile
and the seven Gaussian profiles on exactly the same configurations.

It writes:

* a detailed CSV of V_eff(T), plateau fits, and paired jackknife residuals;
* a per-profile decision CSV;
* a six-panel early-plateau figure;
* a two-panel ratio heatmap with one shared colour scale; and
* a short plain-text interpretation.

No gauge fields or eigenvectors are reread by the C scanner, and no new Slurm
job is required.  The only additional input is the already stored eigenvalue
inventory used by the projection analyser.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

try:
    import inspect_laplace_gaussian_from_delta as projection
except ImportError as error:
    raise SystemExit(
        "Could not import inspect_laplace_gaussian_from_delta.py. "
        "Place this script in the same laplace_static/scripts directory."
    ) from error


PLOT_TIMES = np.array([2, 3, 4, 5], dtype=int)
DEFAULT_CANDIDATE_SIGMA = 0.0500


@dataclass
class ProfileRadius:
    """Central values and delete-one jackknife samples for one profile/R."""

    label: str
    sigma: float | None
    radius: int
    effective: np.ndarray
    effective_error: np.ndarray
    residual: np.ndarray
    residual_error: np.ndarray
    fit_result: projection.FitResult

    @property
    def delta2(self) -> float:
        return float(self.residual[projection.EARLY_TIME - 1])

    @property
    def delta2_error(self) -> float:
        return float(self.residual_error[projection.EARLY_TIME - 1])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an early-plateau profile comparison from a "
            "validated fixed-t0 PROFILE_BASIS=delta scan."
        )
    )
    parser.add_argument(
        "delta_dir",
        type=Path,
        help="directory containing cfg<jobid>_<cfg>.out files",
    )
    parser.add_argument(
        "jobid",
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
        default=None,
        help="output directory (default: the delta scan directory)",
    )
    parser.add_argument(
        "--candidate-sigma",
        type=float,
        default=DEFAULT_CANDIDATE_SIGMA,
        help=(
            "Gaussian width emphasized in the direct plateau plot "
            "(default: %(default).4f)"
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="PNG resolution (default: %(default)s)",
    )
    return parser.parse_args()


def sigma_from_label(label: str) -> float | None:
    if label == "constant":
        return None
    prefix = "sigma="
    if not label.startswith(prefix):
        raise ValueError(f"unrecognized profile label: {label}")
    return float(label[len(prefix) :])


def jackknife_effective(
    correlator: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return central and leave-one-out effective energies."""

    real = np.real(correlator)
    n_cfg = real.shape[0]
    total = real.sum(axis=0)
    mean = total / n_cfg
    leave_one_out = (total[None, :] - real) / (n_cfg - 1)

    effective = np.full(projection.N_T - 1, np.nan)
    valid_mean = (mean[:-1] > 0) & (mean[1:] > 0)
    effective[valid_mean] = np.log(
        mean[:-1][valid_mean] / mean[1:][valid_mean]
    )

    jk_effective = np.full(
        (n_cfg, projection.N_T - 1),
        np.nan,
    )
    valid_jk = (
        (leave_one_out[:, :-1] > 0)
        & (leave_one_out[:, 1:] > 0)
    )
    jk_effective[valid_jk] = np.log(
        leave_one_out[:, :-1][valid_jk]
        / leave_one_out[:, 1:][valid_jk]
    )
    return effective, jk_effective


def analyze(
    labels: list[str],
    correlators: np.ndarray,
) -> list[ProfileRadius]:
    analyses: list[ProfileRadius] = []
    fit_results: list[projection.FitResult] = []

    for profile_index, label in enumerate(labels):
        for radius in range(1, projection.N_R):
            correlator = correlators[:, :, radius, profile_index]
            fit_result = projection.fit_profile(
                label,
                radius,
                correlator,
            )
            effective, jk_effective = jackknife_effective(correlator)

            effective_error = np.full_like(effective, np.nan)
            residual_error = np.full_like(effective, np.nan)
            residual = effective - fit_result.fit

            for time in PLOT_TIMES:
                index = time - 1
                if not (
                    np.isfinite(effective[index])
                    and np.all(np.isfinite(jk_effective[:, index]))
                ):
                    raise ValueError(
                        f"{label}, R={radius}, T={time}: "
                        "non-finite effective energy"
                    )
                effective_error[index] = projection.jackknife_error(
                    jk_effective[:, index]
                )
                residual_error[index] = projection.jackknife_error(
                    jk_effective[:, index] - fit_result.jk_fit
                )

            analyses.append(
                ProfileRadius(
                    label=label,
                    sigma=sigma_from_label(label),
                    radius=radius,
                    effective=effective,
                    effective_error=effective_error,
                    residual=residual,
                    residual_error=residual_error,
                    fit_result=fit_result,
                )
            )
            fit_results.append(fit_result)

    projection.add_constant_comparisons(fit_results)
    return analyses


def index_analyses(
    analyses: list[ProfileRadius],
) -> dict[tuple[str, int], ProfileRadius]:
    return {
        (analysis.label, analysis.radius): analysis
        for analysis in analyses
    }


def write_points_csv(
    path: Path,
    analyses: list[ProfileRadius],
) -> None:
    fields = [
        "profile",
        "sigma",
        "R_over_a",
        "T_over_a",
        "Veff",
        "Veff_jackknife_error",
        "plateau_fit_T3_5",
        "plateau_fit_error",
        "Veff_minus_own_plateau",
        "paired_residual_jackknife_error",
        "z2",
        "plateau_p_value",
        "plateau_energy_shift_from_constant",
        "paired_shift_z_from_constant",
        "in_plateau_fit",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for analysis in analyses:
            fit = analysis.fit_result
            for time in PLOT_TIMES:
                index = time - 1
                writer.writerow(
                    {
                        "profile": analysis.label,
                        "sigma": (
                            ""
                            if analysis.sigma is None
                            else f"{analysis.sigma:.4f}"
                        ),
                        "R_over_a": analysis.radius,
                        "T_over_a": time,
                        "Veff": f"{analysis.effective[index]:.12g}",
                        "Veff_jackknife_error": (
                            f"{analysis.effective_error[index]:.12g}"
                        ),
                        "plateau_fit_T3_5": f"{fit.fit:.12g}",
                        "plateau_fit_error": f"{fit.error:.12g}",
                        "Veff_minus_own_plateau": (
                            f"{analysis.residual[index]:.12g}"
                        ),
                        "paired_residual_jackknife_error": (
                            f"{analysis.residual_error[index]:.12g}"
                        ),
                        "z2": f"{fit.early_z:.12g}",
                        "plateau_p_value": f"{fit.p_value:.12g}",
                        "plateau_energy_shift_from_constant": (
                            f"{fit.delta_constant:.12g}"
                        ),
                        "paired_shift_z_from_constant": (
                            f"{fit.delta_constant_z:.12g}"
                        ),
                        "in_plateau_fit": (
                            "yes"
                            if time in projection.PLATEAU_TIMES
                            else "no"
                        ),
                    }
                )


def summarize_profiles(
    labels: list[str],
    analyses: list[ProfileRadius],
) -> list[dict[str, object]]:
    indexed = index_analyses(analyses)
    radii = range(1, projection.N_R)
    rows: list[dict[str, object]] = []

    for label in labels:
        selected = [indexed[(label, radius)] for radius in radii]
        early_gaps = np.array(
            [abs(item.delta2) for item in selected]
        )
        early_z = np.array(
            [abs(item.fit_result.early_z) for item in selected]
        )
        relative_fit_error = np.array(
            [
                item.fit_result.error / abs(item.fit_result.fit)
                for item in selected
            ]
        )

        if label == "constant":
            early_better: list[int] = []
            precision_better: list[int] = []
            both_better: list[int] = []
            mean_early_gap_ratio = 1.0
            mean_fit_error_ratio = 1.0
        else:
            early_better = []
            precision_better = []
            both_better = []
            early_ratios = []
            fit_error_ratios = []
            for item in selected:
                reference = indexed[("constant", item.radius)]
                early_is_better = (
                    abs(item.delta2) < abs(reference.delta2)
                )
                precision_is_better = (
                    item.fit_result.error
                    < reference.fit_result.error
                )
                if early_is_better:
                    early_better.append(item.radius)
                if precision_is_better:
                    precision_better.append(item.radius)
                if early_is_better and precision_is_better:
                    both_better.append(item.radius)
                early_ratios.append(
                    abs(item.delta2) / abs(reference.delta2)
                )
                fit_error_ratios.append(
                    item.fit_result.error
                    / reference.fit_result.error
                )
            mean_early_gap_ratio = float(np.mean(early_ratios))
            mean_fit_error_ratio = float(
                np.mean(fit_error_ratios)
            )

        rows.append(
            {
                "profile": label,
                "sigma": (
                    ""
                    if label == "constant"
                    else f"{sigma_from_label(label):.4f}"
                ),
                "mean_abs_delta2": float(np.mean(early_gaps)),
                "median_abs_delta2": float(np.median(early_gaps)),
                "median_abs_z2": float(np.median(early_z)),
                "median_relative_plateau_error": float(
                    np.median(relative_fit_error)
                ),
                "mean_paired_early_gap_ratio_to_constant": (
                    mean_early_gap_ratio
                ),
                "mean_paired_fit_error_ratio_to_constant": (
                    mean_fit_error_ratio
                ),
                "early_gap_better_radii": early_better,
                "fit_error_better_radii": precision_better,
                "both_better_radii": both_better,
                "max_abs_paired_plateau_shift_z": float(
                    max(
                        abs(item.fit_result.delta_constant_z)
                        for item in selected
                    )
                ),
            }
        )
    return rows


def radii_text(radii: list[int]) -> str:
    return ";".join(str(radius) for radius in radii)


def write_summary_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    fields = [
        "profile",
        "sigma",
        "mean_abs_delta2",
        "median_abs_delta2",
        "median_abs_z2",
        "median_relative_plateau_error",
        "mean_paired_early_gap_ratio_to_constant",
        "mean_paired_fit_error_ratio_to_constant",
        "early_gap_better_radii",
        "fit_error_better_radii",
        "both_better_radii",
        "max_abs_paired_plateau_shift_z",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            output = dict(row)
            for key in (
                "mean_abs_delta2",
                "median_abs_delta2",
                "median_abs_z2",
                "median_relative_plateau_error",
                "mean_paired_early_gap_ratio_to_constant",
                "mean_paired_fit_error_ratio_to_constant",
                "max_abs_paired_plateau_shift_z",
            ):
                output[key] = f"{float(output[key]):.12g}"
            for key in (
                "early_gap_better_radii",
                "fit_error_better_radii",
                "both_better_radii",
            ):
                output[key] = radii_text(output[key])
            writer.writerow(output)


def choose_candidate_label(
    labels: list[str],
    sigma: float,
) -> str:
    candidates = [
        label
        for label in labels
        if label != "constant"
        and np.isclose(
            sigma_from_label(label),
            sigma,
            rtol=0.0,
            atol=5.0e-5,
        )
    ]
    if len(candidates) != 1:
        available = ", ".join(labels[1:])
        raise ValueError(
            f"candidate sigma={sigma:.4f} not found uniquely; "
            f"available profiles: {available}"
        )
    return candidates[0]


def plot_plateau_approach(
    path: Path,
    labels: list[str],
    analyses: list[ProfileRadius],
    candidate_label: str,
    n_cfg: int,
    dpi: int,
) -> None:
    indexed = index_analyses(analyses)
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(14.5, 8.8),
        sharex=True,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.985,
        bottom=0.085,
        top=0.79,
        wspace=0.18,
        hspace=0.28,
    )

    constant_color = "#15284b"
    candidate_color = "#d55e00"
    other_facecolor = "#8d96a0"
    other_edgecolor = "#69727d"
    plot_indices = PLOT_TIMES - 1

    for radius, axis in zip(range(1, projection.N_R), axes.flat):
        axis.axvspan(
            2.5,
            5.5,
            color="#dfe7ef",
            alpha=0.55,
            zorder=0,
            label="fit window" if radius == 1 else None,
        )
        axis.axhline(0.0, color="#56616f", linewidth=0.9, zorder=1)

        other_residuals = np.vstack(
            [
                indexed[(label, radius)].residual[plot_indices]
                for label in labels[1:]
                if label != candidate_label
            ]
        )
        other_min = np.min(other_residuals, axis=0)
        other_max = np.max(other_residuals, axis=0)
        axis.fill_between(
            PLOT_TIMES,
            other_min,
            other_max,
            facecolor=other_facecolor,
            edgecolor=other_edgecolor,
            alpha=0.34,
            linewidth=1.0,
            zorder=2,
            label=(
                "other Gaussian widths\n(central-value range)"
                if radius == 1
                else None
            ),
        )
        axis.plot(
            PLOT_TIMES,
            other_min,
            color=other_edgecolor,
            alpha=0.80,
            linewidth=0.75,
            zorder=2.2,
        )
        axis.plot(
            PLOT_TIMES,
            other_max,
            color=other_edgecolor,
            alpha=0.80,
            linewidth=0.75,
            zorder=2.2,
        )

        for label, color, marker, display in (
            ("constant", constant_color, "o", r"$\rho_i=1$"),
            (
                candidate_label,
                candidate_color,
                "s",
                rf"$\sigma={sigma_from_label(candidate_label):.4f}$",
            ),
        ):
            item = indexed[(label, radius)]
            axis.errorbar(
                PLOT_TIMES,
                item.residual[plot_indices],
                yerr=item.residual_error[plot_indices],
                color=color,
                marker=marker,
                markersize=4.5,
                linewidth=1.45,
                elinewidth=1.0,
                capsize=2.4,
                label=display if radius == 1 else None,
                zorder=4,
            )

        reference = indexed[("constant", radius)]
        candidate = indexed[(candidate_label, radius)]
        annotation = (
            rf"$\Delta_2$: {reference.delta2:+.4f} "
            rf"({reference.fit_result.early_z:+.2f}$\sigma$)"
            "\n"
            rf"$\Delta_2^{{\sigma={sigma_from_label(candidate_label):.4f}}}$: "
            rf"{candidate.delta2:+.4f} "
            rf"({candidate.fit_result.early_z:+.2f}$\sigma$)"
        )
        axis.text(
            0.035,
            0.045,
            annotation,
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.6,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#ccd2d8",
                "alpha": 0.88,
            },
        )
        axis.set_title(rf"$R/a={radius}$", fontsize=11.5)
        axis.set_xticks(PLOT_TIMES)
        axis.grid(axis="y", color="#e4e7eb", linewidth=0.7)
        axis.tick_params(labelsize=9)

    for axis in axes[:, 0]:
        axis.set_ylabel(
            r"$V_{\rm eff}(T)-\widehat V_{3:5}$",
            fontsize=10.5,
        )
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$T/a$", fontsize=10.5)

    handles, legend_labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.865),
    )
    figure.suptitle(
        "Early approach to the correlated plateau",
        fontsize=16,
        fontweight="semibold",
        y=0.975,
    )
    figure.text(
        0.5,
        0.925,
        (
            rf"Fixed $t_0=0$, $x$-oriented pilot, "
            rf"$N_{{\rm cfg}}={n_cfg}$. Each profile is centred on its "
            rf"own correlated $T/a=3\ldots5$ fit; the gray band spans "
            rf"the central values of the other Gaussian widths."
        ),
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#424b55",
    )
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_ratio_heatmap(
    path: Path,
    labels: list[str],
    analyses: list[ProfileRadius],
    dpi: int,
) -> None:
    indexed = index_analyses(analyses)
    gaussian_labels = labels[1:]
    radii = list(range(1, projection.N_R))
    early_ratio = np.empty((len(gaussian_labels), len(radii)))
    error_ratio = np.empty_like(early_ratio)

    for row, label in enumerate(gaussian_labels):
        for column, radius in enumerate(radii):
            item = indexed[(label, radius)]
            reference = indexed[("constant", radius)]
            early_ratio[row, column] = (
                abs(item.delta2) / abs(reference.delta2)
            )
            error_ratio[row, column] = (
                item.fit_result.error
                / reference.fit_result.error
            )

    values = np.concatenate((early_ratio.ravel(), error_ratio.ravel()))
    vmin = min(float(values.min()), 0.98)
    vmax = max(float(values.max()), 1.02)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(14.8, 6.4),
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.90,
        bottom=0.17,
        top=0.85,
        wspace=0.18,
    )
    titles = (
        (
            r"Early-gap ratio "
            r"$|\Delta_{2,\sigma}|/|\Delta_{2,\rho=1}|$"
        ),
        (
            r"Plateau-error ratio "
            r"$\delta\widehat V_\sigma/"
            r"\delta\widehat V_{\rho=1}$"
        ),
    )

    image = None
    for axis, matrix, title in zip(
        axes,
        (early_ratio, error_ratio),
        titles,
    ):
        image = axis.imshow(
            matrix,
            cmap="RdYlGn_r",
            norm=norm,
            aspect="auto",
        )
        axis.set_title(title, fontsize=12)
        axis.set_xticks(np.arange(len(radii)), radii)
        axis.set_yticks(
            np.arange(len(gaussian_labels)),
            [
                rf"$\sigma={sigma_from_label(label):.4f}$"
                for label in gaussian_labels
            ],
        )
        axis.set_xlabel(r"$R/a$")
        axis.tick_params(labelsize=9)

        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix[row, column]
                axis.text(
                    column,
                    row,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=(
                        "white"
                        if norm(value) < 0.12 or norm(value) > 0.88
                        else "#202124"
                    ),
                )
    if image is None:
        raise RuntimeError("ratio heatmap has no panels")
    colorbar = figure.colorbar(
        image,
        ax=axes.ravel().tolist(),
        shrink=0.82,
        pad=0.025,
    )
    colorbar.set_label("ratio to constant profile")

    figure.suptitle(
        "Gaussian-to-constant diagnostic ratios",
        fontsize=16,
        fontweight="semibold",
        y=0.97,
    )
    figure.text(
        0.5,
        0.04,
        (
            "Values below 1 (green) are numerically smaller than for "
            "the constant profile; this alone does not establish "
            "statistical significance.\nA simultaneous numerical "
            "improvement requires both panels below 1 for the same "
            "width and radius."
        ),
        ha="center",
        va="bottom",
        fontsize=9.5,
        color="#424b55",
        linespacing=1.25,
    )
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def write_interpretation(
    path: Path,
    jobid: str,
    n_cfg: int,
    labels: list[str],
    analyses: list[ProfileRadius],
    summary_rows: list[dict[str, object]],
    candidate_label: str,
) -> None:
    indexed = index_analyses(analyses)
    summary_by_label = {
        str(row["profile"]): row for row in summary_rows
    }
    constant = summary_by_label["constant"]
    candidate = summary_by_label[candidate_label]

    lines = [
        "GAUSSIAN EARLY-PLATEAU COMPARISON",
        f"Source job: {jobid}",
        f"Configurations: {n_cfg}",
        (
            "Scope: fixed t0=0, x-oriented separation, "
            "R/a=1..6, plateau T/a=3..5"
        ),
        "",
        "Definition",
        (
            "Delta2(R) = V_eff(T/a=2,R) - Vhat_3:5(R). "
            "A smaller |Delta2| is a closer raw early-time approach."
        ),
        (
            "z2 = Delta2 / jackknife_error(Delta2). "
            "A smaller |z2| can come either from a smaller gap or "
            "from a larger uncertainty, so it is not used alone."
        ),
        (
            "Statistical precision is represented by the correlated "
            "plateau-fit error delta Vhat."
        ),
        "",
        "Constant profile versus the most competitive Gaussian",
        (
            f"constant: median |z2|="
            f"{float(constant['median_abs_z2']):.3f}, "
            f"median relative plateau error="
            f"{100.0 * float(constant['median_relative_plateau_error']):.3f}%"
        ),
        (
            f"{candidate_label}: median |z2|="
            f"{float(candidate['median_abs_z2']):.3f}, "
            f"median relative plateau error="
            f"{100.0 * float(candidate['median_relative_plateau_error']):.3f}%"
        ),
        "",
        "Per-radius Delta2 and plateau-fit errors",
        (
            "R  Delta2(const)  Delta2(candidate)  "
            "err_fit(const)  err_fit(candidate)  both_better"
        ),
    ]
    for radius in range(1, projection.N_R):
        reference = indexed[("constant", radius)]
        item = indexed[(candidate_label, radius)]
        both = (
            abs(item.delta2) < abs(reference.delta2)
            and item.fit_result.error < reference.fit_result.error
        )
        lines.append(
            f"{radius:d}  "
            f"{reference.delta2:+.8f}  "
            f"{item.delta2:+.8f}  "
            f"{reference.fit_result.error:.8f}  "
            f"{item.fit_result.error:.8f}  "
            f"{'yes' if both else 'no'}"
        )

    lines.extend(
        [
            "",
            "All Gaussian widths",
        ]
    )
    any_systematic = False
    for label in labels[1:]:
        row = summary_by_label[label]
        both = list(row["both_better_radii"])
        if len(both) == projection.N_R - 1:
            any_systematic = True
        lines.append(
            f"{label}: simultaneous numerical improvement at "
            f"R/a={radii_text(both) if both else 'none'}"
        )

    low_radius_improvement = any(
        any(radius <= 4 for radius in row["both_better_radii"])
        for row in summary_rows
        if row["profile"] != "constant"
    )
    lines.extend(
        [
            "",
            "Interpretation",
            (
                "No Gaussian improves both diagnostics at all six radii."
                if not any_systematic
                else (
                    "At least one Gaussian improves both numerical "
                    "diagnostics at all six radii."
                )
            ),
            (
                "No Gaussian gives a simultaneous numerical improvement "
                "at R/a=1..4, where the constant-profile early-time "
                "mismatch is most visible."
                if not low_radius_improvement
                else (
                    "At least one Gaussian gives a simultaneous numerical "
                    "improvement within R/a=1..4."
                )
            ),
            (
                "At R/a=5,6 the constant profile already has |z2|<1, "
                "so small favorable numerical changes there are not "
                "evidence for an earlier plateau."
            ),
            (
                "Conclusion: no clear, systematic, statistically "
                "meaningful simultaneous improvement over rho_i=1 was "
                "found for the tested Nv=10 profiles."
            ),
            "",
            (
                "This is a paired comparison inside the current fixed-t0 "
                "pilot. The earlier t0-averaged, axis-averaged production "
                "is not overlaid because it has different averaging and "
                "therefore is not the correct head-to-head reference."
            ),
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    output_dir = (
        args.delta_dir
        if args.output_dir is None
        else args.output_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = projection.discover_outputs(args.delta_dir, args.jobid)
    data, matrix = projection.read_delta_outputs(outputs)
    eigenvalues = projection.read_eigenvalues(
        args.eigenvalue_root,
        len(outputs),
    )
    labels, correlators = projection.project_profiles(
        data,
        matrix,
        eigenvalues,
    )
    analyses = analyze(labels, correlators)
    candidate_label = choose_candidate_label(
        labels,
        args.candidate_sigma,
    )
    summary_rows = summarize_profiles(labels, analyses)

    stem = f"gaussian_plateau_comparison_{args.jobid}"
    points_csv = output_dir / f"{stem}_points.csv"
    summary_csv = output_dir / f"{stem}_summary.csv"
    plateau_png = output_dir / f"{stem}.png"
    heatmap_png = output_dir / f"{stem}_ratios.png"
    interpretation_txt = output_dir / f"{stem}.txt"

    write_points_csv(points_csv, analyses)
    write_summary_csv(summary_csv, summary_rows)
    plot_plateau_approach(
        plateau_png,
        labels,
        analyses,
        candidate_label,
        len(outputs),
        args.dpi,
    )
    plot_ratio_heatmap(
        heatmap_png,
        labels,
        analyses,
        args.dpi,
    )
    write_interpretation(
        interpretation_txt,
        args.jobid,
        len(outputs),
        labels,
        analyses,
        summary_rows,
        candidate_label,
    )

    print(f"Validated configurations: {len(outputs)}")
    print(f"Candidate profile: {candidate_label}")
    for path in (
        points_csv,
        summary_csv,
        plateau_png,
        heatmap_png,
        interpretation_txt,
    ):
        print(f"Wrote: {path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
