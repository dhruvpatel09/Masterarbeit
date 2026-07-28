#!/usr/bin/env python3
"""Project smooth Gaussian Laplace profiles from a mode-delta scan.

The mode-delta output contains the complete Nv x Nv endpoint-weight kernel.
Together with the stored Laplacian eigenvalues, it therefore determines the
correlator for every Gaussian profile without another gauge-field calculation.

For the constant profile and the seven Gaussian widths used in the weighted
scanner, this program reports jackknife effective energies, correlated
T/a=3..5 plateau fits, early-time displacement, and agreement with the
constant-profile plateau.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re

import h5py
import numpy as np


N_T = 8
N_R = 7
N_V = 10
N_SRC = 24**3
PLATEAU_TIMES = (3, 4, 5)
EARLY_TIME = 2
SIGMAS = np.array(
    [0.0500, 0.0894, 0.1289, 0.1683, 0.2078, 0.2472, 0.2867],
    dtype=float,
)


@dataclass
class FitResult:
    profile: str
    radius: int
    effective: np.ndarray
    fit: float
    error: float
    chi2_dof: float
    p_value: float
    early_z: float
    jk_fit: np.ndarray
    delta_constant: float = np.nan
    delta_constant_z: float = np.nan


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Project the constant and seven Gaussian profiles from a "
            "fixed-t0 PROFILE_BASIS=delta scan."
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
    return parser.parse_args()


def discover_outputs(root: Path, jobid: str) -> dict[int, Path]:
    prefix = f"cfg{jobid}_"
    outputs: dict[int, Path] = {}

    for path in root.glob(f"{prefix}*.out"):
        suffix = path.name[len(prefix) : -len(".out")]
        if not suffix.isdigit():
            continue
        cfg = int(suffix)
        if cfg in outputs:
            raise ValueError(f"duplicate output for configuration {cfg}")
        outputs[cfg] = path

    if not outputs:
        raise ValueError(f"{root}: no outputs found for job {jobid}")

    expected = list(range(1, max(outputs) + 1))
    if sorted(outputs) != expected:
        missing = sorted(set(expected) - set(outputs))
        raise ValueError(
            f"{root}: non-contiguous configuration inventory; "
            f"missing={missing}"
        )

    return outputs


def parse_meta(fields: list[str]) -> dict[str, str]:
    result = {}
    for field in fields:
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        result[key] = value
    return result


def read_delta_outputs(
    outputs: dict[int, Path],
) -> tuple[np.ndarray, np.ndarray]:
    n_cfg = len(outputs)
    data = np.full(
        (n_cfg, N_T, N_R),
        np.nan + 1j * np.nan,
    )
    matrix = np.full(
        (n_cfg, N_T, N_R, N_V, N_V),
        np.nan + 1j * np.nan,
    )

    for cfg in range(1, n_cfg + 1):
        path = outputs[cfg]
        n_meta = n_data = n_weighted = 0

        with path.open() as stream:
            for lineno, line in enumerate(stream, 1):
                fields = line.split()
                if not fields:
                    continue

                if fields[0] == "META":
                    meta = parse_meta(fields[1:])
                    expected_meta = {
                        "cfg": str(cfg),
                        "Nv": str(N_V),
                        "Rmin": "0",
                        "Rmax": str(N_R - 1),
                        "Tmin": "1",
                        "Tmax": str(N_T),
                        "T0Mode": "fixed",
                        "T0Fixed": "0",
                        "ProfileBasis": "delta",
                        "NProfiles": str(N_V),
                    }
                    mismatches = {
                        key: (meta.get(key), wanted)
                        for key, wanted in expected_meta.items()
                        if meta.get(key) != wanted
                    }
                    if mismatches:
                        raise ValueError(
                            f"{path}:{lineno}: incompatible META "
                            f"{mismatches}"
                        )
                    n_meta += 1

                elif fields[0] == "DATA":
                    if len(fields) != 7:
                        raise ValueError(
                            f"{path}:{lineno}: malformed DATA"
                        )
                    _, c, time, radius, sources, real, imag = fields
                    c, time, radius, sources = map(
                        int, (c, time, radius, sources)
                    )
                    if (
                        c != cfg
                        or sources != N_SRC
                        or not (1 <= time <= N_T)
                        or not (0 <= radius < N_R)
                        or np.isfinite(data[c - 1, time - 1, radius])
                    ):
                        raise ValueError(
                            f"{path}:{lineno}: invalid DATA"
                        )
                    data[c - 1, time - 1, radius] = complex(
                        float(real), float(imag)
                    )
                    n_data += 1

                elif fields[0] == "WDATA":
                    if len(fields) != 9:
                        raise ValueError(
                            f"{path}:{lineno}: malformed WDATA"
                        )
                    (
                        _,
                        c,
                        time,
                        radius,
                        left,
                        right,
                        sources,
                        real,
                        imag,
                    ) = fields
                    c, time, radius, left, right, sources = map(
                        int,
                        (c, time, radius, left, right, sources),
                    )
                    if (
                        c != cfg
                        or sources != N_SRC
                        or not (1 <= time <= N_T)
                        or not (0 <= radius < N_R)
                        or not (
                            0 <= left < N_V and 0 <= right < N_V
                        )
                        or np.isfinite(
                            matrix[
                                c - 1,
                                time - 1,
                                radius,
                                left,
                                right,
                            ]
                        )
                    ):
                        raise ValueError(
                            f"{path}:{lineno}: invalid WDATA"
                        )
                    matrix[
                        c - 1,
                        time - 1,
                        radius,
                        left,
                        right,
                    ] = complex(float(real), float(imag))
                    n_weighted += 1

        if (
            n_meta != 1
            or n_data != N_T * N_R
            or n_weighted != N_T * N_R * N_V**2
        ):
            raise ValueError(
                f"{path}: META={n_meta}, DATA={n_data}, "
                f"WDATA={n_weighted}"
            )

    if not np.all(np.isfinite(data)):
        raise ValueError("non-finite or missing DATA value")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("non-finite or missing WDATA value")

    reconstructed = matrix.sum(axis=(-2, -1))
    if not np.allclose(
        reconstructed,
        data,
        rtol=5.0e-12,
        atol=5.0e-14,
    ):
        difference = np.abs(reconstructed - data)
        raise ValueError(
            "constant-channel reconstruction failed; "
            f"max_abs={difference.max():.3e}"
        )

    return data, matrix


def read_eigenvalues(root: Path, n_cfg: int) -> np.ndarray:
    all_values = []

    for cfg in range(1, n_cfg + 1):
        paths = sorted(
            (root / f"n{cfg}" / "eigenvalues").glob("*.h5")
        )
        if len(paths) != 1:
            raise ValueError(
                f"configuration {cfg}: expected one eigenvalue file "
                f"under {root / f'n{cfg}' / 'eigenvalues'}, "
                f"found {len(paths)}"
            )

        path = paths[0]
        with h5py.File(path, "r") as h5:
            dataset = h5["disteigvals"]
            if dataset.ndim != 2 or dataset.shape[1] < N_V:
                raise ValueError(
                    f"{path}: unexpected disteigvals shape "
                    f"{dataset.shape}"
                )
            values = np.asarray(dataset[:, :N_V], dtype=float)
            times = np.asarray(h5["times"][:], dtype=int)
            cfg_id = int(h5["cnfg_id"][()])

        if values.shape[0] <= N_T:
            raise ValueError(
                f"{path}: need times 0..{N_T}, "
                f"found {values.shape[0]} rows"
            )
        if not np.array_equal(times, np.arange(values.shape[0])):
            raise ValueError(f"{path}: unexpected times")
        if cfg_id != cfg:
            raise ValueError(
                f"{path}: cnfg_id={cfg_id}, expected {cfg}"
            )
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{path}: non-finite eigenvalue")
        if np.any(np.diff(values, axis=1) < -1.0e-13):
            raise ValueError(
                f"{path}: eigenvalues are not ordered"
            )

        all_values.append(values[: N_T + 1])

    return np.stack(all_values)


def project_profiles(
    data: np.ndarray,
    matrix: np.ndarray,
    eigenvalues: np.ndarray,
) -> tuple[list[str], np.ndarray]:
    scaled = (
        eigenvalues[:, None, :, :]
        / SIGMAS[None, :, None, None]
    )
    weights = np.exp(-0.5 * scaled * scaled)

    source = weights[:, :, 0, :]
    sink = weights[:, :, 1 : N_T + 1, :]
    gaussian = np.einsum(
        "cpi,ctrij,cptj->ctrp",
        source,
        matrix,
        sink,
        optimize=True,
    )

    correlators = np.empty(
        (*data.shape, len(SIGMAS) + 1),
        dtype=complex,
    )
    correlators[..., 0] = data
    correlators[..., 1:] = gaussian

    labels = ["constant"] + [
        f"sigma={sigma:.4f}" for sigma in SIGMAS
    ]
    return labels, correlators


def jackknife_error(samples: np.ndarray) -> float:
    n_cfg = len(samples)
    centered = samples - samples.mean()
    return float(
        np.sqrt((n_cfg - 1) / n_cfg * np.sum(centered * centered))
    )


def fit_profile(
    profile: str,
    radius: int,
    correlator: np.ndarray,
) -> FitResult:
    n_cfg = correlator.shape[0]
    real = np.real(correlator)
    total = real.sum(axis=0)
    mean = total / n_cfg
    leave_one_out = (
        total[None, :] - real
    ) / (n_cfg - 1)

    required_effective = np.array(
        sorted(
            {
                EARLY_TIME - 1,
                *(time - 1 for time in PLATEAU_TIMES),
            }
        )
    )
    required_correlators = np.unique(
        np.concatenate(
            (required_effective, required_effective + 1)
        )
    )
    if (
        np.any(mean[required_correlators] <= 0)
        or np.any(
            leave_one_out[:, required_correlators] <= 0
        )
    ):
        raise ValueError(
            f"{profile}, R={radius}: non-positive ensemble or "
            "leave-one-out correlator needed for T/a=2..5"
        )

    effective = np.full(N_T - 1, np.nan)
    valid_mean = (mean[:-1] > 0) & (mean[1:] > 0)
    effective[valid_mean] = np.log(
        mean[:-1][valid_mean] / mean[1:][valid_mean]
    )

    jk_effective = np.full((n_cfg, N_T - 1), np.nan)
    valid_jk = (
        (leave_one_out[:, :-1] > 0)
        & (leave_one_out[:, 1:] > 0)
    )
    jk_effective[valid_jk] = np.log(
        leave_one_out[:, :-1][valid_jk]
        / leave_one_out[:, 1:][valid_jk]
    )

    indices = np.array(PLATEAU_TIMES) - 1
    selected = jk_effective[:, indices]
    centered = selected - selected.mean(axis=0)
    covariance = (
        (n_cfg - 1) / n_cfg * centered.T @ centered
    )
    inverse = np.linalg.pinv(covariance, rcond=1.0e-12)
    ones = np.ones(len(indices))
    denominator = float(ones @ inverse @ ones)
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError(
            f"{profile}, R={radius}: invalid plateau covariance"
        )
    fit_weights = inverse @ ones / denominator

    fit = float(effective[indices] @ fit_weights)
    jk_fit = selected @ fit_weights
    fit_error = jackknife_error(jk_fit)

    residual = effective[indices] - fit
    chi2 = float(residual @ inverse @ residual)
    dof = len(indices) - 1
    if dof != 2:
        raise AssertionError("the p-value formula assumes two dof")
    p_value = float(np.exp(-0.5 * chi2))

    early_index = EARLY_TIME - 1
    early_difference = effective[early_index] - fit
    jk_early_difference = (
        jk_effective[:, early_index] - jk_fit
    )
    early_error = jackknife_error(jk_early_difference)
    early_z = float(
        early_difference / early_error
        if early_error > 0
        else np.inf
    )

    return FitResult(
        profile=profile,
        radius=radius,
        effective=effective,
        fit=fit,
        error=fit_error,
        chi2_dof=chi2 / dof,
        p_value=p_value,
        early_z=early_z,
        jk_fit=jk_fit,
    )


def add_constant_comparisons(
    results: list[FitResult],
) -> None:
    constants = {
        result.radius: result
        for result in results
        if result.profile == "constant"
    }

    for result in results:
        reference = constants[result.radius]
        result.delta_constant = result.fit - reference.fit
        difference = result.jk_fit - reference.jk_fit
        error = jackknife_error(difference)
        result.delta_constant_z = (
            result.delta_constant / error
            if error > 0
            else 0.0
        )


def print_results(
    labels: list[str],
    correlators: np.ndarray,
    results: list[FitResult],
    failures: list[str],
) -> None:
    ensemble_mean = correlators.mean(axis=0)
    real_scale = np.maximum(
        np.abs(np.real(ensemble_mean)),
        1.0e-300,
    )
    imag_ratio = np.abs(np.imag(ensemble_mean)) / real_scale

    print(
        "Profiles: constant and "
        + ", ".join(labels[1:])
    )
    print(
        "Maximum |Im <L>|/|Re <L>|: "
        f"{imag_ratio.max():.3e}"
    )
    print(
        "\nCorrelated effective-energy fits "
        "(plateau T/a=3..5; z2 compares T/a=2 with the plateau)"
    )
    print(
        "profile       R       V2       V3       V4       V5 "
        "    fit      err chi2/dof      p      z2 "
        " dEconst  zconst"
    )

    for result in results:
        values = result.effective
        print(
            f"{result.profile:13s} {result.radius:d} "
            f"{values[1]:+8.5f} {values[2]:+8.5f} "
            f"{values[3]:+8.5f} {values[4]:+8.5f} "
            f"{result.fit:+8.5f} {result.error:8.5f} "
            f"{result.chi2_dof:8.3f} {result.p_value:6.3f} "
            f"{result.early_z:+7.2f} "
            f"{result.delta_constant:+8.5f} "
            f"{result.delta_constant_z:+7.2f}"
        )

    print("\nProfile summary over R/a=1..6")
    print(
        "profile       valid p_min p_median "
        "median|z2| max|zconst| median_relerr"
    )
    for label in labels:
        selected = [
            result
            for result in results
            if result.profile == label
        ]
        if not selected:
            print(
                f"{label:13s} {0:5d} "
                "  nan      nan        nan "
                "        nan           nan"
            )
            continue
        p_values = np.array(
            [result.p_value for result in selected]
        )
        early = np.array(
            [abs(result.early_z) for result in selected]
        )
        constant_z = np.array(
            [abs(result.delta_constant_z) for result in selected]
        )
        relative_error = np.array(
            [
                result.error / abs(result.fit)
                for result in selected
            ]
        )
        print(
            f"{label:13s} {len(selected):5d} "
            f"{p_values.min():5.3f} "
            f"{np.median(p_values):8.3f} "
            f"{np.median(early):10.2f} "
            f"{constant_z.max():11.2f} "
            f"{np.median(relative_error):13.3e}"
        )

    if failures:
        print("\nInvalid profile/radius fits")
        for failure in failures:
            print(f"- {failure}")


def main() -> None:
    args = parse_args()
    outputs = discover_outputs(args.delta_dir, args.jobid)
    data, matrix = read_delta_outputs(outputs)
    eigenvalues = read_eigenvalues(
        args.eigenvalue_root,
        len(outputs),
    )
    labels, correlators = project_profiles(
        data,
        matrix,
        eigenvalues,
    )

    results = []
    failures = []
    for profile_index, label in enumerate(labels):
        for radius in range(1, N_R):
            try:
                result = fit_profile(
                    label,
                    radius,
                    correlators[:, :, radius, profile_index],
                )
            except ValueError as error:
                if label == "constant":
                    raise
                failures.append(str(error))
            else:
                results.append(result)

    add_constant_comparisons(results)

    reconstructed = matrix.sum(axis=(-2, -1))
    absolute = np.abs(reconstructed - data)
    relative = absolute / np.maximum(
        np.abs(data), 1.0e-300
    )
    print(f"Validated configurations: {len(outputs)}")
    print(
        f"Validated DATA/WDATA rows: "
        f"{data.size}/{matrix.size}"
    )
    print(
        "Maximum delta reconstruction residual: "
        f"abs={absolute.max():.3e}, rel={relative.max():.3e}"
    )
    print(
        "Validated eigenvalue shape: "
        f"{eigenvalues.shape}"
    )
    print_results(labels, correlators, results, failures)


if __name__ == "__main__":
    main()
