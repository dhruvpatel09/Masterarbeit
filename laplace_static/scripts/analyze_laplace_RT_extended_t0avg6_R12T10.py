#!/usr/bin/env python3

from pathlib import Path

import math

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyerrors as pe


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"

FIGDIR.mkdir(parents=True, exist_ok=True)

INPUT = (
    RESULTS /
    "spatial_avg_RT_R0-12_T1-10_t0avg6_n1-100.csv"
)

OUT_STATUS = (
    RESULTS /
    "laplace_extended_R0-12_T1-10_t0avg6_n1-100_status.csv"
)

OUT_VEFF = (
    RESULTS /
    "laplace_Veff_R0-12_T1-9_t0avg6_n1-100_pyerrors.csv"
)

EXPECTED_CFGS = list(range(1, 101))
EXPECTED_R = list(range(0, 13))
EXPECTED_T = list(range(1, 11))
EXPECTED_NSRC = 6 * 24 ** 3
ENSEMBLE_NAME = "Em1p4_R12T10"
GAMMA_S = 2.0


def validate_input(df):
    required = {
        "cfg",
        "T",
        "R",
        "Nsrc",
        "Re",
        "Im",
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"Missing columns: {sorted(missing)}"
        )

    if len(df) != 13000:
        raise RuntimeError(
            f"Expected 13000 rows, found {len(df)}"
        )

    if df.duplicated(["cfg", "T", "R"]).any():
        raise RuntimeError(
            "Duplicate cfg,T,R entries found"
        )

    if sorted(df["cfg"].unique()) != EXPECTED_CFGS:
        raise RuntimeError(
            "Unexpected cfg range"
        )

    if sorted(df["R"].unique()) != EXPECTED_R:
        raise RuntimeError(
            "Unexpected R range"
        )

    if sorted(df["T"].unique()) != EXPECTED_T:
        raise RuntimeError(
            "Unexpected T range"
        )

    if not (df["Nsrc"] == EXPECTED_NSRC).all():
        raise RuntimeError(
            "Unexpected Nsrc value"
        )


def get_obs_value(obs):
    for attr in ["value", "val"]:
        if hasattr(obs, attr):
            return float(getattr(obs, attr))

    raise RuntimeError(
        "Could not extract value from pyerrors Obs"
    )


def get_obs_error(obs):
    for attr in ["dvalue", "dval", "err"]:
        if hasattr(obs, attr):
            return float(getattr(obs, attr))

    raise RuntimeError(
        "Could not extract error from pyerrors Obs"
    )


def make_obs(samples):
    return pe.Obs(
        [np.asarray(samples, dtype=float)],
        [ENSEMBLE_NAME],
    )


def pyerrors_log_ratio(samples_t, samples_tp1):
    samples_t = np.asarray(samples_t, dtype=float)
    samples_tp1 = np.asarray(samples_tp1, dtype=float)

    ncfg = samples_t.size

    if ncfg != samples_tp1.size:
        raise RuntimeError(
            "Mismatched sample lengths"
        )

    central_t = float(np.mean(samples_t))
    central_tp1 = float(np.mean(samples_tp1))

    if central_t <= 0.0 or central_tp1 <= 0.0:
        return {
            "Veff": np.nan,
            "err_pyerrors": np.nan,
            "relative_err_percent": np.nan,
            "status": "nonpositive_central_mean",
        }

    obs_t = make_obs(samples_t)
    obs_tp1 = make_obs(samples_tp1)

    attempts = []

    try:
        veff_obs = pe.log(obs_t / obs_tp1)
    except Exception as exc:
        attempts.append(
            f"pe.log failed: {exc}"
        )

        try:
            veff_obs = np.log(obs_t / obs_tp1)
        except Exception as exc:
            attempts.append(
                f"np.log failed: {exc}"
            )

            if hasattr(pe, "derived_observable"):
                try:
                    veff_obs = pe.derived_observable(
                        lambda x: np.log(x[0] / x[1]),
                        [obs_t, obs_tp1],
                    )
                except Exception as exc:
                    attempts.append(
                        f"pe.derived_observable failed: {exc}"
                    )
                    raise RuntimeError(
                        "Could not construct pyerrors log-ratio observable:\n"
                        + "\n".join(attempts)
                    ) from exc
            else:
                raise RuntimeError(
                    "Could not construct pyerrors log-ratio observable:\n"
                    + "\n".join(attempts)
                )

    veff_obs.gamma_method(S=GAMMA_S)

    value = get_obs_value(veff_obs)
    error = get_obs_error(veff_obs)

    if not np.isfinite(value) or not np.isfinite(error):
        return {
            "Veff": value,
            "err_pyerrors": error,
            "relative_err_percent": np.nan,
            "status": "pyerrors_nan",
        }

    relative = (
        100.0 * error / abs(value)
        if value != 0.0
        else np.nan
    )

    return {
        "Veff": value,
        "err_pyerrors": error,
        "relative_err_percent": relative,
        "status": "ok",
    }


def jackknife_log_ratio(samples_t, samples_tp1):
    samples_t = np.asarray(samples_t, dtype=float)
    samples_tp1 = np.asarray(samples_tp1, dtype=float)

    ncfg = samples_t.size

    if ncfg != samples_tp1.size:
        raise RuntimeError(
            "Mismatched sample lengths"
        )

    total_t = float(np.sum(samples_t))
    total_tp1 = float(np.sum(samples_tp1))

    central_t = total_t / ncfg
    central_tp1 = total_tp1 / ncfg

    if central_t <= 0.0 or central_tp1 <= 0.0:
        return {
            "err_jackknife": np.nan,
            "relative_err_jackknife_percent": np.nan,
            "n_valid_jackknife": 0,
            "jackknife_status": "nonpositive_central_mean",
        }

    jk_values = []

    for i in range(ncfg):
        mean_t = (total_t - samples_t[i]) / (ncfg - 1)
        mean_tp1 = (
            total_tp1 - samples_tp1[i]
        ) / (ncfg - 1)

        if mean_t <= 0.0 or mean_tp1 <= 0.0:
            continue

        jk_values.append(
            math.log(mean_t / mean_tp1)
        )

    jk_values = np.asarray(jk_values, dtype=float)

    if jk_values.size < ncfg:
        return {
            "err_jackknife": np.nan,
            "relative_err_jackknife_percent": np.nan,
            "n_valid_jackknife": int(jk_values.size),
            "jackknife_status": "invalid_jackknife_samples",
        }

    central = math.log(
        central_t / central_tp1
    )

    err = math.sqrt(
        (ncfg - 1) / ncfg
        * float(
            np.sum(
                (jk_values - jk_values.mean()) ** 2
            )
        )
    )

    relative = (
        100.0 * err / abs(central)
        if central != 0.0
        else np.nan
    )

    return {
        "err_jackknife": err,
        "relative_err_jackknife_percent": relative,
        "n_valid_jackknife": int(jk_values.size),
        "jackknife_status": "ok",
    }


def add_reliability_flags(row):
    if row["status"] != "ok":
        return "invalid"

    if not np.isfinite(row["relative_err_percent"]):
        return "invalid"

    if row["relative_err_percent"] < 10.0:
        return "good"

    if row["relative_err_percent"] < 25.0:
        return "marginal"

    return "diagnostic_only"


def plot_veff_vs_t(veff):
    plt.figure(figsize=(8.2, 5.4))

    for R in range(1, 13):
        group = (
            veff[
                (veff["R"] == R)
                & (veff["status"] == "ok")
            ]
            .sort_values("T")
        )

        if group.empty:
            continue

        plt.errorbar(
            group["T"],
            group["Veff"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=3,
            linewidth=1.0,
            label=f"R/a={R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Extended Laplace effective potential, pyerrors/Gamma method")
    plt.xticks(range(1, 10))
    plt.grid(True, alpha=0.3)
    plt.legend(ncol=3, fontsize=8)
    plt.tight_layout()
    plt.savefig(
        FIGDIR /
        "laplace_Veff_vs_T_pyerrors_t0avg6_R12T10.png",
        dpi=200,
    )
    plt.close()


def plot_selected_t_vs_r(veff):
    plt.figure(figsize=(8.2, 5.4))

    for T in [2, 3, 4, 5]:
        group = (
            veff[
                (veff["T"] == T)
                & (veff["R"] > 0)
                & (veff["status"] == "ok")
            ]
            .sort_values("R")
        )

        if group.empty:
            continue

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
        alpha=0.6,
        label="half lattice: R/a=12",
    )

    plt.xlabel("R/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Extended early-time potential versus R")
    plt.xticks(range(1, 13))
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR /
        "laplace_Veff_vs_R_selected_T_pyerrors_t0avg6_R12T10.png",
        dpi=200,
    )
    plt.close()


def plot_relative_error(veff):
    plt.figure(figsize=(8.2, 5.4))

    for R in range(1, 13):
        group = (
            veff[
                (veff["R"] == R)
                & np.isfinite(
                    veff["relative_err_percent"]
                )
            ]
            .sort_values("T")
        )

        if group.empty:
            continue

        plt.plot(
            group["T"],
            group["relative_err_percent"],
            marker="o",
            linewidth=1.0,
            label=f"R/a={R}",
        )

    plt.axhline(
        10.0,
        linestyle="--",
        linewidth=1.0,
        alpha=0.7,
        label="10%",
    )

    plt.axhline(
        25.0,
        linestyle=":",
        linewidth=1.0,
        alpha=0.7,
        label="25%",
    )

    plt.yscale("log")
    plt.xlabel("T/a")
    plt.ylabel("Relative uncertainty [%]")
    plt.title("Relative uncertainty of extended Veff, pyerrors/Gamma method")
    plt.xticks(range(1, 10))
    plt.grid(True, alpha=0.3, which="both")
    plt.legend(ncol=3, fontsize=8)
    plt.tight_layout()
    plt.savefig(
        FIGDIR /
        "laplace_Veff_relative_error_pyerrors_t0avg6_R12T10.png",
        dpi=200,
    )
    plt.close()


def plot_nonpositive_counts(status):
    pivot = (
        status.pivot(
            index="R",
            columns="T",
            values="nonpositive_configs",
        )
        .sort_index(ascending=False)
    )

    plt.figure(figsize=(8.4, 5.8))
    image = plt.imshow(
        pivot.to_numpy(),
        aspect="auto",
        interpolation="nearest",
    )

    plt.colorbar(
        image,
        label="nonpositive configs out of 100",
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
    plt.title("Nonpositive correlator samples in extended scan")

    plt.tight_layout()
    plt.savefig(
        FIGDIR /
        "laplace_nonpositive_counts_t0avg6_R12T10.png",
        dpi=200,
    )
    plt.close()


def main():
    df = pd.read_csv(INPUT)

    df = df.sort_values(
        ["cfg", "T", "R"]
    ).reset_index(drop=True)

    validate_input(df)

    status_rows = []

    for R in EXPECTED_R:
        for T in EXPECTED_T:
            group = (
                df[
                    (df["R"] == R)
                    & (df["T"] == T)
                ]
                .sort_values("cfg")
            )

            samples = group["Re"].to_numpy(dtype=float)

            mean = float(np.mean(samples))
            sem = float(
                np.std(samples, ddof=1) /
                math.sqrt(samples.size)
            )

            rel = (
                100.0 * abs(sem) / abs(mean)
                if mean != 0.0
                else np.nan
            )

            status_rows.append({
                "R": R,
                "T": T,
                "Ncfg": samples.size,
                "mean_Re": mean,
                "sem_Re": sem,
                "relative_err_Re_percent": rel,
                "nonpositive_configs": int(
                    np.count_nonzero(samples <= 0.0)
                ),
            })

    status = pd.DataFrame(status_rows)

    status.to_csv(
        OUT_STATUS,
        index=False,
        lineterminator="\n",
    )

    veff_rows = []

    for R in EXPECTED_R:
        for T in range(1, 10):
            samples_t = (
                df[
                    (df["R"] == R)
                    & (df["T"] == T)
                ]
                .sort_values("cfg")["Re"]
                .to_numpy(dtype=float)
            )

            samples_tp1 = (
                df[
                    (df["R"] == R)
                    & (df["T"] == T + 1)
                ]
                .sort_values("cfg")["Re"]
                .to_numpy(dtype=float)
            )

            gamma_result = pyerrors_log_ratio(
                samples_t,
                samples_tp1,
            )

            jack_result = jackknife_log_ratio(
                samples_t,
                samples_tp1,
            )

            veff_rows.append({
                "R": R,
                "T": T,
                "Ncfg": samples_t.size,
                **gamma_result,
                **jack_result,
            })

    veff = pd.DataFrame(veff_rows)

    veff["reliability"] = veff.apply(
        add_reliability_flags,
        axis=1,
    )

    veff.to_csv(
        OUT_VEFF,
        index=False,
        lineterminator="\n",
    )

    plot_veff_vs_t(veff)
    plot_selected_t_vs_r(veff)
    plot_relative_error(veff)
    plot_nonpositive_counts(status)

    print("Extended R12T10 diagnostic summary")
    print("Primary errors: pyerrors/Gamma method")
    print(f"Gamma-method S parameter: {GAMMA_S}")
    print()

    print("Veff status counts:")
    print(
        veff["status"]
        .value_counts(dropna=False)
        .to_string()
    )
    print()

    print("Reliability counts:")
    print(
        veff["reliability"]
        .value_counts(dropna=False)
        .to_string()
    )
    print()

    print("Jackknife cross-check status counts:")
    print(
        veff["jackknife_status"]
        .value_counts(dropna=False)
        .to_string()
    )
    print()

    print("Large-R Veff table:")
    print(
        veff[
            (veff["R"].between(7, 12))
            & (veff["T"].between(2, 8))
        ][[
            "R",
            "T",
            "Veff",
            "err_pyerrors",
            "relative_err_percent",
            "err_jackknife",
            "status",
            "jackknife_status",
            "reliability",
        ]]
        .to_string(
            index=False,
            float_format=lambda value: f"{value:.6g}",
        )
    )
    print()

    print("Wrote:")
    print(OUT_STATUS)
    print(OUT_VEFF)
    print(
        FIGDIR /
        "laplace_Veff_vs_T_pyerrors_t0avg6_R12T10.png"
    )
    print(
        FIGDIR /
        "laplace_Veff_vs_R_selected_T_pyerrors_t0avg6_R12T10.png"
    )
    print(
        FIGDIR /
        "laplace_Veff_relative_error_pyerrors_t0avg6_R12T10.png"
    )
    print(
        FIGDIR /
        "laplace_nonpositive_counts_t0avg6_R12T10.png"
    )


if __name__ == "__main__":
    main()
