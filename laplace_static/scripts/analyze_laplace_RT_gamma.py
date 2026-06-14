#!/usr/bin/env python3

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

INPUT_FILE = ROOT / "results/spatial_avg_RT_R0-3_T1-4_n1-100.csv"
JK_FILE = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_jackknife.csv"

OUT_L = ROOT / "results/spatial_avg_RT_R0-3_T1-4_n1-100_gamma.csv"
OUT_V = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_gamma.csv"

FIGDIR = ROOT / "results/figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Keep these consistent with the Wilson-tau Gamma-method workflow.
STAU = 1.5
MIN_W_OPT = 5
MAX_W_DEFAULT = 10


def autocorr_gamma(projected, max_w=None):
    """Estimate Gamma(t) for one mean-zero projected time series."""
    x = np.asarray(projected, dtype=float)
    n = len(x)

    if max_w is None:
        max_w = min(MAX_W_DEFAULT, n // 2 - 1)

    gamma = np.zeros(max_w + 1, dtype=float)

    for lag in range(max_w + 1):
        gamma[lag] = np.dot(x[: n - lag], x[lag:]) / float(n - lag)

    return gamma


def choose_window_wolff_like(gamma, n, stau=STAU, min_w=MIN_W_OPT):
    """Simplified Wolff/UWerr-like automatic summation-window choice."""
    gamma0 = gamma[0]
    max_w = len(gamma) - 1

    if gamma0 <= 0.0 or max_w == 0:
        return 0, "pathological_gamma0"

    rho = gamma / gamma0
    g_int = 0.0
    candidate = None
    min_w_eff = min(min_w, max_w)

    for lag in range(1, max_w + 1):
        g_int += rho[lag]

        if g_int <= 0.0:
            tau_w = np.finfo(float).eps
        else:
            tau_w = stau / math.log((g_int + 1.0) / g_int)

        test_value = (
            math.exp(-lag / tau_w)
            - tau_w / math.sqrt(lag * n)
        )

        if test_value < 0.0 and lag >= min_w_eff:
            candidate = lag
            break

    if candidate is None:
        candidate = min_w_eff
        note = "fallback_fixed_min_window"
    else:
        note = "auto_window"

    # Require a positive autocorrelation sum C(W).
    for window in range(candidate, -1, -1):
        c_window = gamma[0] + 2.0 * np.sum(gamma[1 : window + 1])

        if c_window > 0.0:
            if window != candidate:
                note += "_reduced_for_positive_C"

            return window, note

    return 0, "fallback_W0"


def gamma_analysis(value, projected, max_w=None):
    """Gamma-method error analysis for a projected fluctuation series."""
    x = np.asarray(projected, dtype=float)
    x = x - np.mean(x)
    n = len(x)

    if max_w is None:
        max_w = min(MAX_W_DEFAULT, n // 2 - 1)

    gamma = autocorr_gamma(x, max_w=max_w)
    gamma0 = gamma[0]

    if gamma0 <= 0.0:
        return None

    rho = gamma / gamma0

    tauint_values = []
    error_values = []
    error_uncertainties = []
    tauint_uncertainties = []

    for window in range(max_w + 1):
        c_window = gamma[0] + 2.0 * np.sum(gamma[1 : window + 1])
        tauint = c_window / (2.0 * gamma0)

        if c_window > 0.0:
            error = math.sqrt(c_window / n)

            error_uncertainty = (
                error * math.sqrt((window + 0.5) / n)
            )

            argument = (window - tauint + 0.5) / n

            if argument > 0.0:
                tauint_uncertainty = (
                    2.0 * tauint * math.sqrt(argument)
                )
            else:
                tauint_uncertainty = float("nan")
        else:
            error = float("nan")
            error_uncertainty = float("nan")
            tauint_uncertainty = float("nan")

        tauint_values.append(tauint)
        error_values.append(error)
        error_uncertainties.append(error_uncertainty)
        tauint_uncertainties.append(tauint_uncertainty)

    window_opt, note = choose_window_wolff_like(gamma, n)

    tauint = tauint_values[window_opt]
    neff = n / (2.0 * tauint) if tauint > 0.0 else float("nan")

    return {
        "value": value,
        "error": error_values[window_opt],
        "error_uncertainty": error_uncertainties[window_opt],
        "tauint": tauint,
        "tauint_uncertainty": tauint_uncertainties[window_opt],
        "Wopt": window_opt,
        "Neff": neff,
        "window_note": note,
        "rho": rho,
    }


def main():
    df = pd.read_csv(INPUT_FILE)

    required = {"cfg", "T", "R", "Re", "Im"}
    missing = required.difference(df.columns)

    if missing:
        raise RuntimeError(
            f"Input CSV is missing columns: {sorted(missing)}"
        )

    df = df.sort_values(["cfg", "T", "R"]).reset_index(drop=True)

    cfg_ids = sorted(df["cfg"].unique())
    ncfg = len(cfg_ids)

    print(f"Loaded {len(df)} rows")
    print(f"Configurations: {ncfg}")
    print(f"Configuration range: {cfg_ids[0]}..{cfg_ids[-1]}")

    expected_cfgs = list(range(cfg_ids[0], cfg_ids[-1] + 1))

    if cfg_ids != expected_cfgs:
        raise RuntimeError("Configuration sequence is not contiguous")

    # ------------------------------------------------------------
    # Primary observables: Re L(R,T)
    # ------------------------------------------------------------
    primary_rows = []
    primary_curves = {}

    for T in range(1, 5):
        for R in range(0, 4):
            subset = (
                df[(df["T"] == T) & (df["R"] == R)]
                .sort_values("cfg")
            )

            if subset["cfg"].tolist() != cfg_ids:
                raise RuntimeError(
                    f"Configuration mismatch for T={T}, R={R}"
                )

            values = subset["Re"].to_numpy(dtype=float)
            mean_value = float(np.mean(values))
            naive_error = float(
                np.std(values, ddof=1) / math.sqrt(ncfg)
            )

            result = gamma_analysis(
                mean_value,
                values - mean_value,
            )

            row = {
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "mean_Re": mean_value,
                "err_naive": naive_error,
            }

            if result is None:
                row.update({
                    "err_gamma": float("nan"),
                    "dd_err_gamma": float("nan"),
                    "tauint": float("nan"),
                    "dtauint": float("nan"),
                    "Wopt": -1,
                    "Neff": float("nan"),
                    "window_note": "failed",
                })
            else:
                row.update({
                    "err_gamma": result["error"],
                    "dd_err_gamma": result["error_uncertainty"],
                    "tauint": result["tauint"],
                    "dtauint": result["tauint_uncertainty"],
                    "Wopt": result["Wopt"],
                    "Neff": result["Neff"],
                    "window_note": result["window_note"],
                })

                primary_curves[(T, R)] = result

            primary_rows.append(row)

    primary = pd.DataFrame(primary_rows)
    primary.to_csv(OUT_L, index=False)

    # ------------------------------------------------------------
    # Derived observable:
    #
    # Veff(R,T) = log(<L(R,T)> / <L(R,T+1)>)
    #
    # The projected fluctuation uses the gradient
    #
    # dV/dL_T     =  1/<L_T>
    # dV/dL_{T+1} = -1/<L_{T+1}>
    #
    # preserving covariance between neighbouring T values.
    # ------------------------------------------------------------
    jk = pd.read_csv(JK_FILE) if JK_FILE.exists() else None

    veff_rows = []
    veff_curves = {}

    for T in range(1, 4):
        for R in range(0, 4):
            subset_T = (
                df[(df["T"] == T) & (df["R"] == R)]
                .sort_values("cfg")
            )
            subset_Tp1 = (
                df[(df["T"] == T + 1) & (df["R"] == R)]
                .sort_values("cfg")
            )

            if subset_T["cfg"].tolist() != cfg_ids:
                raise RuntimeError(
                    f"Configuration mismatch for T={T}, R={R}"
                )

            if subset_Tp1["cfg"].tolist() != cfg_ids:
                raise RuntimeError(
                    f"Configuration mismatch for T={T + 1}, R={R}"
                )

            x_T = subset_T["Re"].to_numpy(dtype=float)
            x_Tp1 = subset_Tp1["Re"].to_numpy(dtype=float)

            mean_T = float(np.mean(x_T))
            mean_Tp1 = float(np.mean(x_Tp1))

            if mean_T <= 0.0 or mean_Tp1 <= 0.0:
                raise RuntimeError(
                    f"Non-positive mean for T={T}, R={R}"
                )

            veff = math.log(mean_T / mean_Tp1)

            projected = (
                (x_T - mean_T) / mean_T
                - (x_Tp1 - mean_Tp1) / mean_Tp1
            )

            result = gamma_analysis(veff, projected)

            err_jackknife = float("nan")

            if jk is not None:
                match = jk[
                    (jk["T"] == T)
                    & (jk["R"] == R)
                ]

                if len(match) == 1:
                    err_jackknife = float(
                        match["err_Veff_jackknife"].iloc[0]
                    )

            row = {
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "Veff": veff,
                "err_jackknife": err_jackknife,
            }

            if result is None:
                row.update({
                    "err_gamma": float("nan"),
                    "dd_err_gamma": float("nan"),
                    "tauint": float("nan"),
                    "dtauint": float("nan"),
                    "Wopt": -1,
                    "Neff": float("nan"),
                    "window_note": "failed",
                })
            else:
                row.update({
                    "err_gamma": result["error"],
                    "dd_err_gamma": result["error_uncertainty"],
                    "tauint": result["tauint"],
                    "dtauint": result["tauint_uncertainty"],
                    "Wopt": result["Wopt"],
                    "Neff": result["Neff"],
                    "window_note": result["window_note"],
                })

                veff_curves[(T, R)] = result

            veff_rows.append(row)

    veff_df = pd.DataFrame(veff_rows)

    veff_df["ratio_gamma_over_jackknife"] = (
        veff_df["err_gamma"] / veff_df["err_jackknife"]
    )

    veff_df.to_csv(OUT_V, index=False)

    print()
    print("Re L(R,T) Gamma-method summary:")
    print(primary.to_string(index=False))

    print()
    print("Veff(R,T) Gamma-method summary:")
    print(veff_df.to_string(index=False))

    # ------------------------------------------------------------
    # Plot 1: Gamma-method Veff
    # ------------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        group = veff_df[veff_df["R"] == R]

        plt.errorbar(
            group["T"],
            group["Veff"],
            yerr=group["err_gamma"],
            marker="o",
            capsize=3,
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Laplace effective potential: Gamma-method errors")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_vs_T_gamma.png",
        dpi=200,
    )
    plt.close()

    # ------------------------------------------------------------
    # Plot 2: Gamma/jackknife error ratio
    # ------------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        group = veff_df[veff_df["R"] == R]

        plt.plot(
            group["T"],
            group["ratio_gamma_over_jackknife"],
            marker="o",
            label=f"R/a = {R}",
        )

    plt.axhline(1.0, linestyle="--", linewidth=1)
    plt.xlabel("T/a")
    plt.ylabel("Gamma error / jackknife error")
    plt.title("Autocorrelation inflation of Veff errors")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_gamma_over_jackknife.png",
        dpi=200,
    )
    plt.close()

    # ------------------------------------------------------------
    # Plot 3: tau_int for Re L(R,T)
    # ------------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    for R in range(0, 4):
        group = primary[primary["R"] == R]

        plt.errorbar(
            group["T"],
            group["tauint"],
            yerr=group["dtauint"],
            marker="o",
            capsize=3,
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("tau_int [saved configurations]")
    plt.title("Integrated autocorrelation time of Re L(R,T)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_ReL_tauint_gamma.png",
        dpi=200,
    )
    plt.close()

    # ------------------------------------------------------------
    # Plot 4: tau_int for Veff
    # ------------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        group = veff_df[veff_df["R"] == R]

        plt.errorbar(
            group["T"],
            group["tauint"],
            yerr=group["dtauint"],
            marker="o",
            capsize=3,
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("tau_int [saved configurations]")
    plt.title("Integrated autocorrelation time of Veff")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_tauint_gamma.png",
        dpi=200,
    )
    plt.close()

    # ------------------------------------------------------------
    # Plot 5: representative Veff autocorrelation functions
    # ------------------------------------------------------------
    plt.figure(figsize=(7.0, 4.8))

    representative = [
        (2, 1),
        (2, 2),
        (2, 3),
    ]

    for T, R in representative:
        result = veff_curves.get((T, R))

        if result is None:
            continue

        lags = np.arange(len(result["rho"]))

        plt.plot(
            lags,
            result["rho"],
            marker="o",
            label=f"T/a={T}, R/a={R}",
        )

    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("MC separation [saved configurations]")
    plt.ylabel("rho(lag)")
    plt.title("Normalized autocorrelation of Veff")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_rho_gamma.png",
        dpi=200,
    )
    plt.close()

    print()
    print("Wrote:")
    print(OUT_L)
    print(OUT_V)
    print(FIGDIR / "laplace_Veff_vs_T_gamma.png")
    print(FIGDIR / "laplace_Veff_gamma_over_jackknife.png")
    print(FIGDIR / "laplace_ReL_tauint_gamma.png")
    print(FIGDIR / "laplace_Veff_tauint_gamma.png")
    print(FIGDIR / "laplace_Veff_rho_gamma.png")


if __name__ == "__main__":
    main()
