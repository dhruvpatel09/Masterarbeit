#!/usr/bin/env python3

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
PER_CONFIG_FILE = SCRIPT_DIR / "wilson_tau_per_config.csv"
JK_VEFF_FILE = SCRIPT_DIR / "wilson_tau_Veff_jackknife.csv"

STAU = 1.5
MIN_W_OPT = 5


def autocorr_gamma(projected, max_w=None):
    """
    Estimate Gamma(t) for one projected fluctuation time series.
    projected should already have mean zero.
    """
    x = np.asarray(projected, dtype=float)
    n = len(x)

    if max_w is None:
        max_w = min(10, n // 2 - 1)

    gamma = np.zeros(max_w + 1, dtype=float)

    for t in range(max_w + 1):
        gamma[t] = np.dot(x[: n - t], x[t:]) / float(n - t)

    return gamma


def choose_window_wolff_like(gamma, n, stau=STAU, min_w=MIN_W_OPT):
    """
    Simplified Wolff/UWerr-like automatic window choice.

    For N=30 this should be viewed as diagnostic. We enforce a minimum
    window when possible, but fall back safely if C(W) becomes pathological.
    """
    gamma0 = gamma[0]
    max_w = len(gamma) - 1

    if gamma0 <= 0 or max_w == 0:
        return 0, "pathological_gamma0"

    rho = gamma / gamma0
    g_int = 0.0
    candidate = None

    min_w_eff = min(min_w, max_w)

    for t in range(1, max_w + 1):
        g_int += rho[t]

        if g_int <= 0:
            tau_w = np.finfo(float).eps
        else:
            tau_w = stau / math.log((g_int + 1.0) / g_int)

        g_w = math.exp(-t / tau_w) - tau_w / math.sqrt(t * n)

        if g_w < 0 and t >= min_w_eff:
            candidate = t
            break

    if candidate is None:
        candidate = min_w_eff
        note = "fallback_fixed_min_window"
    else:
        note = "auto_window"

    # Ensure positive C(W). If not, reduce W until positive.
    for w in range(candidate, -1, -1):
        c_w = gamma[0] + 2.0 * np.sum(gamma[1 : w + 1])
        if c_w > 0:
            if w != candidate:
                note += "_reduced_for_positive_C"
            return w, note

    return 0, "fallback_W0"


def gamma_analysis(value, projected, max_w=None):
    x = np.asarray(projected, dtype=float)
    x = x - np.mean(x)
    n = len(x)

    if max_w is None:
        max_w = min(10, n // 2 - 1)

    gamma = autocorr_gamma(x, max_w=max_w)
    gamma0 = gamma[0]

    if gamma0 <= 0:
        return None

    rho = gamma / gamma0

    c_vals = []
    tauint_vals = []
    err_vals = []
    dd_err_vals = []
    dtauint_vals = []

    for w in range(max_w + 1):
        c_w = gamma[0] + 2.0 * np.sum(gamma[1 : w + 1])
        c_vals.append(c_w)

        tauint_w = c_w / (2.0 * gamma0)
        tauint_vals.append(tauint_w)

        if c_w > 0:
            err_w = math.sqrt(c_w / n)
            dd_err_w = err_w * math.sqrt((w + 0.5) / n)
            arg = (w - tauint_w + 0.5) / n
            dtauint_w = tauint_w * 2.0 * math.sqrt(arg) if arg > 0 else float("nan")
        else:
            err_w = float("nan")
            dd_err_w = float("nan")
            dtauint_w = float("nan")

        err_vals.append(err_w)
        dd_err_vals.append(dd_err_w)
        dtauint_vals.append(dtauint_w)

    wopt, note = choose_window_wolff_like(gamma, n)

    dvalue = err_vals[wopt]
    ddvalue = dd_err_vals[wopt]
    tauint = tauint_vals[wopt]
    dtauint = dtauint_vals[wopt]
    neff = n / (2.0 * tauint) if tauint > 0 else float("nan")

    return {
        "value": value,
        "dvalue_gamma": dvalue,
        "ddvalue_gamma": ddvalue,
        "tauint": tauint,
        "dtauint": dtauint,
        "Wopt": wopt,
        "window_note": note,
        "Ncfg": n,
        "Neff": neff,
        "gamma": gamma,
        "rho": rho,
        "tauint_vs_W": np.asarray(tauint_vals),
        "err_vs_W": np.asarray(err_vals),
        "max_w": max_w,
    }


def main():
    df = pd.read_csv(PER_CONFIG_FILE)

    W = df.pivot(index="nc", columns="tau", values="ReW").sort_index()
    taus = list(W.columns)
    ncfg = len(W)

    print(f"Loaded {ncfg} configurations and taus: {taus}")

    primary_rows = []
    primary_curves = {}

    # ------------------------------------------------------------
    # Primary observable: <Re W(tau)>
    # ------------------------------------------------------------
    for tau in taus:
        x = W[tau].to_numpy(dtype=float)
        value = float(np.mean(x))
        naive = float(np.std(x, ddof=1) / math.sqrt(ncfg))

        projected = x - value
        res = gamma_analysis(value, projected)

        row = {
            "tau": tau,
            "Ncfg": ncfg,
            "mean_ReW": value,
            "err_naive": naive,
        }

        if res is not None:
            row.update({
                "err_gamma": res["dvalue_gamma"],
                "dd_err_gamma": res["ddvalue_gamma"],
                "tauint": res["tauint"],
                "dtauint": res["dtauint"],
                "Wopt": res["Wopt"],
                "Neff": res["Neff"],
                "window_note": res["window_note"],
            })
            primary_curves[tau] = res
        else:
            row.update({
                "err_gamma": float("nan"),
                "dd_err_gamma": float("nan"),
                "tauint": float("nan"),
                "dtauint": float("nan"),
                "Wopt": -1,
                "Neff": float("nan"),
                "window_note": "failed",
            })

        primary_rows.append(row)

    df_primary = pd.DataFrame(primary_rows)
    df_primary.to_csv(SCRIPT_DIR / "wilson_tau_ReW_gamma_method.csv", index=False)

    # ------------------------------------------------------------
    # Derived observable:
    # Veff(tau) = log(mean W_tau / mean W_tau_next) / delta_tau
    #
    # Use projected Gamma-method fluctuations:
    # grad = [1/W_tau, -1/W_tau_next] / delta_tau
    # ------------------------------------------------------------
    veff_rows = []
    veff_curves = {}

    jk_df = None
    if JK_VEFF_FILE.exists():
        jk_df = pd.read_csv(JK_VEFF_FILE)

    for tau, tau_next in zip(taus[:-1], taus[1:]):
        data = W[[tau, tau_next]].to_numpy(dtype=float)
        means = np.mean(data, axis=0)

        delta_tau = float(tau_next - tau)
        value = math.log(means[0] / means[1]) / delta_tau

        grad = np.array([1.0 / means[0], -1.0 / means[1]]) / delta_tau
        projected = np.dot(data - means, grad)

        res = gamma_analysis(value, projected)

        err_jk = float("nan")
        if jk_df is not None:
            match = jk_df[(jk_df["tau"] == tau) & (jk_df["tau_next"] == tau_next)]
            if len(match) == 1:
                err_jk = float(match["err_Veff_jackknife"].iloc[0])

        row = {
            "tau": tau,
            "tau_next": tau_next,
            "Ncfg": ncfg,
            "Veff": value,
            "err_jackknife": err_jk,
        }

        if res is not None:
            row.update({
                "err_gamma": res["dvalue_gamma"],
                "dd_err_gamma": res["ddvalue_gamma"],
                "tauint": res["tauint"],
                "dtauint": res["dtauint"],
                "Wopt": res["Wopt"],
                "Neff": res["Neff"],
                "window_note": res["window_note"],
            })
            veff_curves[tau] = res
        else:
            row.update({
                "err_gamma": float("nan"),
                "dd_err_gamma": float("nan"),
                "tauint": float("nan"),
                "dtauint": float("nan"),
                "Wopt": -1,
                "Neff": float("nan"),
                "window_note": "failed",
            })

        veff_rows.append(row)

    df_veff = pd.DataFrame(veff_rows)
    df_veff.to_csv(SCRIPT_DIR / "wilson_tau_Veff_gamma_method.csv", index=False)

    print()
    print("Primary <ReW> Gamma-method summary:")
    print(df_primary.to_string(index=False))

    print()
    print("Veff Gamma-method summary:")
    print(df_veff.to_string(index=False))

    # ------------------------------------------------------------
    # Plot: Veff with jackknife and Gamma-method errors
    # ------------------------------------------------------------
    plt.figure(figsize=(7, 5))

    if "err_jackknife" in df_veff:
        plt.errorbar(
            df_veff["tau"] - 0.04,
            df_veff["Veff"],
            yerr=df_veff["err_jackknife"],
            fmt="o",
            capsize=4,
            label="jackknife",
        )

    plt.errorbar(
        df_veff["tau"] + 0.04,
        df_veff["Veff"],
        yerr=df_veff["err_gamma"],
        fmt="s",
        capsize=4,
        label="Gamma method",
    )

    plt.xlabel(r"$\tau/a$")
    plt.ylabel(r"$V_{\mathrm{eff}}(a,\tau)$")
    plt.title(r"Effective potential: jackknife vs Gamma-method errors")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_Veff_gamma_vs_jackknife.png", dpi=200)

    # ------------------------------------------------------------
    # Plot: tau_int for primary W(tau)
    # ------------------------------------------------------------
    plt.figure(figsize=(7, 5))
    plt.errorbar(
        df_primary["tau"],
        df_primary["tauint"],
        yerr=df_primary["dtauint"],
        fmt="o-",
        capsize=4,
    )
    plt.xlabel(r"$\tau/a$")
    plt.ylabel(r"$\tau_{\mathrm{int}}$")
    plt.title(r"Integrated autocorrelation time for $\mathrm{Re}\,W(a,\tau)$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_ReW_tauint_gamma.png", dpi=200)

    # ------------------------------------------------------------
    # Plot: tau_int for Veff
    # ------------------------------------------------------------
    plt.figure(figsize=(7, 5))
    plt.errorbar(
        df_veff["tau"],
        df_veff["tauint"],
        yerr=df_veff["dtauint"],
        fmt="o-",
        capsize=4,
    )
    plt.xlabel(r"$\tau/a$")
    plt.ylabel(r"$\tau_{\mathrm{int}}$")
    plt.title(r"Integrated autocorrelation time for $V_{\mathrm{eff}}$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_Veff_tauint_gamma.png", dpi=200)

    # ------------------------------------------------------------
    # Plot representative rho(t) curves
    # ------------------------------------------------------------
    selected_primary = [1, 6, 12]
    plt.figure(figsize=(7, 5))
    for tau in selected_primary:
        if tau in primary_curves:
            res = primary_curves[tau]
            t = np.arange(len(res["rho"]))
            plt.plot(t, res["rho"], "o-", label=fr"$\tau={tau}$")
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("MC separation")
    plt.ylabel(r"$\rho(t)$")
    plt.title(r"Normalized autocorrelation: $\mathrm{Re}\,W(a,\tau)$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_ReW_rho_gamma.png", dpi=200)

    selected_veff = [1, 4, 8, 10]
    plt.figure(figsize=(7, 5))
    for tau in selected_veff:
        if tau in veff_curves:
            res = veff_curves[tau]
            t = np.arange(len(res["rho"]))
            plt.plot(t, res["rho"], "o-", label=fr"$\tau={tau}$")
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("MC separation")
    plt.ylabel(r"$\rho(t)$")
    plt.title(r"Normalized autocorrelation: $V_{\mathrm{eff}}$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_Veff_rho_gamma.png", dpi=200)

    print()
    print("Wrote:")
    print(SCRIPT_DIR / "wilson_tau_ReW_gamma_method.csv")
    print(SCRIPT_DIR / "wilson_tau_Veff_gamma_method.csv")
    print(SCRIPT_DIR / "wilson_tau_Veff_gamma_vs_jackknife.png")
    print(SCRIPT_DIR / "wilson_tau_ReW_tauint_gamma.png")
    print(SCRIPT_DIR / "wilson_tau_Veff_tauint_gamma.png")
    print(SCRIPT_DIR / "wilson_tau_ReW_rho_gamma.png")
    print(SCRIPT_DIR / "wilson_tau_Veff_rho_gamma.png")


if __name__ == "__main__":
    main()