#!/usr/bin/env python3

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from scipy.stats import chi2
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


SCRIPT_DIR = Path(__file__).resolve().parent

PER_CONFIG_FILE = SCRIPT_DIR / "wilson_tau_per_config.csv"
PYERRORS_VEFF_FILE = SCRIPT_DIR / "wilson_tau_Veff_pyerrors.csv"
JK_VEFF_FILE = SCRIPT_DIR / "wilson_tau_Veff_jackknife.csv"

OUT_TABLE = SCRIPT_DIR / "wilson_tau_plateau_windows.csv"
OUT_BEST = SCRIPT_DIR / "wilson_tau_plateau_best_windows.csv"

# Candidate plateau windows. tau refers to Veff(tau -> tau+1).
# Example: tmin=4, tmax=9 includes Veff(4->5), ..., Veff(9->10).
MIN_POINTS = 3
TMIN_CANDIDATES = range(2, 8)
TMAX_CANDIDATES = range(5, 12)

# Windows to highlight in the plot
HIGHLIGHT_WINDOWS = [(4, 9), (5, 9)]


def jackknife_veff_matrix(W):
    """
    W: DataFrame, rows=configs, columns=tau, values=ReW

    Returns:
      veff_full: DataFrame with tau, tau_next, Veff from full mean values
      jk_values: ndarray shape (Ncfg, Nveff), leave-one-out Veff estimates
    """
    taus = list(W.columns)
    ncfg = len(W)

    full_means = W.mean(axis=0)

    full_rows = []
    jk_values = []

    for tau, tau_next in zip(taus[:-1], taus[1:]):
        delta_tau = tau_next - tau
        full_v = math.log(full_means[tau] / full_means[tau_next]) / delta_tau

        jk_col = []
        for nc in W.index:
            jk_sample = W.drop(index=nc)
            jk_means = jk_sample.mean(axis=0)
            jk_v = math.log(jk_means[tau] / jk_means[tau_next]) / delta_tau
            jk_col.append(jk_v)

        full_rows.append({
            "tau": tau,
            "tau_next": tau_next,
            "Veff": full_v,
        })
        jk_values.append(jk_col)

    veff_full = pd.DataFrame(full_rows)
    jk_values = np.asarray(jk_values, dtype=float).T  # shape: Ncfg x Nveff

    return veff_full, jk_values


def jackknife_cov(jk_values):
    """
    Jackknife covariance for vector-valued observables.

    jk_values shape: Ncfg x Nobs
    """
    ncfg = jk_values.shape[0]
    mean_jk = jk_values.mean(axis=0)
    centered = jk_values - mean_jk
    cov = (ncfg - 1) / ncfg * centered.T @ centered
    return cov


def invert_covariance(cov):
    """
    SVD pseudo-inverse for numerical stability.
    """
    u, s, vh = np.linalg.svd(cov)
    smax = s.max() if len(s) else 0.0

    if smax <= 0:
        raise RuntimeError("Covariance matrix has no positive singular values.")

    cutoff = smax * 1e-12
    s_inv = np.array([1.0 / x if x > cutoff else 0.0 for x in s])

    cov_inv = (vh.T * s_inv) @ u.T
    cond = smax / max(s[s > cutoff].min(), cutoff)

    return cov_inv, cond


def correlated_constant_fit(y, cov, jk_samples_window):
    """
    Correlated constant fit using jackknife covariance.

    y: full-sample Veff values in window
    cov: jackknife covariance matrix in window
    jk_samples_window: Ncfg x Npoints jackknife estimates in same window
    """
    npts = len(y)
    one = np.ones(npts)

    cov_inv, cond = invert_covariance(cov)

    denom = one @ cov_inv @ one
    weights = cov_inv @ one / denom

    v_fit = weights @ y

    residual = y - v_fit
    chi2_val = float(residual @ cov_inv @ residual)
    dof = npts - 1
    p_value = chi2.sf(chi2_val, dof) if HAVE_SCIPY and dof > 0 else np.nan

    # Formula error from covariance
    err_formula = math.sqrt(1.0 / denom)

    # Jackknife error of fitted plateau value, using fixed correlated-fit weights
    jk_plateau = jk_samples_window @ weights
    mean_jk = jk_plateau.mean()
    ncfg = len(jk_plateau)
    err_jk = math.sqrt((ncfg - 1) / ncfg * np.sum((jk_plateau - mean_jk) ** 2))

    return {
        "V_corr": v_fit,
        "err_corr_formula": err_formula,
        "err_corr_jackknife": err_jk,
        "chi2_corr": chi2_val,
        "dof_corr": dof,
        "chi2_dof_corr": chi2_val / dof if dof > 0 else np.nan,
        "p_value_corr": p_value,
        "cov_condition": cond,
    }


def uncorrelated_constant_fit(y, err):
    """
    Uncorrelated weighted constant fit using diagonal errors.
    """
    w = 1.0 / np.asarray(err) ** 2
    v_fit = np.sum(w * y) / np.sum(w)
    err_fit = math.sqrt(1.0 / np.sum(w))

    residual = y - v_fit
    chi2_val = float(np.sum((residual / err) ** 2))
    dof = len(y) - 1
    p_value = chi2.sf(chi2_val, dof) if HAVE_SCIPY and dof > 0 else np.nan

    return {
        "V_uncorr": v_fit,
        "err_uncorr": err_fit,
        "chi2_uncorr": chi2_val,
        "dof_uncorr": dof,
        "chi2_dof_uncorr": chi2_val / dof if dof > 0 else np.nan,
        "p_value_uncorr": p_value,
    }


def main():
    df = pd.read_csv(PER_CONFIG_FILE)
    W = df.pivot(index="nc", columns="tau", values="ReW").sort_index()

    veff_full, jk_values = jackknife_veff_matrix(W)
    taus_eff = list(veff_full["tau"])
    cov_all = jackknife_cov(jk_values)

    py = pd.read_csv(PYERRORS_VEFF_FILE)
    jk = pd.read_csv(JK_VEFF_FILE)

    veff = veff_full.merge(
        py[["tau", "tau_next", "err_pyerrors"]],
        on=["tau", "tau_next"],
        how="left",
    ).merge(
        jk[["tau", "tau_next", "err_Veff_jackknife"]],
        on=["tau", "tau_next"],
        how="left",
    )

    rows = []

    for tmin in TMIN_CANDIDATES:
        for tmax in TMAX_CANDIDATES:
            if tmax < tmin:
                continue

            mask = (veff["tau"] >= tmin) & (veff["tau"] <= tmax)
            idx = np.where(mask.to_numpy())[0]

            if len(idx) < MIN_POINTS:
                continue

            y = veff.loc[mask, "Veff"].to_numpy(dtype=float)
            err_py = veff.loc[mask, "err_pyerrors"].to_numpy(dtype=float)

            cov_win = cov_all[np.ix_(idx, idx)]
            jk_win = jk_values[:, idx]

            try:
                corr = correlated_constant_fit(y, cov_win, jk_win)
                corr_ok = True
            except Exception as exc:
                corr = {
                    "V_corr": np.nan,
                    "err_corr_formula": np.nan,
                    "err_corr_jackknife": np.nan,
                    "chi2_corr": np.nan,
                    "dof_corr": len(idx) - 1,
                    "chi2_dof_corr": np.nan,
                    "p_value_corr": np.nan,
                    "cov_condition": np.nan,
                }
                corr_ok = False

            uncorr = uncorrelated_constant_fit(y, err_py)

            rows.append({
                "tmin": tmin,
                "tmax": tmax,
                "npts": len(idx),
                "corr_ok": corr_ok,
                **corr,
                **uncorr,
            })

    table = pd.DataFrame(rows)

    # A useful ranking:
    # prefer acceptable p-value, moderate condition number, and more points.
    table["rank_score"] = (
        table["p_value_corr"].fillna(0.0)
        - 0.02 * (table["chi2_dof_corr"].fillna(99.0) - 1.0).abs()
        + 0.005 * table["npts"]
    )

    table = table.sort_values(
        by=["corr_ok", "rank_score", "npts"],
        ascending=[False, False, False],
    )

    table.to_csv(OUT_TABLE, index=False)

    best = table.head(12).copy()
    best.to_csv(OUT_BEST, index=False)

    print("Best candidate plateau windows:")
    cols = [
        "tmin", "tmax", "npts",
        "V_corr", "err_corr_jackknife",
        "chi2_dof_corr", "p_value_corr",
        "V_uncorr", "err_uncorr",
    ]
    print(best[cols].to_string(index=False))

    print()
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_BEST}")

    # ------------------------------------------------------------
    # Plot Veff with highlighted plateau windows
    # ------------------------------------------------------------
    plt.figure(figsize=(8, 5.5))

    plt.errorbar(
        veff["tau"],
        veff["Veff"],
        yerr=veff["err_pyerrors"],
        fmt="o",
        capsize=4,
        label=r"$V_{\mathrm{eff}}$ pyerrors",
    )

    colors = ["C1", "C2", "C3", "C4"]

    for i, (tmin, tmax) in enumerate(HIGHLIGHT_WINDOWS):
        match = table[(table["tmin"] == tmin) & (table["tmax"] == tmax)]
        if len(match) != 1:
            continue

        row = match.iloc[0]
        v = row["V_corr"]
        e = row["err_corr_jackknife"]
        color = colors[i % len(colors)]

        plt.axhspan(
            v - e,
            v + e,
            alpha=0.18,
            color=color,
            label=fr"corr. fit $\tau={tmin}\ldots{tmax}$",
        )
        plt.hlines(
            v,
            tmin - 0.15,
            tmax + 0.15,
            colors=color,
            linewidth=2,
        )

    plt.xlabel(r"$\tau/a$")
    plt.ylabel(r"$V_{\mathrm{eff}}(a,\tau)$")
    plt.title(r"Plateau fit for fixed spatial extent $a$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_plateau_fit.png", dpi=200)

    # ------------------------------------------------------------
    # Plot plateau value versus fit window
    # ------------------------------------------------------------
    plot_table = table[table["corr_ok"]].copy()
    plot_table = plot_table.sort_values(["tmin", "tmax"])

    labels = [f"{int(r.tmin)}-{int(r.tmax)}" for r in plot_table.itertuples()]
    x = np.arange(len(plot_table))

    plt.figure(figsize=(max(9, 0.35 * len(plot_table)), 5.5))
    plt.errorbar(
        x,
        plot_table["V_corr"],
        yerr=plot_table["err_corr_jackknife"],
        fmt="o",
        capsize=3,
    )
    plt.xticks(x, labels, rotation=60, ha="right")
    plt.xlabel(r"fit window $\tau_{\min}$-$\tau_{\max}$")
    plt.ylabel(r"$V(a)$")
    plt.title(r"Plateau stability over fit windows")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_plateau_window_scan.png", dpi=200)

    print()
    print("Wrote plots:")
    print(SCRIPT_DIR / "wilson_tau_plateau_fit.png")
    print(SCRIPT_DIR / "wilson_tau_plateau_window_scan.png")


if __name__ == "__main__":
    main()