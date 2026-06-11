#!/usr/bin/env python3

import math
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

OUT_DIR = Path("analysis_wilson_tau")
INFILE = OUT_DIR / "wilson_tau_per_config.csv"

df = pd.read_csv(INFILE)

# Pivot to matrix form:
# rows = configuration number, columns = tau, values = ReW
W = df.pivot(index="nc", columns="tau", values="ReW").sort_index()
taus = list(W.columns)
N = len(W)

if N < 2:
    raise RuntimeError("Need at least two configurations for jackknife.")

# ------------------------------------------------------------
# Jackknife helper
# ------------------------------------------------------------
def jackknife_stats(samples):
    """
    samples: array-like jackknife estimates, length N
    returns mean, jackknife error
    """
    samples = pd.Series(samples, dtype=float)
    mean_jk = samples.mean()
    err_jk = math.sqrt((N - 1) / N * ((samples - mean_jk) ** 2).sum())
    return mean_jk, err_jk


# ------------------------------------------------------------
# <W(tau)> jackknife
# ------------------------------------------------------------
rows_W = []

for tau in taus:
    full_mean = W[tau].mean()

    jk_means = []
    for nc in W.index:
        jk_sample = W.drop(index=nc)
        jk_means.append(jk_sample[tau].mean())

    mean_jk, err_jk = jackknife_stats(jk_means)

    naive_err = W[tau].std(ddof=1) / math.sqrt(N)

    rows_W.append({
        "tau": tau,
        "Ncfg": N,
        "mean_ReW": full_mean,
        "err_ReW_naive": naive_err,
        "err_ReW_jackknife": err_jk,
    })

df_W = pd.DataFrame(rows_W)


# ------------------------------------------------------------
# Veff(tau) = log(W(tau)/W(tau+1)) jackknife
# ------------------------------------------------------------
rows_V = []

for tau, tau_next in zip(taus[:-1], taus[1:]):
    full_W_tau = W[tau].mean()
    full_W_next = W[tau_next].mean()
    full_veff = math.log(full_W_tau / full_W_next) / (tau_next - tau)

    jk_veffs = []
    for nc in W.index:
        jk_sample = W.drop(index=nc)
        m_tau = jk_sample[tau].mean()
        m_next = jk_sample[tau_next].mean()
        jk_veffs.append(math.log(m_tau / m_next) / (tau_next - tau))

    mean_jk, err_jk = jackknife_stats(jk_veffs)

    rows_V.append({
        "tau": tau,
        "tau_next": tau_next,
        "Ncfg": N,
        "Veff": full_veff,
        "err_Veff_jackknife": err_jk,
    })

df_V = pd.DataFrame(rows_V)


# ------------------------------------------------------------
# Save CSVs
# ------------------------------------------------------------
df_W.to_csv(OUT_DIR / "wilson_tau_ReW_jackknife.csv", index=False)
df_V.to_csv(OUT_DIR / "wilson_tau_Veff_jackknife.csv", index=False)

print("Wrote:")
print(OUT_DIR / "wilson_tau_ReW_jackknife.csv")
print(OUT_DIR / "wilson_tau_Veff_jackknife.csv")

print()
print("Veff jackknife summary:")
print(df_V.to_string(index=False))


# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

plt.figure(figsize=(7, 5))
plt.errorbar(
    df_W["tau"],
    df_W["mean_ReW"],
    yerr=df_W["err_ReW_jackknife"],
    fmt="o-",
    capsize=4,
)
plt.xlabel(r"$\tau/a$")
plt.ylabel(r"$\langle \mathrm{Re}\,W(a,\tau)\rangle$")
plt.title(r"Wilson loop with jackknife errors")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_ReW_vs_tau_jackknife.png", dpi=200)


plt.figure(figsize=(7, 5))
logW = df_W["mean_ReW"].apply(math.log)
logW_err = df_W["err_ReW_jackknife"] / df_W["mean_ReW"]

plt.errorbar(
    df_W["tau"],
    logW,
    yerr=logW_err,
    fmt="o-",
    capsize=4,
)
plt.xlabel(r"$\tau/a$")
plt.ylabel(r"$\log\langle \mathrm{Re}\,W(a,\tau)\rangle$")
plt.title(r"Log Wilson loop with jackknife errors")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_log_ReW_vs_tau_jackknife.png", dpi=200)


plt.figure(figsize=(7, 5))
plt.errorbar(
    df_V["tau"],
    df_V["Veff"],
    yerr=df_V["err_Veff_jackknife"],
    fmt="o-",
    capsize=4,
)
plt.xlabel(r"$\tau/a$")
plt.ylabel(r"$V_{\mathrm{eff}}(a,\tau)$")
plt.title(r"Effective potential with jackknife errors")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_Veff_jackknife.png", dpi=200)

print()
print("Wrote plots:")
print(OUT_DIR / "wilson_tau_ReW_vs_tau_jackknife.png")
print(OUT_DIR / "wilson_tau_log_ReW_vs_tau_jackknife.png")
print(OUT_DIR / "wilson_tau_Veff_jackknife.png")