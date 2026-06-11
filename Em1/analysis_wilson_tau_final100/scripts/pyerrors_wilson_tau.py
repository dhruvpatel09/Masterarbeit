#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import pyerrors as pe
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
PER_CONFIG = SCRIPT_DIR / "wilson_tau_per_config.csv"
OUT_W = SCRIPT_DIR / "wilson_tau_ReW_pyerrors.csv"
OUT_V = SCRIPT_DIR / "wilson_tau_Veff_pyerrors.csv"

ENSEMBLE = "Em1p4"
S_VALUE = 2.0  # pyerrors default is commonly S=2; keep explicit.

df = pd.read_csv(PER_CONFIG)

Wmat = df.pivot(index="nc", columns="tau", values="ReW").sort_index()
cfg_ids = list(Wmat.index)
taus = list(Wmat.columns)

print(f"Loaded {len(cfg_ids)} configurations and taus: {taus}")

obs = {}

for tau in taus:
    samples = Wmat[tau].to_numpy(dtype=float)

    # idl tells pyerrors which configuration numbers these samples belong to.
    # For our current 1..30 data this is regular, but this also keeps things explicit.
    o = pe.Obs([samples], [ENSEMBLE], idl=[cfg_ids])
    o.gamma_method(S=S_VALUE)

    obs[tau] = o

rows_W = []

for tau in taus:
    o = obs[tau]

    rows_W.append({
        "tau": tau,
        "Ncfg": len(cfg_ids),
        "mean_ReW": o.value,
        "err_pyerrors": o.dvalue,
        "dd_err_pyerrors": o.ddvalue,
    })

rows_V = []

for tau, tau_next in zip(taus[:-1], taus[1:]):
    delta_tau = tau_next - tau

    v = np.log(obs[tau] / obs[tau_next]) / delta_tau
    v.gamma_method(S=S_VALUE)

    rows_V.append({
        "tau": tau,
        "tau_next": tau_next,
        "Ncfg": len(cfg_ids),
        "Veff": v.value,
        "err_pyerrors": v.dvalue,
        "dd_err_pyerrors": v.ddvalue,
    })

df_W = pd.DataFrame(rows_W)
df_V = pd.DataFrame(rows_V)

df_W.to_csv(OUT_W, index=False)
df_V.to_csv(OUT_V, index=False)

print("\n<W> from pyerrors:")
print(df_W.to_string(index=False))

print("\nVeff from pyerrors:")
print(df_V.to_string(index=False))

# Compare Veff against jackknife if available
jk_file = SCRIPT_DIR / "wilson_tau_Veff_jackknife.csv"
if jk_file.exists():
    jk = pd.read_csv(jk_file)
    cmp = df_V.merge(
        jk[["tau", "tau_next", "err_Veff_jackknife"]],
        on=["tau", "tau_next"],
        how="left",
    )
    cmp["ratio_pyerrors_over_jackknife"] = cmp["err_pyerrors"] / cmp["err_Veff_jackknife"]
    cmp.to_csv(SCRIPT_DIR / "wilson_tau_Veff_pyerrors_vs_jackknife.csv", index=False)

    print("\nVeff error comparison:")
    print(cmp.to_string(index=False))

    plt.figure(figsize=(7, 5))
    plt.errorbar(
        cmp["tau"] - 0.04,
        cmp["Veff"],
        yerr=cmp["err_Veff_jackknife"],
        fmt="o",
        capsize=4,
        label="jackknife",
    )
    plt.errorbar(
        cmp["tau"] + 0.04,
        cmp["Veff"],
        yerr=cmp["err_pyerrors"],
        fmt="s",
        capsize=4,
        label="pyerrors / Gamma method",
    )
    plt.xlabel(r"$\tau/a$")
    plt.ylabel(r"$V_{\mathrm{eff}}(a,\tau)$")
    plt.title(r"$V_{\mathrm{eff}}$: jackknife vs pyerrors")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / "wilson_tau_Veff_pyerrors_vs_jackknife.png", dpi=200)

print("\nWrote:")
print(OUT_W)
print(OUT_V)
if jk_file.exists():
    print(SCRIPT_DIR / "wilson_tau_Veff_pyerrors_vs_jackknife.csv")
    print(SCRIPT_DIR / "wilson_tau_Veff_pyerrors_vs_jackknife.png")