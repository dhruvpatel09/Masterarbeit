#!/usr/bin/env python3

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

OUT_DIR = Path("analysis_wilson_tau")
SUMMARY_FILE = OUT_DIR / "wilson_tau_summary.csv"

tau = []
mean_ReW = []
err_ReW = []
mean_ImW = []
err_ImW = []
Veff = []
Veff_tau = []

with open(SUMMARY_FILE, "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        t = int(row["tau"])
        tau.append(t)
        mean_ReW.append(float(row["mean_ReW"]))
        err_ReW.append(float(row["err_ReW"]))
        mean_ImW.append(float(row["mean_ImW"]))
        err_ImW.append(float(row["err_ImW"]))

        if row["Veff_to_next"].strip() != "":
            Veff_tau.append(t)
            Veff.append(float(row["Veff_to_next"]))

# --------------------------------------------------
# Plot 1: Re W vs tau
# --------------------------------------------------
plt.figure(figsize=(7,5))
plt.errorbar(tau, mean_ReW, yerr=err_ReW, fmt='o-', capsize=4)
plt.xlabel(r'$\tau/a$')
plt.ylabel(r'$\langle \mathrm{Re}\,W(a,\tau)\rangle$')
plt.title(r'Wilson loop at fixed spatial extent $a$')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_ReW_vs_tau.png", dpi=200)

# --------------------------------------------------
# Plot 2: log(Re W) vs tau
# --------------------------------------------------
log_ReW = [math.log(x) for x in mean_ReW]
log_err = [e/x for x, e in zip(mean_ReW, err_ReW)]

plt.figure(figsize=(7,5))
plt.errorbar(tau, log_ReW, yerr=log_err, fmt='o-', capsize=4)
plt.xlabel(r'$\tau/a$')
plt.ylabel(r'$\log\langle \mathrm{Re}\,W(a,\tau)\rangle$')
plt.title(r'Log-decay of Wilson loop')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_log_ReW_vs_tau.png", dpi=200)

# --------------------------------------------------
# Plot 3: Im W vs tau
# --------------------------------------------------
plt.figure(figsize=(7,5))
plt.errorbar(tau, mean_ImW, yerr=err_ImW, fmt='o-', capsize=4)
plt.xlabel(r'$\tau/a$')
plt.ylabel(r'$\langle \mathrm{Im}\,W(a,\tau)\rangle$')
plt.title(r'Imaginary part of Wilson loop')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_ImW_vs_tau.png", dpi=200)

# --------------------------------------------------
# Plot 4: Veff vs tau
# --------------------------------------------------
plt.figure(figsize=(7,5))
plt.plot(Veff_tau, Veff, 'o-')
plt.xlabel(r'$\tau/a$')
plt.ylabel(r'$V_{\mathrm{eff}}(a,\tau)$')
plt.title(r'Effective potential estimate')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / "wilson_tau_Veff.png", dpi=200)

print("Wrote plots:")
print(OUT_DIR / "wilson_tau_ReW_vs_tau.png")
print(OUT_DIR / "wilson_tau_log_ReW_vs_tau.png")
print(OUT_DIR / "wilson_tau_ImW_vs_tau.png")
print(OUT_DIR / "wilson_tau_Veff.png")
