#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyerrors as pe


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

INPUT = ROOT / "results/spatial_avg_RT_R0-3_T1-4_n1-100.csv"
JK_FILE = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_jackknife.csv"
GAMMA_FILE = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_gamma.csv"

OUT_L = ROOT / "results/spatial_avg_RT_R0-3_T1-4_n1-100_pyerrors.csv"
OUT_V = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_pyerrors.csv"
OUT_CMP = ROOT / "results/laplace_Veff_R0-3_T1-3_n1-100_error_comparison.csv"

FIGDIR = ROOT / "results/figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

ENSEMBLE = "Em1p4_laplace"
S_VALUE = 2.0


def main():
    df = pd.read_csv(INPUT)
    df = df.sort_values(["cfg", "T", "R"]).reset_index(drop=True)

    cfg_ids = sorted(df["cfg"].unique())
    ncfg = len(cfg_ids)

    obs = {}
    rows_L = []

    for T in range(1, 5):
        for R in range(0, 4):
            g = df[
                (df["T"] == T) &
                (df["R"] == R)
            ].sort_values("cfg")

            if g["cfg"].tolist() != cfg_ids:
                raise RuntimeError(
                    f"Configuration mismatch for T={T}, R={R}"
                )

            samples = g["Re"].to_numpy(dtype=float)

            o = pe.Obs(
                [samples],
                [ENSEMBLE],
                idl=[cfg_ids],
            )
            o.gamma_method(S=S_VALUE)
            obs[(T, R)] = o

            rows_L.append({
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "mean_Re": o.value,
                "err_pyerrors": o.dvalue,
                "dd_err_pyerrors": o.ddvalue,
                "S": S_VALUE,
            })

    rows_V = []

    for T in range(1, 4):
        for R in range(0, 4):
            v = np.log(
                obs[(T, R)] / obs[(T + 1, R)]
            )
            v.gamma_method(S=S_VALUE)

            rows_V.append({
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "Veff": v.value,
                "err_pyerrors": v.dvalue,
                "dd_err_pyerrors": v.ddvalue,
                "S": S_VALUE,
            })

    out_L = pd.DataFrame(rows_L)
    out_V = pd.DataFrame(rows_V)

    out_L.to_csv(OUT_L, index=False)
    out_V.to_csv(OUT_V, index=False)

    jk = pd.read_csv(JK_FILE)
    gamma = pd.read_csv(GAMMA_FILE)

    cmp = (
        out_V
        .merge(
            jk[[
                "T", "R",
                "err_Veff_jackknife",
            ]],
            on=["T", "R"],
            how="left",
        )
        .merge(
            gamma[[
                "T", "R",
                "err_gamma",
                "tauint",
                "Wopt",
            ]],
            on=["T", "R"],
            how="left",
        )
    )

    cmp = cmp.rename(columns={
        "err_gamma": "err_custom_gamma",
        "tauint": "tauint_custom_gamma",
        "Wopt": "Wopt_custom_gamma",
    })

    cmp["ratio_custom_gamma_over_jackknife"] = (
        cmp["err_custom_gamma"]
        / cmp["err_Veff_jackknife"]
    )
    cmp["ratio_pyerrors_over_jackknife"] = (
        cmp["err_pyerrors"]
        / cmp["err_Veff_jackknife"]
    )
    cmp["ratio_pyerrors_over_custom_gamma"] = (
        cmp["err_pyerrors"]
        / cmp["err_custom_gamma"]
    )

    cmp.to_csv(OUT_CMP, index=False)

    print("pyerrors L(R,T):")
    print(out_L.to_string(index=False))

    print()
    print("pyerrors Veff(R,T):")
    print(out_V.to_string(index=False))

    print()
    print("Error comparison:")
    print(cmp.to_string(index=False))

    # Final physics plot: pyerrors only
    plt.figure(figsize=(7.0, 4.8))

    for R in range(1, 4):
        g = cmp[cmp["R"] == R]

        plt.errorbar(
            g["T"],
            g["Veff"],
            yerr=g["err_pyerrors"],
            marker="o",
            capsize=4,
            elinewidth=1.3,
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Laplace effective potential: pyerrors")
    plt.xticks([1, 2, 3])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_vs_T_pyerrors.png",
        dpi=200,
    )
    plt.close()

    # Potential versus R with pyerrors
    plt.figure(figsize=(7.0, 4.8))

    for T in range(1, 4):
        g = cmp[
            (cmp["T"] == T) &
            (cmp["R"] > 0)
        ]

        plt.errorbar(
            g["R"],
            g["Veff"],
            yerr=g["err_pyerrors"],
            marker="o",
            capsize=4,
            elinewidth=1.3,
            label=f"T/a = {T}",
        )

    plt.xlabel("R/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Laplace effective potential versus separation")
    plt.xticks([1, 2, 3])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_vs_R_pyerrors.png",
        dpi=200,
    )
    plt.close()

    # Offset comparison: one uncluttered figure per R
    offsets = {
        "jackknife": -0.055,
        "custom_gamma": 0.0,
        "pyerrors": 0.055,
    }

    for R in range(1, 4):
        g = cmp[cmp["R"] == R]

        plt.figure(figsize=(7.0, 4.8))

        plt.errorbar(
            g["T"] + offsets["jackknife"],
            g["Veff"],
            yerr=g["err_Veff_jackknife"],
            fmt="o",
            capsize=4,
            elinewidth=1.3,
            label="Jackknife",
        )

        plt.errorbar(
            g["T"] + offsets["custom_gamma"],
            g["Veff"],
            yerr=g["err_custom_gamma"],
            fmt="s",
            capsize=4,
            elinewidth=1.3,
            label="Custom Gamma",
        )

        plt.errorbar(
            g["T"] + offsets["pyerrors"],
            g["Veff"],
            yerr=g["err_pyerrors"],
            fmt="^",
            capsize=4,
            elinewidth=1.3,
            label="pyerrors",
        )

        plt.xlabel("T/a")
        plt.ylabel("a Veff(R,T)")
        plt.title(f"Error-method comparison, R/a = {R}")
        plt.xticks([1, 2, 3])
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            FIGDIR / f"laplace_Veff_methods_R{R}.png",
            dpi=200,
        )
        plt.close()

    print()
    print("Wrote:")
    print(OUT_L)
    print(OUT_V)
    print(OUT_CMP)
    print(FIGDIR / "laplace_Veff_vs_T_pyerrors.png")
    print(FIGDIR / "laplace_Veff_vs_R_pyerrors.png")
    for R in range(1, 4):
        print(FIGDIR / f"laplace_Veff_methods_R{R}.png")


if __name__ == "__main__":
    main()
