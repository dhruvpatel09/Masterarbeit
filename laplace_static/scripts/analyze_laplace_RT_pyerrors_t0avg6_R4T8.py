#!/usr/bin/env python3

from pathlib import Path

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
    "spatial_avg_RT_R0-4_T1-8_t0avg6_n1-100.csv"
)

JK_FILE = (
    RESULTS /
    "laplace_Veff_R0-4_T1-7_t0avg6_n1-100_jackknife.csv"
)

OUT_L = (
    RESULTS /
    "spatial_avg_RT_R0-4_T1-8_t0avg6_n1-100_pyerrors.csv"
)

OUT_V = (
    RESULTS /
    "laplace_Veff_R0-4_T1-7_t0avg6_n1-100_pyerrors.csv"
)

OUT_CMP = (
    RESULTS /
    "laplace_Veff_R0-4_T1-7_t0avg6_n1-100_error_comparison.csv"
)

ENSEMBLE = "Em1p4_laplace_t0avg6_R4T8"
S_VALUE = 2.0

EXPECTED_CFGS = list(range(1, 101))
EXPECTED_T = list(range(1, 9))
EXPECTED_R = list(range(0, 5))
EXPECTED_NSRC = 6 * 24 ** 3


def validate_input(df):
    required = {"cfg", "T", "R", "Nsrc", "Re", "Im"}

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"Missing columns: {sorted(missing)}"
        )

    if df.duplicated(["cfg", "T", "R"]).any():
        raise RuntimeError(
            "Duplicate (cfg,T,R) rows found"
        )

    cfgs = sorted(df["cfg"].unique().tolist())
    t_values = sorted(df["T"].unique().tolist())
    r_values = sorted(df["R"].unique().tolist())

    if cfgs != EXPECTED_CFGS:
        raise RuntimeError(
            f"Unexpected configuration IDs: {cfgs}"
        )

    if t_values != EXPECTED_T:
        raise RuntimeError(
            f"Unexpected T values: {t_values}"
        )

    if r_values != EXPECTED_R:
        raise RuntimeError(
            f"Unexpected R values: {r_values}"
        )

    if len(df) != 4000:
        raise RuntimeError(
            f"Expected 4000 rows, found {len(df)}"
        )

    if not (df["Nsrc"] == EXPECTED_NSRC).all():
        raise RuntimeError(
            "Unexpected Nsrc values found"
        )

    return cfgs


def finite_group(frame, columns):
    cleaned = frame.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    return cleaned.dropna(subset=columns)


def main():
    df = pd.read_csv(INPUT)
    df = df.sort_values(
        ["cfg", "T", "R"]
    ).reset_index(drop=True)

    cfg_ids = validate_input(df)
    ncfg = len(cfg_ids)

    obs = {}
    rows_L = []

    for T in EXPECTED_T:
        for R in EXPECTED_R:
            group = (
                df[
                    (df["T"] == T) &
                    (df["R"] == R)
                ]
                .sort_values("cfg")
            )

            if group["cfg"].tolist() != cfg_ids:
                raise RuntimeError(
                    f"Configuration mismatch at T={T}, R={R}"
                )

            samples = group[
                "Re"
            ].to_numpy(dtype=float)

            observable = pe.Obs(
                [samples],
                [ENSEMBLE],
                idl=[cfg_ids],
            )

            observable.gamma_method(S=S_VALUE)
            obs[(T, R)] = observable

            relative_error = (
                100.0
                * observable.dvalue
                / abs(observable.value)
                if observable.value != 0.0
                else np.nan
            )

            snr = (
                abs(observable.value)
                / observable.dvalue
                if observable.dvalue > 0.0
                else np.inf
            )

            rows_L.append({
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "mean_Re": observable.value,
                "err_pyerrors": observable.dvalue,
                "dd_err_pyerrors": observable.ddvalue,
                "relative_err_Re_percent": relative_error,
                "snr_Re": snr,
                "n_nonpositive_Re": int(
                    np.count_nonzero(samples <= 0.0)
                ),
                "S": S_VALUE,
            })

    rows_V = []

    for T in EXPECTED_T[:-1]:
        for R in EXPECTED_R:
            numerator = obs[(T, R)]
            denominator = obs[(T + 1, R)]

            if (
                numerator.value <= 0.0 or
                denominator.value <= 0.0
            ):
                rows_V.append({
                    "T": T,
                    "R": R,
                    "Ncfg": ncfg,
                    "Veff": np.nan,
                    "err_pyerrors": np.nan,
                    "dd_err_pyerrors": np.nan,
                    "relative_err_percent": np.nan,
                    "S": S_VALUE,
                    "status": "nonpositive_central_mean",
                })
                continue

            try:
                veff = np.log(
                    numerator / denominator
                )
                veff.gamma_method(S=S_VALUE)

                relative_error = (
                    100.0
                    * veff.dvalue
                    / abs(veff.value)
                    if veff.value != 0.0
                    else np.nan
                )

                rows_V.append({
                    "T": T,
                    "R": R,
                    "Ncfg": ncfg,
                    "Veff": veff.value,
                    "err_pyerrors": veff.dvalue,
                    "dd_err_pyerrors": veff.ddvalue,
                    "relative_err_percent": relative_error,
                    "S": S_VALUE,
                    "status": "ok",
                })

            except Exception as error:
                rows_V.append({
                    "T": T,
                    "R": R,
                    "Ncfg": ncfg,
                    "Veff": np.nan,
                    "err_pyerrors": np.nan,
                    "dd_err_pyerrors": np.nan,
                    "relative_err_percent": np.nan,
                    "S": S_VALUE,
                    "status": (
                        "pyerrors_failure_"
                        + type(error).__name__
                    ),
                })

                print(
                    "WARNING: pyerrors Veff failed "
                    f"for T={T}, R={R}: {error}"
                )

    out_L = pd.DataFrame(rows_L)
    out_V = pd.DataFrame(rows_V)

    out_L.to_csv(
        OUT_L,
        index=False,
        lineterminator="\n",
    )

    out_V.to_csv(
        OUT_V,
        index=False,
        lineterminator="\n",
    )

    jk = pd.read_csv(JK_FILE)

    comparison = out_V.merge(
        jk[[
            "T",
            "R",
            "Veff_log_ratio_of_means",
            "err_Veff_jackknife",
            "relative_err_Veff_percent",
            "n_valid_jackknife",
            "status",
        ]].rename(columns={
            "Veff_log_ratio_of_means":
                "Veff_jackknife_central",
            "relative_err_Veff_percent":
                "relative_err_jackknife_percent",
            "status":
                "status_jackknife",
        }),
        on=["T", "R"],
        how="left",
        validate="one_to_one",
    )

    valid_ratio = (
        np.isfinite(comparison["err_pyerrors"]) &
        np.isfinite(
            comparison["err_Veff_jackknife"]
        ) &
        (
            comparison["err_Veff_jackknife"] > 0.0
        )
    )

    comparison[
        "ratio_pyerrors_over_jackknife"
    ] = np.where(
        valid_ratio,
        comparison["err_pyerrors"]
        / comparison["err_Veff_jackknife"],
        np.nan,
    )

    comparison.to_csv(
        OUT_CMP,
        index=False,
        lineterminator="\n",
    )

    # Correlator versus T
    plt.figure(figsize=(7.4, 5.0))

    for R in EXPECTED_R:
        group = finite_group(
            out_L[
                (out_L["R"] == R) &
                (out_L["mean_Re"] > 0.0)
            ].sort_values("T"),
            ["mean_Re", "err_pyerrors"],
        )

        plt.errorbar(
            group["T"],
            group["mean_Re"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=3,
            label=f"R/a = {R}",
        )

    plt.yscale("log")
    plt.xlabel("T/a")
    plt.ylabel("Re L(R,T)")
    plt.title("Laplace correlator: six-source temporal average")
    plt.xticks(EXPECTED_T)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_ReL_vs_T_pyerrors_t0avg6_R4T8.png",
        dpi=200,
    )
    plt.close()

    # Effective potential versus T
    plt.figure(figsize=(7.4, 5.0))

    for R in EXPECTED_R[1:]:
        group = finite_group(
            comparison[
                comparison["R"] == R
            ].sort_values("T"),
            ["Veff", "err_pyerrors"],
        )

        plt.errorbar(
            group["T"],
            group["Veff"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=4,
            elinewidth=1.3,
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Laplace effective potential: six-source temporal average")
    plt.xticks(EXPECTED_T[:-1])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_vs_T_pyerrors_t0avg6_R4T8.png",
        dpi=200,
    )
    plt.close()

    # Effective potential versus R
    plt.figure(figsize=(7.4, 5.0))

    for T in EXPECTED_T[:-1]:
        group = finite_group(
            comparison[
                (comparison["T"] == T) &
                (comparison["R"] > 0)
            ].sort_values("R"),
            ["Veff", "err_pyerrors"],
        )

        plt.errorbar(
            group["R"],
            group["Veff"],
            yerr=group["err_pyerrors"],
            marker="o",
            capsize=3,
            label=f"T/a = {T}",
        )

    plt.xlabel("R/a")
    plt.ylabel("a Veff(R,T)")
    plt.title("Laplace effective potential versus separation (six-source average)")
    plt.xticks(EXPECTED_R[1:])
    plt.grid(True, alpha=0.3)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(
        FIGDIR / "laplace_Veff_vs_R_pyerrors_t0avg6_R4T8.png",
        dpi=200,
    )
    plt.close()

    # Relative uncertainty
    plt.figure(figsize=(7.4, 5.0))

    for R in EXPECTED_R[1:]:
        group = finite_group(
            comparison[
                comparison["R"] == R
            ].sort_values("T"),
            ["relative_err_percent"],
        )

        plt.plot(
            group["T"],
            group["relative_err_percent"],
            marker="o",
            label=f"R/a = {R}",
        )

    plt.xlabel("T/a")
    plt.ylabel("Relative uncertainty [%]")
    plt.title("Relative uncertainty of Veff (six-source temporal average)")
    plt.xticks(EXPECTED_T[:-1])
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        FIGDIR /
        "laplace_Veff_relative_error_pyerrors_t0avg6_R4T8.png",
        dpi=200,
    )
    plt.close()

    # Per-R zoom and method comparison
    for R in EXPECTED_R[1:]:
        group = comparison[
            comparison["R"] == R
        ].sort_values("T")

        py_group = finite_group(
            group,
            ["Veff", "err_pyerrors"],
        )

        plt.figure(figsize=(7.0, 4.8))
        plt.errorbar(
            py_group["T"],
            py_group["Veff"],
            yerr=py_group["err_pyerrors"],
            marker="o",
            capsize=5,
            elinewidth=1.4,
        )
        plt.xlabel("T/a")
        plt.ylabel("a Veff(R,T)")
        plt.title(
            f"Laplace effective potential, R/a = {R}"
        )
        plt.xticks(EXPECTED_T[:-1])
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(
            FIGDIR /
            f"laplace_Veff_pyerrors_zoom_R{R}_t0avg6_R4T8.png",
            dpi=200,
        )
        plt.close()

        jk_group = finite_group(
            group,
            [
                "Veff",
                "err_Veff_jackknife",
            ],
        )

        plt.figure(figsize=(7.0, 4.8))

        plt.errorbar(
            jk_group["T"] - 0.04,
            jk_group["Veff"],
            yerr=jk_group[
                "err_Veff_jackknife"
            ],
            fmt="o",
            capsize=4,
            label="Jackknife",
        )

        plt.errorbar(
            py_group["T"] + 0.04,
            py_group["Veff"],
            yerr=py_group["err_pyerrors"],
            fmt="^",
            capsize=4,
            label="Gamma method (pyerrors)",
        )

        plt.xlabel("T/a")
        plt.ylabel("a Veff(R,T)")
        plt.title(
            f"Error-method comparison, R/a = {R}"
        )
        plt.xticks(EXPECTED_T[:-1])
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            FIGDIR /
            f"laplace_Veff_methods_R{R}_t0avg6_R4T8.png",
            dpi=200,
        )
        plt.close()

    print("pyerrors Veff status counts:")
    print(
        out_V["status"]
        .value_counts(dropna=False)
        .to_string()
    )

    print()
    print("Physical Veff results:")
    print(
        comparison[
            comparison["R"] > 0
        ][[
            "T",
            "R",
            "Veff",
            "err_pyerrors",
            "relative_err_percent",
            "err_Veff_jackknife",
            "ratio_pyerrors_over_jackknife",
            "status",
            "status_jackknife",
        ]].to_string(index=False)
    )

    print()
    print("Wrote:")
    print(OUT_L)
    print(OUT_V)
    print(OUT_CMP)

    for path in sorted(
        FIGDIR.glob("*t0avg6_R4T8.png")
    ):
        print(path)


if __name__ == "__main__":
    main()
