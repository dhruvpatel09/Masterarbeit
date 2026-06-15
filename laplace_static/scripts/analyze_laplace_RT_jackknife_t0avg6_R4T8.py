#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

INPUT = RESULTS / "spatial_avg_RT_R0-4_T1-8_t0avg6_n1-100.csv"

OUT_L = (
    RESULTS /
    "spatial_avg_RT_R0-4_T1-8_t0avg6_n1-100_jackknife.csv"
)

OUT_V = (
    RESULTS /
    "laplace_Veff_R0-4_T1-7_t0avg6_n1-100_jackknife.csv"
)

EXPECTED_CFGS = list(range(1, 101))
EXPECTED_T = list(range(1, 9))
EXPECTED_R = list(range(0, 5))
EXPECTED_NSRC = 6 * 24 ** 3


def validate_input(df):
    required = {"cfg", "T", "R", "Nsrc", "Re", "Im"}

    missing_columns = required - set(df.columns)

    if missing_columns:
        raise RuntimeError(
            f"Missing columns: {sorted(missing_columns)}"
        )

    if df.duplicated(["cfg", "T", "R"]).any():
        duplicates = df[
            df.duplicated(
                ["cfg", "T", "R"],
                keep=False,
            )
        ]
        raise RuntimeError(
            "Duplicate (cfg,T,R) rows:\n"
            + duplicates.to_string(index=False)
        )

    cfgs = sorted(df["cfg"].unique().tolist())
    t_values = sorted(df["T"].unique().tolist())
    r_values = sorted(df["R"].unique().tolist())

    if cfgs != EXPECTED_CFGS:
        raise RuntimeError(
            f"Configuration IDs differ: {cfgs}"
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

    bad_nsrc = df[df["Nsrc"] != EXPECTED_NSRC]

    if not bad_nsrc.empty:
        raise RuntimeError(
            "Unexpected Nsrc values:\n"
            + bad_nsrc.to_string(index=False)
        )

    expected_grid = {
        (T, R)
        for T in EXPECTED_T
        for R in EXPECTED_R
    }

    for cfg, group in df.groupby("cfg"):
        grid = set(zip(group["T"], group["R"]))

        if grid != expected_grid:
            raise RuntimeError(
                f"Grid mismatch for cfg={cfg}: "
                f"missing={sorted(expected_grid - grid)}, "
                f"extra={sorted(grid - expected_grid)}"
            )

    return cfgs


def jackknife_error(samples):
    samples = np.asarray(samples, dtype=float)
    n = len(samples)

    if n < 2:
        return np.nan

    center = samples.mean()

    return np.sqrt(
        (n - 1) / n
        * np.sum((samples - center) ** 2)
    )


def main():
    df = pd.read_csv(INPUT)
    df = df.sort_values(
        ["cfg", "T", "R"]
    ).reset_index(drop=True)

    cfgs = validate_input(df)
    ncfg = len(cfgs)

    arrays_re = {}
    arrays_im = {}

    for T in EXPECTED_T:
        for R in EXPECTED_R:
            group = (
                df[
                    (df["T"] == T) &
                    (df["R"] == R)
                ]
                .sort_values("cfg")
            )

            if group["cfg"].tolist() != cfgs:
                raise RuntimeError(
                    f"Configuration mismatch at T={T}, R={R}"
                )

            arrays_re[(T, R)] = (
                group["Re"].to_numpy(dtype=float)
            )
            arrays_im[(T, R)] = (
                group["Im"].to_numpy(dtype=float)
            )

    rows_L = []

    for T in EXPECTED_T:
        for R in EXPECTED_R:
            re_values = arrays_re[(T, R)]
            im_values = arrays_im[(T, R)]

            sum_re = re_values.sum()
            sum_im = im_values.sum()

            jk_re = (
                sum_re - re_values
            ) / (ncfg - 1)

            jk_im = (
                sum_im - im_values
            ) / (ncfg - 1)

            mean_re = re_values.mean()
            mean_im = im_values.mean()

            err_re = jackknife_error(jk_re)
            err_im = jackknife_error(jk_im)

            relative_error = (
                100.0 * err_re / abs(mean_re)
                if mean_re != 0.0
                else np.nan
            )

            snr = (
                abs(mean_re) / err_re
                if err_re > 0.0
                else np.inf
            )

            rows_L.append({
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "mean_Re": mean_re,
                "err_Re_jackknife": err_re,
                "relative_err_Re_percent": relative_error,
                "snr_Re": snr,
                "n_nonpositive_Re": int(
                    np.count_nonzero(re_values <= 0.0)
                ),
                "mean_Im": mean_im,
                "err_Im_jackknife": err_im,
            })

    rows_V = []

    for T in EXPECTED_T[:-1]:
        for R in EXPECTED_R:
            values_T = arrays_re[(T, R)]
            values_Tp1 = arrays_re[(T + 1, R)]

            mean_T = values_T.mean()
            mean_Tp1 = values_Tp1.mean()

            jk_mean_T = (
                values_T.sum() - values_T
            ) / (ncfg - 1)

            jk_mean_Tp1 = (
                values_Tp1.sum() - values_Tp1
            ) / (ncfg - 1)

            valid_jk = (
                (jk_mean_T > 0.0) &
                (jk_mean_Tp1 > 0.0)
            )

            n_valid = int(valid_jk.sum())

            if mean_T <= 0.0 or mean_Tp1 <= 0.0:
                veff = np.nan
                err_veff = np.nan
                status = "nonpositive_central_mean"

            elif n_valid != ncfg:
                veff = np.log(mean_T / mean_Tp1)
                err_veff = np.nan
                status = "nonpositive_leave_one_out_mean"

            else:
                veff = np.log(mean_T / mean_Tp1)

                jk_veff = np.log(
                    jk_mean_T / jk_mean_Tp1
                )

                err_veff = jackknife_error(jk_veff)
                status = "ok"

            relative_error = (
                100.0 * err_veff / abs(veff)
                if (
                    np.isfinite(veff) and
                    np.isfinite(err_veff) and
                    veff != 0.0
                )
                else np.nan
            )

            rows_V.append({
                "T": T,
                "R": R,
                "Ncfg": ncfg,
                "mean_Re_T": mean_T,
                "mean_Re_Tp1": mean_Tp1,
                "Veff_log_ratio_of_means": veff,
                "err_Veff_jackknife": err_veff,
                "relative_err_Veff_percent": relative_error,
                "n_valid_jackknife": n_valid,
                "status": status,
            })

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

    print(f"Input rows: {len(df)}")
    print(f"Configurations: {ncfg}")
    print(f"L(R,T) rows: {len(out_L)}")
    print(f"Veff rows: {len(out_V)}")

    print()
    print("Veff status counts:")
    print(
        out_V["status"]
        .value_counts(dropna=False)
        .to_string()
    )

    invalid = out_V[out_V["status"] != "ok"]

    if not invalid.empty:
        print()
        print("Non-standard Veff rows:")
        print(
            invalid[[
                "T",
                "R",
                "mean_Re_T",
                "mean_Re_Tp1",
                "n_valid_jackknife",
                "status",
            ]].to_string(index=False)
        )

    print()
    print("Wrote:")
    print(OUT_L)
    print(OUT_V)


if __name__ == "__main__":
    main()
