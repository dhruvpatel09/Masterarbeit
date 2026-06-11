#!/usr/bin/env python3

import csv
import math
from collections import defaultdict

inp = "results/spatial_avg_RT_R0-3_T1-4_n1-100.csv"

rows = []
with open(inp, newline="") as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append({
            "cfg": int(r["cfg"]),
            "T": int(r["T"]),
            "R": int(r["R"]),
            "Nsrc": int(r["Nsrc"]),
            "Re": float(r["Re"]),
            "Im": float(r["Im"]),
        })

cfgs = sorted({r["cfg"] for r in rows})
N = len(cfgs)

by_TR = defaultdict(list)
by_cfg_TR = {}

for r in rows:
    by_TR[(r["T"], r["R"])].append(r)
    by_cfg_TR[(r["cfg"], r["T"], r["R"])] = r

def mean(xs):
    return sum(xs) / len(xs)

def jackknife_error(samples):
    n = len(samples)
    m = mean(samples)
    return math.sqrt((n - 1) / n * sum((x - m)**2 for x in samples))

# Jackknife for mean L(R,T)
out_L = "results/spatial_avg_RT_R0-3_T1-4_n1-100_jackknife.csv"
with open(out_L, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "T", "R", "Ncfg",
        "mean_Re", "err_Re_jackknife",
        "mean_Im", "err_Im_jackknife",
    ])

    for T in range(1, 5):
        for R in range(0, 4):
            re_all = [by_cfg_TR[(cfg, T, R)]["Re"] for cfg in cfgs]
            im_all = [by_cfg_TR[(cfg, T, R)]["Im"] for cfg in cfgs]

            mean_re_all = mean(re_all)
            mean_im_all = mean(im_all)

            jk_re = []
            jk_im = []

            for omit in cfgs:
                re_loo = [by_cfg_TR[(cfg, T, R)]["Re"] for cfg in cfgs if cfg != omit]
                im_loo = [by_cfg_TR[(cfg, T, R)]["Im"] for cfg in cfgs if cfg != omit]
                jk_re.append(mean(re_loo))
                jk_im.append(mean(im_loo))

            w.writerow([
                T, R, N,
                f"{mean_re_all:.16e}", f"{jackknife_error(jk_re):.16e}",
                f"{mean_im_all:.16e}", f"{jackknife_error(jk_im):.16e}",
            ])

# Jackknife for Veff from ratio of leave-one-out means
out_V = "results/laplace_Veff_R0-3_T1-3_n1-100_jackknife.csv"
with open(out_V, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "T", "R", "Ncfg",
        "Veff_log_ratio_of_means",
        "err_Veff_jackknife",
    ])

    for T in range(1, 4):
        for R in range(0, 4):
            re_T_all = [by_cfg_TR[(cfg, T, R)]["Re"] for cfg in cfgs]
            re_Tp1_all = [by_cfg_TR[(cfg, T + 1, R)]["Re"] for cfg in cfgs]

            mean_T = mean(re_T_all)
            mean_Tp1 = mean(re_Tp1_all)

            if mean_T <= 0 or mean_Tp1 <= 0:
                veff_all = float("nan")
            else:
                veff_all = math.log(mean_T / mean_Tp1)

            jk_veff = []

            for omit in cfgs:
                re_T_loo = [by_cfg_TR[(cfg, T, R)]["Re"] for cfg in cfgs if cfg != omit]
                re_Tp1_loo = [by_cfg_TR[(cfg, T + 1, R)]["Re"] for cfg in cfgs if cfg != omit]

                mT = mean(re_T_loo)
                mTp1 = mean(re_Tp1_loo)

                if mT <= 0 or mTp1 <= 0:
                    continue

                jk_veff.append(math.log(mT / mTp1))

            w.writerow([
                T, R, N,
                f"{veff_all:.16e}",
                f"{jackknife_error(jk_veff):.16e}",
            ])

print("Rows read:", len(rows))
print("Configs:", N)
print("Wrote:")
print(" ", out_L)
print(" ", out_V)
