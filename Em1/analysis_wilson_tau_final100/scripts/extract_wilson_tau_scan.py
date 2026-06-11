#!/usr/bin/env python3

import csv
import math
import re
import subprocess
import os
from pathlib import Path
from statistics import mean, stdev

DAT_DIR = Path("dat")
OUT_DIR = Path("analysis_wilson_tau")
OUT_DIR.mkdir(exist_ok=True)

pattern = "Em1p4_wilson_mu0_tau*.wilson_tau.bdio"


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

def strip_ansi(text):
    return ANSI_RE.sub("", text)

def run(cmd):
    env = os.environ.copy()
    env["TERM"] = "dumb"
    env["NO_COLOR"] = "1"
    out = subprocess.check_output(cmd, text=True, env=env)
    return strip_ansi(out)


def parse_tau_from_metadata(filename):
    meta = run(["lsbdio", "-c", "0", "-d", "2", str(filename)])
    m = re.search(r"^TAU=(\d+)$", meta, re.MULTILINE)
    if not m:
        raise RuntimeError(f"Could not find TAU in metadata of {filename}")
    return int(m.group(1))


def find_data_record_ids(filename):
    table = run(["lsbdio", str(filename)])
    rec_ids = []
    for line in table.splitlines():
        parts = line.split()
        if len(parts) >= 5:
            # Example line:
            # 4 record f64 le 32 byte 5 ...
            if parts[1:4] == ["record", "f64", "le"]:
                rec_ids.append(int(parts[0]))
    return rec_ids


def dump_record(filename, rec_id):
    txt = run(["lsbdio", "-c", "0", "-d", str(rec_id), str(filename)])
    vals = [float(x) for x in txt.split()]
    if len(vals) != 4:
        raise RuntimeError(
            f"Expected 4 doubles in record {rec_id} of {filename}, got {len(vals)}"
        )
    nc, corr_id, re_w, im_w = vals
    return int(round(nc)), int(round(corr_id)), re_w, im_w


all_rows = []
summary_rows = []

files = sorted(DAT_DIR.glob(pattern))

if not files:
    raise SystemExit(f"No files found matching dat/{pattern}")

for filename in files:
    tau = parse_tau_from_metadata(filename)
    rec_ids = find_data_record_ids(filename)

    rows = []
    for rec_id in rec_ids:
        nc, corr_id, re_w, im_w = dump_record(filename, rec_id)
        abs_w = math.sqrt(re_w * re_w + im_w * im_w)
        rows.append(
            {
                "tau": tau,
                "nc": nc,
                "corr_id": corr_id,
                "ReW": re_w,
                "ImW": im_w,
                "AbsW": abs_w,
            }
        )

    rows = sorted(rows, key=lambda r: r["nc"])
    all_rows.extend(rows)

    re_vals = [r["ReW"] for r in rows]
    im_vals = [r["ImW"] for r in rows]
    abs_vals = [r["AbsW"] for r in rows]

    n = len(rows)
    re_mean = mean(re_vals)
    im_mean = mean(im_vals)
    abs_mean = mean(abs_vals)

    re_err = stdev(re_vals) / math.sqrt(n) if n > 1 else 0.0
    im_err = stdev(im_vals) / math.sqrt(n) if n > 1 else 0.0
    abs_err = stdev(abs_vals) / math.sqrt(n) if n > 1 else 0.0

    summary_rows.append(
        {
            "tau": tau,
            "Ncfg": n,
            "mean_ReW": re_mean,
            "err_ReW": re_err,
            "mean_ImW": im_mean,
            "err_ImW": im_err,
            "mean_AbsW": abs_mean,
            "err_AbsW": abs_err,
            "filename": str(filename),
        }
    )

summary_rows = sorted(summary_rows, key=lambda r: r["tau"])

# Effective potential from neighboring available tau values:
# V_eff(tau) = log(W(tau) / W(tau_next)) / (tau_next - tau)
for i, row in enumerate(summary_rows):
    if i + 1 < len(summary_rows):
        tau = row["tau"]
        tau_next = summary_rows[i + 1]["tau"]
        w = row["mean_ReW"]
        w_next = summary_rows[i + 1]["mean_ReW"]
        if w > 0 and w_next > 0:
            row["tau_next"] = tau_next
            row["Veff_to_next"] = math.log(w / w_next) / (tau_next - tau)
        else:
            row["tau_next"] = ""
            row["Veff_to_next"] = ""
    else:
        row["tau_next"] = ""
        row["Veff_to_next"] = ""

with open(OUT_DIR / "wilson_tau_per_config.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=["tau", "nc", "corr_id", "ReW", "ImW", "AbsW"]
    )
    writer.writeheader()
    writer.writerows(all_rows)

with open(OUT_DIR / "wilson_tau_summary.csv", "w", newline="") as f:
    fieldnames = [
        "tau",
        "Ncfg",
        "mean_ReW",
        "err_ReW",
        "mean_ImW",
        "err_ImW",
        "mean_AbsW",
        "err_AbsW",
        "tau_next",
        "Veff_to_next",
        "filename",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"Wrote {OUT_DIR / 'wilson_tau_per_config.csv'}")
print(f"Wrote {OUT_DIR / 'wilson_tau_summary.csv'}")

print()
print("Summary:")
for r in summary_rows:
    print(
        f"tau={r['tau']:2d}  "
        f"N={r['Ncfg']:2d}  "
        f"<ReW>={r['mean_ReW']:.8e} +/- {r['err_ReW']:.2e}  "
        f"<ImW>={r['mean_ImW']:.3e}  "
        f"Veff={r['Veff_to_next'] if r['Veff_to_next'] != '' else '---'}"
    )