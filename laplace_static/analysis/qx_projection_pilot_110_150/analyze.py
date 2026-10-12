#!/usr/bin/env python3
"""Recompute paired piter diagnostics from the campaign's unchanged raw CSVs.

Usage: python3 analyze.py REVIEW.zip --output NEW_DIRECTORY [--plot]
The connected quantity uses configuration means before ratio/subtraction.
Only the Python standard library is required unless --plot is supplied.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from decimal import Decimal, getcontext
from pathlib import Path
import zipfile

getcontext().prec = 50
D = Decimal
CONFIGS = (110, 120, 130, 140, 150)
KEY_FIELDS = ("tau", "axis", "R", "long", "tr1", "tr2")
RAW_FIELDS = ("L_RE", "L_IM", "S_RAW", "LS_RE", "LS_IM")
PREFIX = "laplace_static/validation/hyp_projection_configs_qx_queue/collection_mixed_5918419/"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def key(row):
    return tuple(int(row[k]) for k in KEY_FIELDS)


def mean(values):
    return sum(values, D(0)) / D(len(values))


def connected(rows):
    values = {f: mean([D(r[f]) for r in rows]) for f in RAW_FIELDS}
    require(values["L_RE"] != 0, "Zero mean real loop denominator")
    return values, values["LS_RE"] / values["L_RE"] - values["S_RAW"]


def jk_se(values):
    # Descriptive only: five adjacent configurations, no autocorrelation blocking.
    centre = mean(values)
    return ((D(len(values) - 1) / D(len(values))) *
            sum(((v - centre) ** 2 for v in values), D(0))).sqrt()


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def maximum(rows, field, absolute=False):
    row = max(rows, key=lambda r: abs(D(r[field])) if absolute else D(r[field]))
    return {"value": str(abs(D(row[field])) if absolute else D(row[field])),
            "point": {k: row[k] for k in KEY_FIELDS}}


def analyze(archive):
    with zipfile.ZipFile(archive) as z:
        require(z.testzip() is None, "ZIP CRC failure")
        read_csv = lambda name: list(csv.DictReader(io.StringIO(z.read(PREFIX + name).decode())))
        raw = read_csv("profiles_raw.csv")
        paired = read_csv("profiles_paired.csv")
        windows = read_csv("su3_windows.csv")
        campaign = json.loads(z.read(PREFIX + "SUMMARY.json"))
    require(campaign["success"] and not campaign["errors"], "Campaign did not succeed")
    require(campaign["validated_cases"] == 15, "Expected 15 cases")
    require(len(raw) == 9960 and len(paired) == 4980 and len(windows) == 30,
            "Unexpected campaign row counts")
    require({(c["configuration"], c["tau"]) for c in campaign["cases"]} ==
            {(c, t) for c in CONFIGS for t in (2, 4, 6)}, "Incomplete case set")

    indexed = {}
    points = defaultdict(dict)
    for row in raw:
        cfg, piter = int(row["configuration"]), int(row["piter"])
        point = key(row)
        require(cfg in CONFIGS and piter in (4, 40), "Unexpected configuration or piter")
        require(int(row["Nv"]) == 100 and int(row["tmid"]) == 98, "Changed measurement setup")
        require(int(row["t0"]) == 98 - point[0] // 2 and
                int(row["t1"]) == 98 + point[0] // 2, "Changed temporal window")
        require(all(D(row[f]).is_finite() for f in RAW_FIELDS), "Nonfinite raw value")
        token = (cfg, piter) + point
        require(token not in indexed, "Duplicate raw record")
        indexed[token] = row
        points[point][(cfg, piter)] = row
    require(len(points) == 996, "Unexpected profile geometry")

    # Validate exported per-configuration diagnostics against raw records.
    for row in paired:
        cfg = int(row["configuration"])
        point = key(row)
        values = []
        for p in (4, 40):
            source = indexed[(cfg, p) + point]
            for field in RAW_FIELDS:
                require(D(row[f"{field}_p{p}"]) == D(source[field]), "Paired raw column mismatch")
            _, c = connected([source])
            require(abs(c - D(row[f"C_p{p}"])) < D("1e-42"), "Paired connected export mismatch")
            values.append(c)
        require(abs(values[0] - values[1] - D(row["C_difference"])) < D("1e-42"),
                "Paired difference mismatch")

    ensemble = []
    for point, records in sorted(points.items()):
        require(set(records) == {(c, p) for c in CONFIGS for p in (4, 40)},
                "A point lacks a configuration or prescription")
        for cfg in CONFIGS:
            require(records[(cfg, 4)]["S_RAW"] == records[(cfg, 40)]["S_RAW"],
                    "Original-field clover changed between prescriptions")
        result = dict(zip(KEY_FIELDS, point))
        result["configurations"] = len(CONFIGS)
        cs, replicates = {}, {}
        for p in (4, 40):
            values, cs[p] = connected([records[(c, p)] for c in CONFIGS])
            for field in RAW_FIELDS:
                result[f"mean_{field}_p{p}"] = str(values[field])
            result[f"C_p{p}"] = str(cs[p])
            replicates[p] = [connected([records[(c, p)] for c in CONFIGS if c != omit])[1]
                             for omit in CONFIGS]
        delta = cs[4] - cs[40]
        require(cs[40] != 0, "Zero reference connected diagnostic")
        result["delta_C_p4_minus_p40"] = str(delta)
        result["absolute_relative_delta_C"] = str(abs(delta / cs[40]))
        result["unblocked_jackknife_se_C_p40"] = str(jk_se(replicates[40]))
        result["unblocked_paired_jackknife_se_delta_C"] = str(jk_se(
            [a - b for a, b in zip(replicates[4], replicates[40])]))
        # Independent double-precision calculation detects gross ratio/averaging mistakes.
        f40 = [records[(c, 40)] for c in CONFIGS]
        independent = (math.fsum(float(r["LS_RE"]) for r in f40) /
                       math.fsum(float(r["L_RE"]) for r in f40) -
                       math.fsum(float(r["S_RAW"]) for r in f40) / len(f40))
        require(abs(independent - float(cs[40])) < 3e-15, "Independent float mean calculation disagrees")
        ensemble.append(result)

    case_metrics = []
    for case in campaign["cases"]:
        cfg, tau = case["configuration"], case["tau"]
        p4 = [r for r in paired if int(r["configuration"]) == cfg and int(r["tau"]) == tau]
        require(len(p4) == 332, "Incomplete case profile")
        case_metrics.append({
            "configuration": cfg, "tau": tau, "job": case["array_job_id"],
            "p4_unit_max": case["p4_window"]["unit_max"],
            "p4_unit_rms": case["p4_window"]["unit_rms"],
            "p4_det_max": case["p4_window"]["det_max"],
            "p4_det_rms": case["p4_window"]["det_rms"],
            "p40_unit_max": case["p40_window"]["unit_max"],
            "p40_det_max": case["p40_window"]["det_max"],
            "max_abs_C_delta_all_axes": str(max(abs(D(r["C_difference"])) for r in p4)),
            "max_relative_C_delta_combined": str(max(abs(D(r["C_relative_difference"]))
                for r in p4 if int(r["axis"]) == -1)),
            "mpi_seconds": case["mpi_seconds"],
        })
    combined = [r for r in ensemble if r["axis"] == -1]
    w4 = [r for r in windows if int(r["piter"]) == 4]
    w40 = [r for r in windows if int(r["piter"]) == 40]
    summary = {
        "source": {"archive": archive.name, "sha256": hashlib.sha256(archive.read_bytes()).hexdigest()},
        "configurations": list(CONFIGS), "tau": [2, 4, 6], "midpoint": 98,
        "raw_rows": len(raw), "paired_rows": len(paired), "cases": 15,
        "mean_profile_points": len(ensemble), "combined_axis_points": len(combined),
        "p4_max_window_unit_frobenius": str(max(D(r["unit_max"]) for r in w4)),
        "p4_max_window_unit_rms": str(max(D(r["unit_rms"]) for r in w4)),
        "p4_max_window_determinant_defect": str(max(D(r["det_max"]) for r in w4)),
        "p4_window_unit_links_above_1e_minus2": sum(int(r["unit_above_1em2"]) for r in w4),
        "p4_window_det_links_above_1e_minus2": sum(int(r["det_above_1em2"]) for r in w4),
        "window_count_scope": "Window counters overlap between tau values; not counts of unique ensemble links",
        "p40_max_window_unit_or_det": str(max(D(r[f]) for r in w40 for f in ("unit_max", "det_max"))),
        "max_per_config_relative_C_combined": max(c["max_relative_difference"]["C_combined"]
                                                   for c in campaign["cases"]),
        "max_per_config_relative_C_all_axes": max(c["max_relative_difference"]["C_all_axes"]
                                                   for c in campaign["cases"]),
        "mean_max_abs_delta_combined": maximum(combined, "delta_C_p4_minus_p40", True),
        "mean_max_relative_delta_combined": maximum(combined, "absolute_relative_delta_C"),
        "mean_max_abs_delta_all_axes": maximum(ensemble, "delta_C_p4_minus_p40", True),
        "mean_max_relative_delta_all_axes": maximum(ensemble, "absolute_relative_delta_C"),
        "estimator": "mean(LS_RE)/mean(L_RE)-mean(S_RAW); equal weight per configuration",
        "uncertainty_scope": "Exploratory delete-one-configuration jackknife on five successive configurations; no autocorrelation correction or final ensemble precision claim",
        "recommendation": "piter4 working prescription for tested central-window pilot measurements; retain piter40 controls and recheck new configurations/windows",
        "excluded_scope": "Open-boundary regions, other windows/configurations, independent HYP/operator implementation, physical normalization, full saved-mode residual coverage",
    }
    return ensemble, case_metrics, summary


def report(summary, cases):
    rows = []
    for cfg in CONFIGS:
        group = [c for c in cases if c["configuration"] == cfg]
        rows.append("| {} | {:.3e} | {:.3e} | {:.6f}% |".format(
            cfg, max(float(c["p4_unit_max"]) for c in group),
            max(float(c["p4_det_max"]) for c in group),
            100 * max(float(c["max_relative_C_delta_combined"]) for c in group)))
    absolute = float(summary["mean_max_abs_delta_combined"]["value"])
    relative = 100 * float(summary["mean_max_relative_delta_combined"]["value"])
    return """# qX central-window projection pilot: n110–n150

All 15 cases completed and passed structural/profile checks, unchanged original-field clover checks, and the piter40 SU(3) control. The n100 tau4 regression also reproduced its previous outputs exactly. These checks report piter4 defects without imposing a physics acceptance threshold.

## Setup and conventions

Five configurations (110,120,130,140,150), tau=2,4,6 with midpoint98, separations R=1,2,4, Nv=100, rho_i=1. One temporal HYP step with alpha=(1,1,0.5), compared at piter=4 and40. Projection iterations are separate from HYP smearing steps and from the spatial APE20,alpha0.5 Laplace-eigenvector operator.

The clover probe uses the original gauge field: F=(Q-Q†)/8 and S_RAW=sum over the six mu<nu planes of Tr(F†F), with no color-trace subtraction and no extra1/2. All reported connected quantities are dimensionless diagnostics in this raw convention, without a physical normalization or continuum interpretation.

Physical temporal links in the actual transport windows t0..t1-1 are audited. The open-boundary regions are outside this central-window claim. Observable jobs use grid1x4x4x4,64MPI ranks; this is distinct from the eigenvector-production grid6x4x4x4. FUNNELED initialization is retained.

## Projection diagnostics

Maxima over the three tau values for each configuration:

| Configuration | p4 max Frobenius unitarity defect | p4 max determinant defect | Max per-configuration combined-profile relative change |
|---|---:|---:|---:|
""" + "\n".join(rows) + """

The worst p4 window RMS unitarity defect is9.227e-7. No audited central-window links exceed1e-2 in either defect. The p40 control's largest unitarity/determinant defect is7.792e-16. A 1e-2 counter is a diagnostic threshold, not an established physics acceptance criterion.

The largest individual-axis relative connected change is5.0795%, at a small reference connected value (5.5601e-7); its absolute change is2.8243e-8. Small denominators can inflate these percentages. The largest per-configuration absolute connected change across all axes is8.3494e-8.

## Five-configuration raw-average comparison

For each matched profile point and prescription, compute

    C = mean_configuration(LS_RE)/mean_configuration(L_RE) - mean_configuration(S_RAW).

Each configuration has equal weight; spatial sites, axes and tau values are not extra independent configurations. This calculation does not average the per-configuration ratios.

The largest combined-axis absolute change after raw averaging is **{absolute:.6e}**, and the largest combined-axis relative change is **{relative:.6f}%**. Both maxima occur at tau6,R4,long5,tr1=2,tr2=0. There are996 matched mean-profile points, including249 combined-axis points. Largest all-axis absolute change is3.0020e-8; largest all-axis relative change is0.139968%.

`ensemble_profiles.csv` retains the raw means, both connected diagnostics, their paired difference, and exploratory delete-one-configuration jackknife scales. The same configuration is omitted from both prescriptions for each paired replicate. These five successive configurations have not been assessed for autocorrelation; the jackknife columns are not final errors or evidence of ensemble-wide significance. The plotted example cut is illustrative; the maxima above cover all exported points.

## Working decision and next work

The small measured profile sensitivity supports using piter4 as a working central-window pilot prescription, with piter40 retained as a control. This does not establish exact SU(3) projection or ensemble-wide acceptance. Recheck defects and paired observable sensitivity when adding configurations or changing windows; previous large open-boundary defects are not covered by this result.

Proceed with the next four eigenvector configurations n160,n170,n180,n190 using the already successful packed32-node launcher. Temporal HYP piter does not enter that eigenvector operator. After completion, audit saved residuals/normalization/Gram matrices on selected slices and extend the paired profile pilot. Existing saved-mode quality checks cover nine slices of n110–n150; tau2/6 endpoints are outside the independently audited slice set.

## Evidence and reproducibility

Source: `{source}`; SHA256 `{sha}`. The uploaded archive was CRC checked,26 archived source/checker hashes were verified, and the original15-case checkers were re-run on copied stdout. Re-exported individual raw/paired/SU(3) CSVs matched their uploaded counterparts byte for byte. This replay verifies archived numerical checks; it does not re-run MPI gauge computations or independently validate the absent libraries, gauge fields or eigenvectors.

`analyze.py` separately validates coverage, paired export arithmetic, unchanged S_RAW, and raw-before-ratio ensemble arithmetic using50-digit Decimal calculations checked against an independent double-precision mean. Original inputs are never modified. Reproduce into a new directory:

    python3 analyze.py hyp_projection_configs_review_mixed_5918419.zip --output reproduced --plot

The optional plot needs matplotlib. Other analysis uses Python's standard library. No HDF5/eigenvector data are included in this package.
""".format(absolute=absolute, relative=relative, source=summary["source"]["archive"],
           sha=summary["source"]["sha256"])


def plot(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 6.4), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    for tau, color in zip((2, 4, 6), ("#0072B2", "#D55E00", "#009E73")):
        cut = sorted([r for r in rows if (r["tau"], r["axis"], r["R"], r["long"], r["tr2"])
                      == (tau, -1, 2, 1, 0)], key=lambda r: r["tr1"])
        x = [r["tr1"] for r in cut]
        axes[0].errorbar(x, [1e4 * float(r["C_p40"]) for r in cut],
                        yerr=[1e4 * float(r["unblocked_jackknife_se_C_p40"]) for r in cut],
                        color=color, fmt="o-", capsize=3, label=f"tau={tau}, p40")
        axes[0].plot(x, [1e4 * float(r["C_p4"]) for r in cut],
                     color=color, linestyle="--", marker="x", label=f"tau={tau}, p4")
        axes[1].plot(x, [1e9 * float(r["delta_C_p4_minus_p40"]) for r in cut],
                     color=color, marker="o", label=f"tau={tau}")
    axes[0].set_ylabel("Raw connected diagnostic C (×10⁻⁴)")
    axes[1].set_ylabel("C(p4) − C(p40) (×10⁻⁹)")
    axes[1].set_xlabel("Transverse displacement tr1 (lattice sites)")
    axes[1].axhline(0, color="#555555", linewidth=.7)
    axes[0].legend(ncol=3, fontsize=9)
    axes[1].legend(ncol=3, fontsize=9)
    for ax in axes:
        ax.grid(alpha=.2)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("qX n110–n150: raw means before ratio\nCombined axes; R=2, long=1, tr2=0; midpoint=98", fontsize=12)
    fig.text(.12, .018, "Five configurations only. Bars: unblocked jackknife; autocorrelation unassessed.\nOne HYP step, alpha=(1,1,0.5). Raw clover convention; no physical normalization.", fontsize=9)
    fig.tight_layout(rect=(0, .085, 1, .92))
    fig.savefig(output / "projection_pilot_profiles.png", dpi=180)
    fig.savefig(output / "projection_pilot_profiles.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    require(not args.output.exists(), "Output exists; preserving it")
    rows, cases, summary = analyze(args.archive)
    args.output.mkdir(parents=True)
    write_csv(args.output / "ensemble_profiles.csv", rows)
    write_csv(args.output / "case_metrics.csv", cases)
    (args.output / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output / "REPORT.md").write_text(report(summary, cases))
    if args.plot:
        plot(rows, args.output)
    print(json.dumps({"success": True, "points": len(rows),
                      "combined_max_absolute_delta": summary["mean_max_abs_delta_combined"]}, indent=2))


if __name__ == "__main__":
    main()
