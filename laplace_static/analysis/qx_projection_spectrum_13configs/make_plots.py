#!/usr/bin/env python3
"""Reproduce the matched qX piter=4,10,20,40 analysis from preserved exports."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PITERS = (4, 10, 20, 40)
KEEP = (
    'profiles_all_piters_13configs.png',
    'profiles_all_piters_13configs.svg',
    'su3_central_windows_all_piters_13configs.png',
    'su3_central_windows_all_piters_13configs.svg',
    'profiles_ensemble_13configs.csv',
    'su3_central_windows_13configs.csv',
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def verify_inputs():
    sources = json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name, expected in sources['preserved_inputs'].items():
        if sha256(HERE/'inputs'/name) != expected:
            raise ValueError('Preserved input changed: '+name)
    for name, expected in sources['upstream_files'].items():
        if sha256(ROOT/name) != expected:
            raise ValueError('Upstream numerical source changed: '+name)
    validation = json.loads((HERE/'inputs/job5926668_validation.json').read_text())
    if not validation['success'] or (validation['cases'], validation['mpi_launches']) != (39, 78):
        raise ValueError('Incomplete intermediate-iteration validation')
    return sources


def comparison_summary(profiles, windows):
    combined = [row for row in profiles if int(row['axis']) == -1]
    result = []
    with localcontext() as context:
        context.prec = 50
        for p in PITERS:
            differences = [(abs(Decimal(r['C_p'+str(p)])-Decimal(r['C_p40'])), r)
                           for r in combined]
            maximum, location = max(differences, key=lambda value: value[0])
            relative = max(value/abs(Decimal(row['C_p40']))*100
                           for value, row in differences)
            selected = [r for r in windows if int(r['piter']) == p]
            unit = max(selected, key=lambda r: Decimal(r['unit_max']))
            det = max(selected, key=lambda r: Decimal(r['det_max']))
            result.append(dict(
                piter=p, configurations=13, combined_profile_positions=len(combined),
                max_abs_delta_C_vs_p40=str(maximum),
                max_relative_delta_percent_vs_p40=str(relative),
                max_delta_tau=location['tau'], max_delta_R=location['R'],
                max_delta_long=location['long'], max_delta_tr1=location['tr1'],
                max_delta_tr2=location['tr2'],
                largest_window_unit_max=unit['unit_max'],
                unit_max_configuration=unit['configuration'], unit_max_tau=unit['tau'],
                largest_window_det_max=det['det_max'],
                det_max_configuration=det['configuration'], det_max_tau=det['tau']))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'results')
    args = parser.parse_args()
    sources = verify_inputs()
    args.output.mkdir(parents=True, exist_ok=True)
    generator = ROOT/'laplace_static/validation/qx_projection_spectrum/plot_exports.py'
    with tempfile.TemporaryDirectory(prefix='qx_spectrum13_') as temporary:
        generated = Path(temporary)
        subprocess.run([
            sys.executable, str(generator),
            '--extra-raw', str(HERE/'inputs/job5926668_profiles_raw.csv'),
            '--extra-windows', str(HERE/'inputs/job5926668_su3_windows.csv'),
            '--output', str(generated)], check=True)
        summary = json.loads((generated/'SUMMARY.json').read_text())
        if summary['configurations'] != list(range(110, 231, 10)) or summary['profile_piters'] != list(PITERS):
            raise ValueError('Changed matched profile cohort')
        for name in KEEP:
            shutil.copyfile(generated/name, args.output/name)
    comparisons = comparison_summary(
        rows(args.output/'profiles_ensemble_13configs.csv'),
        rows(args.output/'su3_central_windows_13configs.csv'))
    with (args.output/'piter_comparison_13configs.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(comparisons)
    summary.update(
        completed_scan_job=5926668,
        profile_prescriptions=list(PITERS),
        piter_comparison=comparisons,
        intermediate_launches=78,
        intermediate_raw_records=25896,
        intermediate_window_records=78,
        total_window_records=156,
        scope='13 matched configurations and three central time windows. Raw clover convention. '
              'The p4/p40 jackknife errors remain preliminary and are not autocorrelation corrected. '
              'Intermediate-piter profile plots have no final statistical error bars.',
        preserved_inputs=sources['preserved_inputs'],
        review_archive_sha256=sources['review_archive_sha256'])
    (args.output/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('PASS matched four-piter results:', args.output)


if __name__ == '__main__':
    main()

