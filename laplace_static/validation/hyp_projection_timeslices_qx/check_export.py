#!/usr/bin/env python3
"""Validate and summarize the time-resolved physical temporal-link diagnostics."""
import csv
import json
import math
import re
import sys
from pathlib import Path

STAGES = ('ORIGINAL', 'PITER4', 'PITER40')
METRICS = ('UNIT_MAX_ENTRY', 'UNIT_MAX_FROB', 'UNIT_RMS_FROB', 'DET_MAX', 'DET_RMS')
THRESHOLDS = (1e-12, 1e-8, 1e-4, 1e-2)
COUNTS = tuple(f'{family}_COUNT_{k}' for k in range(4) for family in ('UNIT', 'DET'))


def fields(line):
    return dict(re.findall(r'(\w+)=([^\s]+)', line))


def near(a, b, label):
    if not math.isclose(a, b, rel_tol=2e-12, abs_tol=0.0):
        raise ValueError('Global/per-time inconsistency: '+label)


def read(text):
    by_time, globals_, exceed, differences = {}, {}, {}, {}
    global_difference = None
    for line in text.splitlines():
        if line.startswith('SU3_TIME '):
            stage = line.split()[1]
            values = fields(line)
            t = int(values['T'])
            if stage not in STAGES or (stage, t) in by_time:
                raise ValueError('Unknown or duplicate per-time row')
            row = {'stage': stage, 't': t, 'N': int(values['N'])}
            row.update({k: float(values[k]) for k in METRICS})
            row.update({k: int(values[k]) for k in COUNTS})
            if row['N'] != 64**3:
                raise ValueError('Incorrect per-time physical link count')
            if any(not math.isfinite(row[k]) or row[k] < 0 for k in METRICS):
                raise ValueError('Invalid per-time SU(3) diagnostic')
            for family in ('UNIT', 'DET'):
                counts = [row[f'{family}_COUNT_{k}'] for k in range(4)]
                if any(not 0 <= c <= row['N'] for c in counts) or counts != sorted(counts, reverse=True):
                    raise ValueError('Invalid exceedance counts')
            by_time[(stage, t)] = row
        elif line.startswith('SU3 '):
            stage = line.split()[1]
            if stage not in STAGES or stage in globals_:
                raise ValueError('Duplicate or unknown global SU(3) summary')
            globals_[stage] = fields(line)
        elif line.startswith('SU3_EXCEED '):
            stage = line.split()[1]
            values = fields(line)
            threshold = float(values['threshold'])
            k = THRESHOLDS.index(threshold)
            if (stage, k) in exceed:
                raise ValueError('Duplicate global exceedance summary')
            exceed[(stage, k)] = values
        elif line.startswith('PITER_TIME '):
            values = fields(line)
            t, value = int(values['T']), float(values['COMPONENT_MAX'])
            if t in differences or not math.isfinite(value) or value < 0:
                raise ValueError('Invalid or duplicate per-time difference')
            differences[t] = value
        elif line.startswith('PITER_4_VS_40 '):
            if global_difference is not None:
                raise ValueError('Duplicate global difference')
            global_difference = float(fields(line)['temporal'])
    expected = {(stage, t) for stage in STAGES for t in range(191)}
    if set(by_time) != expected or set(differences) != set(range(191)):
        raise ValueError('Missing time slices, unexpected upper-boundary record, or wrong stage set')
    if set(globals_) != set(STAGES) or len(exceed) != 12:
        raise ValueError('Incomplete global diagnostics')
    for stage in STAGES:
        rows = [by_time[(stage, t)] for t in range(191)]
        n = sum(r['N'] for r in rows)
        if int(globals_[stage]['N']) != n:
            raise ValueError('Global physical-link count mismatch')
        for metric in ('UNIT_MAX_ENTRY', 'UNIT_MAX_FROB', 'DET_MAX'):
            near(max(r[metric] for r in rows), float(globals_[stage][metric]), stage+' '+metric)
        for metric in ('UNIT_RMS_FROB', 'DET_RMS'):
            reconstructed = math.sqrt(sum(r['N'] * r[metric]**2 for r in rows)/n)
            near(reconstructed, float(globals_[stage][metric]), stage+' '+metric)
        for k in range(4):
            for family in ('UNIT', 'DET'):
                if sum(r[f'{family}_COUNT_{k}'] for r in rows) != int(exceed[(stage,k)][family+'_COUNT']):
                    raise ValueError('Global/per-time exceedance count mismatch')
    if global_difference is None or not math.isfinite(global_difference):
        raise ValueError('Missing or invalid global difference')
    near(max(differences.values()), global_difference, 'piter4 versus 40')
    return by_time, differences


def summarize(by_time, differences, times, stage):
    rows = [by_time[(stage, t)] for t in times]
    n = sum(r['N'] for r in rows)
    out = {'time_start': min(times), 'time_end_inclusive': max(times), 'N': n}
    for metric in ('UNIT_MAX_ENTRY', 'UNIT_MAX_FROB', 'DET_MAX'):
        peak = max(rows, key=lambda r: r[metric])
        out[metric], out[metric+'_time'] = peak[metric], peak['t']
    for metric in ('UNIT_RMS_FROB', 'DET_RMS'):
        out[metric] = math.sqrt(sum(r['N']*r[metric]**2 for r in rows)/n)
    for name in COUNTS:
        out[name] = sum(r[name] for r in rows)
    peak_t = max(times, key=lambda t: differences[t])
    out['PITER_4_VS_40_COMPONENT_MAX'] = differences[peak_t]
    out['PITER_4_VS_40_COMPONENT_MAX_time'] = peak_t
    return out


def export(by_time, differences, destination):
    windows = {'full': list(range(191)), 'tau2': [97,98],
               'tau4': list(range(96,100)), 'tau6': list(range(95,101))}
    summary = {
        'configuration': 'qX_p2n100', 'grid': [6,4,4,4], 'alpha': [1,1,0.5],
        'hyp_steps': 1, 'excluded': 'Absent outgoing temporal link at t=191',
        'scope': 'Final physical temporal links; per-time reduction/inventory checks; no ensemble accuracy claim',
        'threshold_index': dict(enumerate(THRESHOLDS)),
        'windows': {name: {stage: summarize(by_time,differences,times,stage) for stage in STAGES}
                    for name,times in windows.items()},
        'piter4_times_above_1e-2': [
            {'t': t, 'UNIT_COUNT': by_time[('PITER4',t)]['UNIT_COUNT_3'],
             'DET_COUNT': by_time[('PITER4',t)]['DET_COUNT_3']}
            for t in range(191) if by_time[('PITER4',t)]['UNIT_COUNT_3']
            or by_time[('PITER4',t)]['DET_COUNT_3']],
    }
    names = ('su3_by_time.csv', 'piter_difference_by_time.csv', 'SUMMARY.json', 'SUMMARY.md')
    if any((destination/name).exists() for name in names):
        raise ValueError('Derived output already exists; preserving it')
    with (destination/names[0]).open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['stage','t','N',*METRICS,*COUNTS])
        writer.writeheader()
        writer.writerows(by_time[(stage,t)] for stage in STAGES for t in range(191))
    with (destination/names[1]).open('x', newline='') as f:
        writer = csv.writer(f); writer.writerow(['t','COMPONENT_MAX'])
        writer.writerows((t, differences[t]) for t in range(191))
    (destination/names[2]).write_text(json.dumps(summary, indent=2)+'\n')
    lines = ['# Final temporal-link projection by time on qX_p2n100', '',
             'One temporal HYP step, alpha=(1,1,0.5). Physical temporal links only; t=191 excluded.', '',
             '| Link window | piter | Max unitarity Frobenius defect | Max determinant defect | Links with unitarity defect > 0.01 | Max p4-vs-p40 component difference |',
             '|---|---:|---:|---:|---:|---:|']
    for name in windows:
        for stage in ('PITER4','PITER40'):
            s = summary['windows'][name][stage]
            lines.append('| {} (t={}..{}) | {} | {:.8e} | {:.8e} | {} / {} | {:.8e} |'.format(
                name,s['time_start'],s['time_end_inclusive'],stage[5:],s['UNIT_MAX_FROB'],
                s['DET_MAX'],s['UNIT_COUNT_3'],s['N'],s['PITER_4_VS_40_COMPONENT_MAX']))
    lines += ['', 'Counts, maxima and RMS defects reconstruct the separate global diagnostics.',
              'This localizes the measured tails. It does not independently validate HYP, intermediate decorated links, or ensemble observables.', '']
    (destination/names[3]).write_text('\n'.join(lines))
    print('PASS 573 per-time SU(3) records, 191 per-time differences and global reconstruction')
    for name in ('full','tau2','tau4','tau6'):
        s = summary['windows'][name]['PITER4']
        print('PITER4_WINDOW', name, 'T={}..{}'.format(s['time_start'],s['time_end_inclusive']),
              'UNIT_MAX_FROB={:.16e}'.format(s['UNIT_MAX_FROB']),
              'DET_MAX={:.16e}'.format(s['DET_MAX']),
              'UNIT_GT_1E_MINUS2={}/{}'.format(s['UNIT_COUNT_3'],s['N']),
              'P4_VS_40={:.16e}'.format(s['PITER_4_VS_40_COMPONENT_MAX']))
    print('CSV:',destination/names[0])
    print('CSV:',destination/names[1])
    print('SUMMARY:',destination/names[3])


if __name__ == '__main__':
    by_time, differences = read(Path(sys.argv[1]).read_text())
    export(by_time, differences, Path(sys.argv[2]))
