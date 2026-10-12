"""Paired ratio-of-raw-means analysis; Decimal avoids cancellation loss.

Example: python3 analyze_nine_configurations.py OLD_COLLECTION NEW_COLLECTION OUT
Uses the Python standard library. Plotting is a separate optional command.
"""
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from decimal import Decimal, localcontext
from pathlib import Path

CONFIGS = tuple(range(110, 200, 10))
KEY = ('tau', 'axis', 'R', 'long', 'tr1', 'tr2')
FIELDS = ('L_RE', 'LS_RE', 'S_RAW')


def load(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def means(rows):
    n = Decimal(len(rows))
    return {f: sum((Decimal(r[f]) for r in rows), Decimal(0))/n for f in FIELDS}


def connected(rows):
    m = means(rows)
    if m['L_RE'] == 0:
        raise ValueError('Undefined connected ratio')
    return m['LS_RE']/m['L_RE'] - m['S_RAW']


def jackknife(paired, block_size=1):
    # All channels and both prescriptions lose exactly the same configurations.
    n = len(paired[4])
    assert n == len(paired[40]) and n % block_size == 0
    groups = n//block_size
    assert groups >= 2
    replicates = []
    for group in range(groups):
        indices = set(range(group*block_size, (group+1)*block_size))
        a = connected([r for i,r in enumerate(paired[4]) if i not in indices])
        b = connected([r for i,r in enumerate(paired[40]) if i not in indices])
        replicates.append((a,b,a-b))
    errors = []
    for k in range(3):
        avg = sum((r[k] for r in replicates), Decimal(0))/groups
        variance = Decimal(groups-1)/groups*sum(((r[k]-avg)**2 for r in replicates), Decimal(0))
        errors.append(float(variance.sqrt()))
    return errors


def summarize_rows(rows):
    def maximum(field):
        row = max(rows, key=lambda r: abs(r[field]))
        return dict(value=abs(row[field]), key={k:row[k] for k in KEY},
                    C_p4=row['C_p4'], C_p40=row['C_p40'], delta=row['delta_p4_minus_p40'])
    return dict(points=len(rows), max_absolute_delta=maximum('delta_p4_minus_p40'),
                max_relative_delta=maximum('relative_delta'),
                max_delta_over_p40_unblocked_jk_se=maximum('delta_over_p40_unblocked_jk_se'))


def main(old, new, output):
    old, new, output = Path(old), Path(new), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    raw = load(old/'profiles_raw.csv') + load(new/'profiles_raw.csv')
    windows = load(old/'su3_windows.csv') + load(new/'su3_windows.csv')
    times = load(old/'su3_times.csv') + load(new/'su3_times.csv')
    per_configuration = load(old/'profiles_paired.csv') + load(new/'profiles_paired.csv')
    assert len(raw) == 17928 and len(per_configuration) == 8964
    assert len(windows) == 54 and len(times) == 216
    groups = defaultdict(dict)
    for row in raw:
        cfg, piter = int(row['configuration']), int(row['piter'])
        assert cfg in CONFIGS and piter in (4,40)
        assert int(row['Nv']) == 100 and int(row['tmid']) == 98
        tau = int(row['tau'])
        assert tau in (2,4,6)
        assert (int(row['t0']),int(row['t1'])) == (98-tau//2,98+tau//2)
        key = tuple(int(row[k]) for k in KEY)
        assert int(row['Nsrc']) == 64**3*(3 if key[1] == -1 else 1)
        assert all(math.isfinite(float(row[f])) for f in ('L_RE','L_IM','LS_RE','LS_IM','S_RAW'))
        assert (cfg,piter) not in groups[key], ('Duplicate',key,cfg,piter)
        groups[key][cfg,piter] = row
    assert len(groups) == 996
    all_rows = []
    with localcontext() as context:
        context.prec = 50
        for key, entries in sorted(groups.items()):
            assert set(entries) == {(cfg,p) for cfg in CONFIGS for p in (4,40)}
            paired = {p:[entries[cfg,p] for cfg in CONFIGS] for p in (4,40)}
            for cfg in CONFIGS:
                assert entries[cfg,4]['S_RAW'] == entries[cfg,40]['S_RAW']
            a,b = connected(paired[4]),connected(paired[40])
            delta = a-b
            jk = jackknife(paired)
            block3 = jackknife(paired,3)
            row = dict(zip(KEY,key))
            row.update(nconfigs=9, C_p4=float(a), C_p40=float(b),
                       delta_p4_minus_p40=float(delta),
                       relative_delta=float(abs(delta/b)) if b else None,
                       C_p4_unblocked_jk_se=jk[0], C_p40_unblocked_jk_se=jk[1],
                       paired_delta_unblocked_jk_se=jk[2],
                       C_p4_block3_jk_se=block3[0], C_p40_block3_jk_se=block3[1],
                       paired_delta_block3_jk_se=block3[2],
                       delta_over_p40_unblocked_jk_se=float(abs(delta))/jk[1] if jk[1] else None)
            for p in (4,40):
                row.update({f'mean_{f}_p{p}':float(v) for f,v in means(paired[p]).items()})
            all_rows.append(row)
    assert all(r['relative_delta'] is not None and r['delta_over_p40_unblocked_jk_se'] is not None for r in all_rows)
    # Use the combined-axis measurements already reduced at the raw level.
    # The three axis records are diagnostics, not additional ensemble samples.
    combined = [r for r in all_rows if r['axis'] == -1]
    assert len(combined) == 249
    summary = dict(configurations=list(CONFIGS), nconfigurations=9,
                   taus=[2,4,6], midpoint=98, separations=[1,2,4],
                   HYP_steps=1, HYP_alpha=[1,1,0.5], piters=[4,40],
                   raw_rows=len(raw), per_configuration_paired_rows=len(per_configuration),
                   averaged_rows=len(all_rows), combined_axis_averaged_rows=len(combined),
                   estimator='mean(LS_RE)/mean(L_RE) - mean(S_RAW), equal configuration weights',
                   uncertainty='Paired configuration jackknife; block_size=1 and exploratory 3 consecutive configurations per block. Nine ordered configurations cannot establish autocorrelation or final statistical errors.',
                   convention='Original-field clover F=(Q-Q_dagger)/8, sum over 6 planes Tr(F_dagger F); no color-trace subtraction or extra factor 1/2. Raw lattice units.',
                   combined_axes=summarize_rows(combined),
                   all_axes_diagnostic=summarize_rows(all_rows),
                   by_tau={str(t):summarize_rows([r for r in combined if r['tau']==t]) for t in (2,4,6)},
                   by_tau_R={f'{t},{r}':summarize_rows([x for x in combined if (x['tau'],x['R'])==(t,r)]) for t in (2,4,6) for r in (1,2,4)})
    p4windows = [r for r in windows if int(r['piter'])==4]
    p40windows = [r for r in windows if int(r['piter'])==40]
    def worst(rows,field):
        r=max(rows,key=lambda r:float(r[field]))
        return dict(configuration=int(r['configuration']),tau=int(r['tau']),value=float(r[field]))
    summary['SU3'] = dict(p4_unit_max=worst(p4windows,'unit_max'),
                          p4_det_max=worst(p4windows,'det_max'),
                          p4_unit_rms_max=worst(p4windows,'unit_rms'),
                          p4_det_rms_max=worst(p4windows,'det_rms'),
                          p40_unit_max=worst(p40windows,'unit_max'),
                          p40_det_max=worst(p40windows,'det_max'),
                          p4_windows_with_unit_or_det_above_1em2=sum(int(r['unit_above_1em2'])>0 or int(r['det_above_1em2'])>0 for r in p4windows),
                          scope='Outgoing temporal links t=95..100 at midpoint 98, with nested tau windows. Counts across windows must not be added as independent links; no boundary/global acceptance claim.')
    summary['single_configuration_max_relative_C_combined'] = max(float(r['C_relative_difference']) for r in per_configuration if int(r['axis'])==-1)
    summary['input_sha256'] = {f'{label}/{name}':hashlib.sha256((p/name).read_bytes()).hexdigest()
                             for label,p in [('campaign_5918419',old),('campaign_5920223',new)]
                             for name in ('profiles_raw.csv','profiles_paired.csv','su3_windows.csv','su3_times.csv','SUMMARY.json')}
    write_csv(output/'ensemble_profiles.csv',all_rows)
    write_csv(output/'ensemble_profiles_combined_axes.csv',combined)
    write_csv(output/'raw_profiles_n110_190.csv',raw)
    write_csv(output/'su3_windows_n110_190.csv',windows)
    write_csv(output/'su3_times_n110_190.csv',times)
    (output/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k] for k in ('configurations','combined_axes','by_tau','SU3')},indent=2))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    main(*sys.argv[1:])
