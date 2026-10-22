#!/usr/bin/env python3
"""Reparse job5942805 and reproduce its single-configuration clover comparison."""
import argparse
import csv
import hashlib
import importlib.util
import json
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INPUT = HERE/'inputs/job5942805'
D = Decimal
DERIVED = ['C_RAW', 'C_TRACE_COMPONENT', 'C_TL']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


pilot = module('qx_trace_pilot_analysis', ROOT/'laplace_static/validation/qx_clover_trace/campaign.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def derive(row):
    values = pilot.connected_from_primaries(row)
    return dict(row, **{k: str(values[k]) for k in DERIVED})


def validated_rows():
    sources = json.loads((HERE/'inputs/SOURCES.json').read_text())
    for name, expected in sources['files'].items():
        if sha(HERE/name) != expected:
            raise ValueError('Pilot input changed: '+name)
    pilot.verify_sources()
    if (INPUT/'source_INPUTS.json').read_bytes() != (pilot.HERE/'INPUTS.json').read_bytes():
        raise ValueError('Pilot source pins differ from preserved implementation')
    prepared = json.loads((INPUT/'PREPARED.json').read_text())
    if sha(INPUT/'PREPARED.json') != (INPUT/'prepared.sha256').read_text().strip():
        raise ValueError('Pilot preparation digest differs')
    if sha(INPUT/'build_manifest.sha256') != prepared['build_manifest_sha256']:
        raise ValueError('Pilot build manifest differs')
    launches = json.loads((INPUT/'LAUNCHES.json').read_text())
    if (prepared['cases'] != pilot.case_plan() or
            [(r['configuration'], r['tau'], r['piter'], r['exit_code'], r['job_id']) for r in launches] !=
            [(100, 4, p, 0, '5942805') for p in (4, 40)]):
        raise ValueError('Pilot launch inventory differs')
    check = pilot.verify_sources().checker()
    scans, result = {}, []
    for p in (4, 40):
        path = INPUT/'outputs'/('p%d.out' % p)
        scans[p] = check.read_measurement(path, 100, 4, p)[1]
        trace = pilot.trace_records(path.read_text(), scans[p], check.avg)
        for row in check.raw_rows(100, 4, p, scans[p]):
            key = tuple(int(row[k]) for k in ('axis', 'R', 'long', 'tr1', 'tr2'))
            row.update({k: trace[key][k] for k in pilot.TRACE_FIELDS})
            result.append({k: str(v) for k, v in row.items()})
    check.regression_check(scans)
    archived = read_csv(INPUT/'profiles_raw_and_trace.csv')
    if result != archived or len(result) != 664:
        raise ValueError('Reparsed primary table differs')
    by_key = {(int(r['piter']), *(int(r[k]) for k in ('axis', 'R', 'long', 'tr1', 'tr2'))): r for r in result}
    for key, row in by_key.items():
        if key[0] == 4:
            other = by_key[(40, *key[1:])]
            for field in ('S_RAW', 'T_TRACE'):
                if D(row[field]) != D(other[field]):
                    raise ValueError('Original clover changed across prescriptions')
    return [derive(r) for r in result], launches


def summary(rows, launches):
    combined = [r for r in rows if int(r['axis']) == -1]
    comparisons = []
    with localcontext() as context:
        context.prec = 50
        for p in (4, 40):
            selected = [r for r in combined if int(r['piter']) == p]
            maximum = max(selected, key=lambda r: abs(D(r['C_TRACE_COMPONENT'])))
            relative = max((r for r in selected if D(r['C_RAW']) != 0),
                           key=lambda r: abs(D(r['C_TRACE_COMPONENT'])/D(r['C_RAW'])))
            midpoint = next(r for r in selected if tuple(int(r[k]) for k in ('R', 'long', 'tr1', 'tr2')) == (4, 2, 0, 0))
            comparisons.append(dict(piter=p, combined_positions=len(selected),
                                    mean_S_RAW=selected[0]['S_RAW'], mean_T_TRACE=selected[0]['T_TRACE'],
                                    maximum_absolute_connected_trace_component=str(abs(D(maximum['C_TRACE_COMPONENT']))),
                                    maximum_relative_change_percent=str(100*abs(D(relative['C_TRACE_COMPONENT'])/D(relative['C_RAW']))),
                                    maximum_relative_change_geometry={k: int(relative[k]) for k in ('R', 'long', 'tr1', 'tr2')},
                                    R4_midpoint={k: midpoint[k] for k in DERIVED},
                                    R4_midpoint_relative_change_percent=str(100*abs(D(midpoint['C_TRACE_COMPONENT'])/D(midpoint['C_RAW'])))))
    return dict(success=True, job_id=5942805, configuration=100, tau=4,
                piters=[4, 40], raw_records=len(rows), combined_positions_per_piter=83,
                raw_regression_maximum_difference=0.0,
                scope='Single-configuration diagnostic; no ensemble uncertainty or physical-normalization claim',
                convention=json.loads((INPUT/'SUMMARY.json').read_text())['convention'],
                launches=launches, comparisons=comparisons)


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export(output):
    rows, launches = validated_rows()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output/'profiles_raw_vs_traceless_n100_tau4.csv', rows)
    data = summary(rows, launches)
    (output/'SUMMARY.json').write_text(json.dumps(data, indent=2)+'\n')
    return rows, data


def render(rows, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.hashsalt': 'qx-clover-trace-pilot'})
    selected = {p: sorted([r for r in rows if (int(r['piter']), int(r['axis']), int(r['R']),
                                             int(r['long']), int(r['tr2'])) == (p, -1, 4, 2, 0)],
                           key=lambda r: int(r['tr1'])) for p in (4, 40)}
    x = [int(r['tr1']) for r in selected[40]]
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for field, label, color in [('C_RAW', 'Raw clover, p40', '#176b8e'),
                                ('C_TL', 'Traceless clover, p40', '#278260')]:
        axes[0].plot(x, [float(r[field]) for r in selected[40]], 'o-', color=color, label=label)
    for p, marker in [(4, 'x--'), (40, 'o-')]:
        axes[1].plot(x, [float(r['C_TRACE_COMPONENT']) for r in selected[p]], marker, label='piter=%d' % p)
    axes[0].set_ylabel(r'$C$')
    axes[1].set_ylabel(r'$C_{\rm RAW}-C_{\rm TL}$')
    axes[1].set_xlabel(r'$r_\perp/a$')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
        ax.set_xticks(x)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    fig.suptitle('n100 only; '+r'$\tau/a=4$, $R/a=4$, longitudinal midpoint'+'\nRaw and colour-traceless six-plane clover')
    fig.subplots_adjust(left=.14, right=.97, bottom=.10, top=.85, hspace=.23)
    for ext in ('png', 'svg'):
        fig.savefig(str(Path(output)/('clover_raw_vs_traceless_n100_tau4.'+ext)), dpi=180,
                    **({'metadata': {'Date': None}} if ext == 'svg' else {}))
    plt.close(fig)


def verify():
    rows, launches = validated_rows()
    if read_csv(HERE/'results/profiles_raw_vs_traceless_n100_tau4.csv') != rows:
        raise ValueError('Derived pilot CSV differs')
    if json.loads((HERE/'results/SUMMARY.json').read_text()) != summary(rows, launches):
        raise ValueError('Derived pilot summary differs')
    print('PASS pilot: 664 primaries re-parsed; raw regression exact; trace identities and derived values reproduced')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--plots', action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify()
    else:
        rows, _ = export(HERE/'results')
        if args.plots:
            render(rows, HERE/'results')
        print('PASS pilot numerical exports')


if __name__ == '__main__':
    main()
