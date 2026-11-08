#!/usr/bin/env python3
"""Generate the qX radius and plateau catalogue from pinned published inputs."""
import argparse
import csv
import json
from pathlib import Path

import analysis
import approach_model
from data import HERE, load


def write_json(path,value):
    temporary = Path(str(path)+'.pending')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def write_csv(path,rows):
    temporary = Path(str(path)+'.pending')
    with temporary.open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    temporary.replace(path)


def generate(destination,replace=False):
    destination = Path(destination)
    if destination.exists() and any(destination.iterdir()) and not replace:
        raise ValueError('Output is not empty; use a new directory or --replace')
    destination.mkdir(parents=True,exist_ok=True)
    data = load()
    tables,summary = analysis.compute(data)
    plan = json.loads((HERE/'ANALYSIS_PLAN.json').read_text())
    rows,solutions = approach_model.compute(data,plan)
    tables['approach_fits'] = rows
    summary['table_rows']['approach_fits'] = len(rows)
    summary['approach_status_counts'] = {s:sum(r['status']==s for r in rows) for s in sorted(set(r['status'] for r in rows))}
    summary['approach_diagnostic_screen_count'] = sum(r['diagnostic_screen'] for r in rows)
    summary['approach_replay'] = 'Archived variable-projection candidates are re-evaluated; optimizer searches are not repeated by --verify.'
    for name,values in tables.items():
        write_csv(destination/(name+'.csv'),values)
    write_json(destination/'approach_solutions.json',solutions)
    write_json(destination/'SUMMARY.json',summary)
    from render_plots import render
    render(destination)
    print('Wrote {} tables and complete radius/time diagnostics to {}'.format(len(tables),destination))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'results')
    parser.add_argument('--replace',action='store_true')
    args = parser.parse_args()
    generate(args.output,args.replace)
