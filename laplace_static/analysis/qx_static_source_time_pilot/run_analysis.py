"""Recompute all descriptive results into a separate empty destination."""
import argparse
from pathlib import Path
from analysis import compute
from common import csv_write, write_json
from replay import cold_replay


def run(output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be empty')
    output.mkdir(parents=True, exist_ok=True)
    data = cold_replay(progress=True)
    tables, summary = compute(data)
    for name, rows in tables.items():
        csv_write(output / (name + '.csv'), rows)
    write_json(output / 'SUMMARY.json', summary)
    write_json(output / 'VALIDATION.json', data['report'])
    print('PASS all descriptive source-time tables; no new measurements or campaign submission.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    run(parser.parse_args().output)
