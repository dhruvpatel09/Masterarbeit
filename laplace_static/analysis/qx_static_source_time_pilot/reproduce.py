"""Rebuild native replay and analysis; check labels, masks and numerical tables."""
import argparse
import math
from analysis import compute
from common import HERE, csv_read, read
from replay import cold_replay, close_summary
from verify_files import verify as verify_files


def compare_table(path, expected):
    actual = csv_read(path)
    if len(actual) != len(expected):
        raise ValueError('Table length: ' + path.name)
    for i, (a, b) in enumerate(zip(actual, expected)):
        if set(a) != set(b):
            raise ValueError('Table columns: ' + path.name)
        for key, value in b.items():
            if value is None:
                valid = a[key] == ''
            elif isinstance(value, float):
                number = float(a[key])
                valid = math.isfinite(number) and math.isfinite(value) and math.isclose(number, value, rel_tol=1e-10, abs_tol=0.)
            else:
                valid = a[key] == str(value)
            if not valid:
                raise ValueError('Numerical/metadata mismatch %s row %d field %s' % (path.name, i, key))


def verify():
    verify_files()
    data = cold_replay(progress=True)
    tables, summary = compute(data)
    for name, rows in tables.items():
        compare_table(HERE / 'results' / (name + '.csv'), rows)
    close_summary(read(HERE / 'results/SUMMARY.json'), summary)
    close_summary(read(HERE / 'results/VALIDATION.json'), data['report'])
    print('PASS source-time pilot replay: 3 native audits; 8 scans; 67,584 matrix rows; 147,456 rank blocks; 13,068 frozen projections.')
    print('PASS 68-configuration source-training replay, exact primary baseline and all descriptive position/control tables.')
    print('Three technical cases only; no production measurement or submission command.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', required=True, action='store_true')
    parser.parse_args()
    verify()
