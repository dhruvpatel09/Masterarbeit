"""Fixed pilot scope and deterministic exports. No ensemble inference."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONFIGS = (170, 640, 650)
RADII = (4, 12, 24, 32)
CUTOFFS = (40, 60, 80, 100)
POSITIONS = (58, 90, 122)
BASES = ('pair', 'triple', 'full')
ANCHORS = (4, 6, 8)
LABELS = {58: 'shift58_p10', 90: 'reused_t90', 122: 'shift122_p10'}


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def csv_read(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def csv_write(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def frobenius(matrix):
    return math.sqrt(math.fsum(abs(z) ** 2 for row in matrix for z in row))


def normalized_matrix(matrix, diagonal):
    if len(diagonal) != 4 or any(x <= 0 for x in diagonal):
        raise ValueError('Invalid fixed equal-time diagonal')
    return [[matrix[i][j] / math.sqrt(diagonal[i] * diagonal[j])
             for j in range(4)] for i in range(4)]


def delta_matrix(a, b):
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def phase(value):
    return math.atan2(value.imag, value.real) if value != 0 else None


def ratio(numerator, denominator):
    return numerator / denominator if denominator != 0 else None
