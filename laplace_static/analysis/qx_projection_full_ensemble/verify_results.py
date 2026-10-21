#!/usr/bin/env python3
"""Recompute and compare every preserved full70 numerical export."""
import tempfile
from pathlib import Path
import reproduce as study


def main():
    study.verify_inputs()
    expected = sorted(p.name for p in study.RESULTS.iterdir()
                      if p.suffix in ('.csv', '.json'))
    if len(expected) != 8:
        raise ValueError('Expected eight preserved numerical exports')
    with tempfile.TemporaryDirectory() as name:
        output = Path(name)
        study.reproduce(output, make_plots=False)
        actual = sorted(p.name for p in output.iterdir())
        if actual != expected:
            raise ValueError('Numerical output inventory changed')
        for filename in expected:
            if (output/filename).read_bytes() != (study.RESULTS/filename).read_bytes():
                raise ValueError('Numerical output differs: '+filename)
    for name in study.analysis.PLOTS+study.EXTRA_PLOTS:
        for extension in ('.png', '.svg'):
            path = study.RESULTS/(name+extension)
            if not path.is_file() or not path.stat().st_size:
                raise ValueError('Missing figure: '+path.name)
    print('PASS: 70 configurations; 210 cases; four piters; 278880 raw rows; '
          '840 window rows; all estimates, paired errors and cohort diagnostics reproduced exactly.')


if __name__ == '__main__':
    main()
