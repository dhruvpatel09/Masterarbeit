"""Estimator identities and input integrity checks; no MPI or HDF5 work."""
import csv
import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('profile_analysis', HERE/'analyze.py')
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class Tests(unittest.TestCase):
    def test_ratio_of_means_is_not_mean_of_ratios(self):
        value, means, _ = analysis.estimator([(D(1), D(2), D(0)), (D(3), D(12), D(0))])
        self.assertEqual(value, D('3.5'))
        self.assertNotEqual(value, (D(2)/D(1)+D(12)/D(3))/2)

    def test_two_value_jackknife_has_analytical_error(self):
        self.assertEqual(analysis.jk_error([D(2), D(4)]), D(1))

    def test_paired_jackknife_preserves_shared_noise(self):
        with localcontext() as context:
            context.prec = 50
            rows = {}
            for cfg, loop, ls, probe in [(110, 1, 2, '.1'), (120, 2, 6, '.3'), (130, 4, 20, '.2')]:
                rows[cfg, 40] = (D(loop), D(ls), D(probe))
                rows[cfg, 4] = (D(loop), D(ls)+D('.005')*D(loop), D(probe))
            result = analysis.paired_estimate(rows, (110, 120, 130))
            self.assertLess(abs(result['delta_C_p4_minus_p40']-D('.005')), D('1e-40'))
            self.assertLess(result['delta_C_paired_jackknife_se'], D('1e-40'))
            self.assertGreater(result['C_p40_jackknife_se'], D(0))

    def change_table(self, base, action, update_digest=True):
        meta = json.loads((base/'SOURCES.json').read_text())
        source = meta['profile_sources'][0]
        path = base/source['raw_table']
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle)
            header, rows = reader.fieldnames, list(reader)
        action(rows)
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            writer.writerows(rows)
        if update_digest:
            source['tables'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            (base/'SOURCES.json').write_text(json.dumps(meta))

    def check_bad_input(self, action, message, update_digest=True):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)/'inputs'
            shutil.copytree(HERE/'inputs', base)
            self.change_table(base, action, update_digest)
            with self.assertRaisesRegex(ValueError, message):
                analysis.read_inputs(base)

    def test_reject_duplicate(self):
        self.check_bad_input(lambda rows: rows.append(rows[0].copy()), 'Duplicate')

    def test_reject_missing_pair(self):
        self.check_bad_input(lambda rows: rows.pop(0), 'Missing paired')

    def test_reject_changed_clover(self):
        def edit(rows):
            rows[0]['S_RAW'] = '9.0'
        self.check_bad_input(edit, 'clover differs')

    def test_reject_changed_raw_bytes(self):
        self.check_bad_input(lambda rows: rows.pop(0), 'digest changed', False)


if __name__ == '__main__':
    unittest.main()
