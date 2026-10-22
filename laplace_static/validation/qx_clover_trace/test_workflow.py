#!/usr/bin/env python3
"""Check clover algebra, covariance bookkeeping, parsers and evidence export."""
import copy
import ctypes
import importlib.util
import math
import random
import shutil
import subprocess
import tempfile
import unittest
import zipfile
from decimal import Decimal, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location('trace_test_'+name, str(HERE/(name+'.py')))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


workflow = load('campaign')
collector = load('collect_n640')


class CloverAlgebra(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('cc') or shutil.which('gcc')
        if not compiler:
            raise unittest.SkipTest('C compiler unavailable')
        cls.temp = tempfile.TemporaryDirectory()
        path = Path(cls.temp.name)
        (path/'test.c').write_text('''#include "clover_trace.h"
void evaluate(const double *r,const double *i,double *out) {
    qx_clover_components c=qx_clover_plane((const double (*)[3])r,(const double (*)[3])i);
    out[0]=c.raw;out[1]=c.trace_component;out[2]=c.traceless_direct;
}
''')
        subprocess.run([compiler, '-std=c11', '-shared', '-fPIC', '-O2', '-Wall', '-Wextra', '-Werror',
                        '-I'+str(HERE), str(path/'test.c'), '-o', str(path/'test.so')], check=True)
        cls.library = ctypes.CDLL(str(path/'test.so'))
        cls.evaluate = cls.library.evaluate
        cls.evaluate.argtypes = [ctypes.POINTER(ctypes.c_double)]*3
        cls.evaluate.restype = None

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'temp'):
            cls.temp.cleanup()

    def components(self, q):
        r = (ctypes.c_double*9)(*[q[a][b].real for a in range(3) for b in range(3)])
        i = (ctypes.c_double*9)(*[q[a][b].imag for a in range(3) for b in range(3)])
        out = (ctypes.c_double*3)()
        self.evaluate(r, i, out)
        return list(out)

    def test_zero_pure_trace_and_traceless(self):
        for diag in ([0j]*3, [4j]*3, [4j, -4j, 0j]):
            q = [[diag[a] if a == b else 0j for b in range(3)] for a in range(3)]
            raw, trace, tl = self.components(q)
            self.assertAlmostEqual(raw, trace+tl)
            if diag == [4j]*3:
                self.assertAlmostEqual(tl, 0)
                self.assertAlmostEqual(raw, 3)
            else:
                self.assertAlmostEqual(trace, 0)

    def test_random_matrices_against_independent_complex_algebra(self):
        rng = random.Random(640)
        for scale in (1e-9, 1e-3, 1, 1e3):
            for _ in range(64):
                q = [[scale*complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
                      for b in range(3)] for a in range(3)]
                f = [[(q[a][b]-q[b][a].conjugate())/8 for b in range(3)] for a in range(3)]
                trace = sum(f[a][a] for a in range(3))
                expected = [sum(abs(f[a][b])**2 for a in range(3) for b in range(3)),
                            abs(trace)**2/3,
                            sum(abs(f[a][b]-(trace/3 if a == b else 0))**2
                                for a in range(3) for b in range(3))]
                actual = self.components(q)
                for a, b in zip(actual, expected):
                    self.assertTrue(math.isclose(a, b, rel_tol=3e-14, abs_tol=1e-28))
                self.assertTrue(math.isclose(actual[0]-actual[1], actual[2], rel_tol=3e-14, abs_tol=1e-28))

    def test_covariance_under_unitary_colour_rotation(self):
        q = [[complex(a+b, a-b+.7) for b in range(3)] for a in range(3)]
        c, s = math.cos(.7), math.sin(.7)
        u = [[c, 1j*s, 0], [1j*s, c, 0], [0, 0, 1]]
        rotated = [[sum(u[a][i]*q[i][j]*complex(u[b][j]).conjugate()
                        for i in range(3) for j in range(3)) for b in range(3)] for a in range(3)]
        for a, b in zip(self.components(q), self.components(rotated)):
            self.assertAlmostEqual(a, b, places=13)


class WorkflowTests(unittest.TestCase):
    def test_source_pins_and_pilot_plan(self):
        workflow.verify_sources()
        self.assertEqual(workflow.case_plan(), [dict(configuration=100, tau=4, piter=p) for p in (4, 40)])

    def test_connected_trace_identity_and_ratio_of_means(self):
        row = dict(L_RE='2', LS_RE='10', S_RAW='3', T_TRACE='.5', LT_RE='2')
        got = workflow.connected_from_primaries(row)
        self.assertEqual(got, dict(C_RAW=Decimal(2), C_TRACE_COMPONENT=Decimal('.5'), C_TL=Decimal('1.5')))
        samples = [dict(L_RE='1', LS_RE='3', S_RAW='1', T_TRACE='.2', LT_RE='.5'),
                   dict(L_RE='3', LS_RE='17', S_RAW='5', T_TRACE='.8', LT_RE='3.5')]
        with localcontext() as context:
            context.prec = 50
            means = {k: str(sum(Decimal(r[k]) for r in samples)/2) for k in row}
            mean_ratios = sum(workflow.connected_from_primaries(r)['C_TL'] for r in samples)/2
        self.assertEqual(workflow.connected_from_primaries(means), got)
        self.assertNotEqual(mean_ratios, got['C_TL'])

    def fixture(self):
        raw, lines = {}, ['CLOVER_TRACE_CHECK count=262144 identity_max=1e-15 status=PASS']
        for kind in ('local', 'cross_boundary', 'spatial_average'):
            lines += ['PROBE kind=%s S_RAW=2' % kind,
                      'CLOVER_TRACE_PROBE kind=%s t=98 T_TRACE=.3 S_TL=1.7' % kind]
        for axis in (0, 1, 2, -1):
            for r in (1, 2, 4):
                offsets = {(0, 1, 2)} | {(l, t, 0) for l in range(-1, r+2) for t in range(-2, 3)}
                for l, t, v in sorted(offsets):
                    count = 64**3*(3 if axis == -1 else 1)
                    raw[(axis, r, l, t, v)] = dict(Nsrc=str(count), S_RAW='2')
                    prefix = 'TRACE_COMBINED' if axis == -1 else 'TRACE_AXIS axis=%d' % axis
                    lt = .4 if axis == -1 else .2*(axis+1)
                    lines.append('%s R=%d long=%d tr1=%d tr2=%d Nsrc=%d T_TRACE=.3 LT_RE=%.17e LT_IM=0' %
                                 (prefix, r, l, t, v, count, lt))
        return raw, '\n'.join(lines)

    def test_trace_parser_matches_all_axes_and_rejects_missing_rows(self):
        raw, text = self.fixture()
        check = workflow.verify_sources().checker()
        result = workflow.trace_records(text, raw, check.avg)
        self.assertEqual(len(result), 332)
        with self.assertRaises(ValueError):
            workflow.trace_records('\n'.join(text.splitlines()[:-1]), raw, check.avg)
        with self.assertRaises(RuntimeError):
            workflow.trace_records(text.replace('LT_RE=4.00000000000000022e-01', 'LT_RE=8e-1'), raw, check.avg)

    def test_names_and_nonfinite_primaries_rejected(self):
        with self.assertRaises(ValueError):
            workflow.folder('../escape')
        with self.assertRaises(ValueError):
            workflow.connected_from_primaries(dict(L_RE='0', LS_RE='1', S_RAW='1', T_TRACE='0', LT_RE='0'))
        with self.assertRaises(ValueError):
            workflow.connected_from_primaries(dict(L_RE='1', LS_RE='nan', S_RAW='1', T_TRACE='0', LT_RE='0'))

    def test_n640_collector_exports_logs_and_skips_binary_and_symlinks(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            log = root/'mental/validation/qx_retry_n640_5925063/log'
            log.mkdir(parents=True)
            (log/'mental.log').write_text('Error! Lanczos on time-slice 1 did not converge!\n')
            (log/'mental.in').write_text('maxRestarts 60\n')
            (log/'evec.h5').write_bytes(b'binary')
            (log/'binary.out').write_bytes(b'\x00')
            outside = root/'private.txt'
            outside.write_text('do not include')
            (log/'external.out').symlink_to(outside)
            output = root/'review.zip'
            manifest = collector.export(root, output)
            self.assertEqual(len(manifest['files']), 2)
            self.assertEqual(len(manifest['skipped']), 1)
            with zipfile.ZipFile(str(output)) as archive:
                self.assertEqual(len(archive.namelist()), 3)
                self.assertIsNone(archive.testzip())
            with self.assertRaises(FileExistsError):
                collector.export(root, output)

    def test_empty_collector_writes_no_archive(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'absent.zip'
            with self.assertRaises(ValueError):
                collector.export(Path(temp), output)
            self.assertFalse(output.exists())

    def test_n640_logs_outside_repository_are_exported(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            root, runs = base/'repository', base/'production/mental_runs'
            root.mkdir()
            run = runs/'n640_evec_Nv100_mr60/log'
            run.mkdir(parents=True)
            (run/'mental.log').write_text('Lanczos on time-slice 1 did not converge\n')
            (run/'eigenvectors').mkdir()
            (run/'eigenvectors'/'n640.out').write_text('never archive mode data')
            self.assertEqual(collector.candidates(root, [run]), [run/'mental.log'])
            output = base/'review.zip'
            manifest = collector.export(root, output, [runs])
            self.assertEqual(len(manifest['files']), 1)
            entry = manifest['files'][0]
            self.assertTrue(entry['archive_path'].startswith('external/'))
            self.assertEqual(entry['path'], str(run/'mental.log'))
            with zipfile.ZipFile(str(output)) as archive:
                self.assertEqual(archive.read(entry['archive_path']),
                                 b'Lanczos on time-slice 1 did not converge\n')

    def test_collector_finds_job_reference_in_generically_named_controller(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            control = root/'mental/production/controller'
            control.mkdir(parents=True)
            (control/'SUBMISSION.json').write_text('{"job_id": "5925063"}\n')
            (control/'unrelated.json').write_text('{"configuration": 650}\n')
            self.assertEqual(collector.candidates(root), [control/'SUBMISSION.json'])


if __name__ == '__main__':
    unittest.main()
