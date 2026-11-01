#!/usr/bin/env python3
"""Independent covariance and cutoff-kernel checks for pilot diagnostics."""
import ast
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import reproduce as r


class Tests(unittest.TestCase):
    def test_paired_normalization_cancels_common_fluctuations(self):
        xs=[1.,2.,4.]
        value,error=r.ratio(xs,xs)
        self.assertEqual(value,1.);self.assertEqual(error,0.)
        value,error=r.ratio([2*x for x in xs],xs)
        self.assertEqual(value,2.);self.assertEqual(error,0.)
        self.assertGreater(r.se(xs),0)

    def test_nonpositive_means_and_invalid_jackknife_are_retained(self):
        self.assertIsNone(r.log_ratio([-2.,-1.,1.],[1.,1.,1.])[0])
        value,error,jk,count,influence=r.log_ratio([-.9,.1,1.1],[1.,1.,1.])
        self.assertAlmostEqual(value,math.log(.1));self.assertGreater(error,0)
        self.assertIsNone(jk);self.assertEqual(count,2);self.assertEqual(len(influence),3)

    def test_log_covariance_cancels_and_python36_syntax(self):
        value,error,jk,count,influence=r.log_ratio([2.,4.,8.],[1.,2.,4.])
        self.assertAlmostEqual(value,math.log(2.));self.assertEqual(error,0.)
        self.assertEqual(jk,0.);self.assertEqual(count,3)
        for path in r.HERE.glob('*.py'):
            options=dict(feature_version=(3,6)) if sys.version_info>=(3,8) else {}
            ast.parse(path.read_text(),filename=str(path),**options)

    @unittest.skipUnless(shutil.which('gcc'),'gcc is required for the native kernel replay')
    def test_captured_cutoff_kernels_and_complete_periodic_shifts(self):
        with tempfile.TemporaryDirectory() as tmp:
            exe=Path(tmp)/'kernels'
            subprocess.run(['gcc','-std=c11','-O2','-I'+str(r.SOURCE/'tests'),
                '-I'+str(r.SOURCE/'build_provenance'),str(r.SOURCE/'tests/kernels_and_shifts.c'),
                '-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)


if __name__=='__main__':unittest.main()
