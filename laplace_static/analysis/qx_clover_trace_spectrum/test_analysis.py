#!/usr/bin/env python3
"""Check matched cohorts, nonlinear ratios, paired errors and compressed evidence."""
import gzip
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import tempfile
import unittest
from decimal import Decimal as D
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('qx_spectrum_analysis',str(HERE/'reproduce.py'))
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)


class AnalysisTests(unittest.TestCase):
    def test_manifest_inputs_survive_git_add_and_checkout(self):
        sources=json.loads((HERE/'inputs/SOURCES.json').read_text())
        relative=HERE.relative_to(a.ROOT)
        with tempfile.TemporaryDirectory() as temp:
            staged=Path(temp)/'staged';staged.mkdir()
            checked=Path(temp)/'checked';checked.mkdir()
            (staged/'.gitignore').write_bytes((a.ROOT/'.gitignore').read_bytes())
            target=staged/relative;target.mkdir(parents=True)
            for name in ('.gitignore','.gitattributes'):
                (target/name).write_bytes((HERE/name).read_bytes())
            for name in sources['files']:
                dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(str(HERE/name),str(dest))
            def git(*args):
                return subprocess.run(['git','-C',str(staged),'-c','core.excludesfile=/dev/null',*args],
                    check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
            git('init','-q')
            git('add','--','.gitignore',relative.as_posix())
            tracked=set(git('ls-files','-z').decode().split('\0'))
            missing=[name for name in sources['files'] if (relative/name).as_posix() not in tracked]
            self.assertEqual(missing,[],'Manifest inputs excluded by Git ignore rules')
            git('checkout-index','--all','--prefix='+str(checked)+'/')
            for name,expected in sources['files'].items():
                self.assertEqual(a.sha(checked/relative/name),expected,name)

    def synthetic(self):
        raw,trace={},{}
        for c,row in enumerate(((1,2,1,'.1','.3'),(2,8,1,'.2','.6'),(3,18,1,'.3','.9'),(4,32,1,'.4','1.2'))):
            values=tuple(D(str(x)) for x in row)
            for p in (4,40,80):raw[c,p]=values[:3];trace[c,p]=values
        return raw,trace

    def test_ratio_of_means(self):
        raw,trace=self.synthetic();r=a.paired_statistics(raw,trace,[0,1,2,3])
        self.assertAlmostEqual(float(r['C_RAW_p40']),60/10-1,places=14)
        self.assertNotAlmostEqual(float(r['C_RAW_p40']),(1+3+5+7)/4,places=10)

    def test_identical_controls_have_zero_paired_difference_error(self):
        raw,trace=self.synthetic();r=a.paired_statistics(raw,trace,[0,1,2,3])
        self.assertGreater(r['C_RAW_p40_jackknife_se'],0)
        for definition in ('RAW','TL'):
            self.assertEqual(r['delta_C_%s_p80_minus_p40'%definition],0)
            self.assertEqual(r['delta_C_%s_p80_minus_p40_jackknife_se'%definition],0)

    def test_jackknife_against_independent_float_calculation(self):
        raw,trace=self.synthetic()
        for c in range(4):
            v=raw[c,4];raw[c,4]=(v[0],v[1]+D(c+1)/100,v[2])
            v=trace[c,4];trace[c,4]=(v[0],v[1]+D(c+1)/100,v[2],v[3],v[4])
        result=a.paired_statistics(raw,trace,[0,1,2,3])
        samples=[]
        for omitted in range(4):
            def estimate(p):
                selected=[tuple(float(v) for v in raw[c,p]) for c in range(4) if c!=omitted]
                return sum(v[1] for v in selected)/sum(v[0] for v in selected)-sum(v[2] for v in selected)/3
            samples.append(100*(estimate(4)-estimate(40))/estimate(40))
        center=sum(samples)/4
        expected=math.sqrt(3/4*sum((v-center)**2 for v in samples))
        self.assertAlmostEqual(float(result['relative_delta_percent_RAW_p4_jackknife_se']),expected,places=11)

    def test_unmatched_configuration_is_rejected(self):
        raw,trace=self.synthetic();del trace[2,80]
        with self.assertRaises(ValueError):a.paired_statistics(raw,trace,[0,1,2,3])

    def test_missing_duplicate_and_unmeasured_traceless_prescription_rejected(self):
        row=a.read_csv(a.JOB/'profiles_raw_and_trace.csv')[0]
        for rows in ([],[row,row],[dict(row,piter='10')],[dict(row,configuration='640')]):
            with self.assertRaises(ValueError):a.group_rows(rows,a.TL_PITERS,('L_RE','LS_RE','S_RAW','T_TRACE','LT_RE'))

    def test_compressed_log_digest_detects_changed_content(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'log.gz';p.write_bytes(gzip.compress(b'original'))
            digest=hashlib.sha256(b'original').hexdigest()
            self.assertEqual(a.log_bytes(p,digest),b'original')
            p.write_bytes(gzip.compress(b'changed'))
            with self.assertRaises(ValueError):a.log_bytes(p,digest)

    def test_published_controls_reproduce_independent_subset_mean(self):
        sources,_=a.verify_inputs()
        baseline=a.selected_rows(a.ROOT/sources['raw_reference'])
        selected=[r for r in baseline if tuple(int(r[k]) for k in a.KEYS)==(4,-1,4,2,0,0) and int(r['piter'])==40]
        self.assertEqual(len(selected),4)
        from decimal import localcontext
        with localcontext() as context:
            context.prec=50
            independent=sum((D(r['LS_RE']) for r in selected),D(0))/sum((D(r['L_RE']) for r in selected),D(0))-sum((D(r['S_RAW']) for r in selected),D(0))/4
        rows=a.read_csv(HERE/'results/profiles_ensemble_4configs.csv')
        current=next(r for r in rows if tuple(int(r[k]) for k in a.KEYS)==(4,-1,4,2,0,0))
        self.assertEqual(D(current['C_RAW_p40']),independent)
        self.assertNotIn('C_TL_p10',current)
        self.assertNotIn('C_TL_p20',current)


if __name__=='__main__':unittest.main()
