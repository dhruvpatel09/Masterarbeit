#!/usr/bin/env python3
"""Scientific controls for covariance propagation and matched samples."""
import unittest
from decimal import Decimal as D,localcontext
import reproduce as a


class Tests(unittest.TestCase):
    def test_correlated_cancellation(self):
        samples=[tuple(map(D,v)) for v in [(1,6,2,0,0),(2,12,2,0,0),(3,18,2,0,0),(4,24,2,0,0)]]
        point,delta=a.channel(samples,'RAW')
        self.assertEqual(point,D(4));self.assertEqual(a.se(delta),0)

    def test_sample_mean_standard_error(self):
        x=[D(v) for v in (1,2,3,4)]
        with localcontext() as ctx:
            ctx.prec=50
            self.assertEqual(a.covariance(x,x),D(5)/12)

    def test_raw_trace_identity_with_covariance(self):
        values=[tuple(map(D,v)) for v in [(2,10,3,1,2),(3,16,4,2,5),(4,20,4,1,4),(5,26,5,2,8)]]
        with localcontext() as ctx:
            ctx.prec=50
            raw,dr=a.channel(values,'RAW');trace,dt=a.channel(values,'TRACE');tl,dl=a.channel(values,'TL')
            self.assertLess(abs(tl-(raw-trace)),D('1e-48'))
            self.assertTrue(all(abs(x-y+z)<D('1e-48') for x,y,z in zip(dl,dr,dt)))
            self.assertLess(abs(a.covariance(dl,dl)-(a.covariance(dr,dr)+a.covariance(dt,dt)-2*a.covariance(dr,dt))),D('1e-48'))

    def test_covariance_is_paired(self):
        x=list(map(D,(1,2,3,4)));y=list(map(D,(4,3,2,1)))
        self.assertEqual(a.covariance(x,y),a.covariance(x,x).copy_negate())
        self.assertEqual(a.se([v+w for v,w in zip(x,y)]),0)

    def test_invalid_samples(self):
        with self.assertRaises(ValueError):a.channel([(D(0),)*5]*3,'RAW')
        with self.assertRaises(ValueError):a.covariance([D(0)],[D(0)])


if __name__=='__main__':unittest.main()
