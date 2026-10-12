"""Small independent checks of the nonlinear estimator and paired error."""
from decimal import Decimal, localcontext
import math
from analyze_nine_configurations import connected, jackknife


def row(l,ls,s=0):
    return {'L_RE':str(l),'LS_RE':str(ls),'S_RAW':str(s)}


with localcontext() as context:
    context.prec = 50
    sample = [row(1,1),row(3,9)]
    assert connected(sample) == Decimal('2.5')  # mean of ratios would give 2
    # For constant denominator, jackknife SE is the ordinary SE of the mean.
    p40 = [row(1,v) for v in range(1,10)]
    p4 = [row(1,v+Decimal('0.001')) for v in range(1,10)]
    e4,e40,edelta = jackknife({4:p4,40:p40})
    assert math.isclose(e40,math.sqrt(7.5/9),rel_tol=1e-14)
    assert e4 == e40 and edelta < 1e-45  # pairing cancels common fluctuations
    assert connected(p4)-connected(p40) == Decimal('0.001')
    # A varying paired difference must not cancel from the uncertainty.
    independent = [row(1,Decimal(v)+Decimal(v)/100) for v in range(1,10)]
    assert math.isclose(jackknife({4:independent,40:p40})[2],e40/100,rel_tol=1e-14)
    eblock = jackknife({4:p4,40:p40},3)
    assert eblock[0] == eblock[1] and eblock[2] < 1e-45
print('PASS ratio-of-means, analytical mean error, paired covariance cancellation, varying shift and matched block deletion')
