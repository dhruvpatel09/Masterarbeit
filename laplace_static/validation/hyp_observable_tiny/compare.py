#!/usr/bin/env python3
"""Compare complete printed numerical records; no external Python packages."""
import math
from pathlib import Path
import re
import sys

number = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')

def read(path):
    lines = [s for s in path.read_text().splitlines() if not s.startswith('HYP_INTEGRATION ')]
    if not any(s.startswith('COMBINED_DIAGNOSTIC ') for s in lines):
        raise RuntimeError(f'Missing observable summary: {path}')
    if re.search(r'(?<![A-Za-z_])[-+]?(?:nan|inf)(?![A-Za-z_])', '\n'.join(lines), re.I):
        raise RuntimeError(f'Non-finite output: {path}')
    return lines

def compare(a, b, label, atol=1e-13, rtol=1e-12):
    if len(a) != len(b):
        raise RuntimeError(f'{label}: record count differs')
    worst = 0.0
    for i, (x,y) in enumerate(zip(a,b),1):
        if number.sub('#',x) != number.sub('#',y):
            raise RuntimeError(f'{label}: record layout differs at line {i}')
        for sx,sy in zip(number.findall(x),number.findall(y)):
            vx,vy=float(sx),float(sy)
            d=abs(vx-vy)
            worst=max(worst,d)
            if not math.isfinite(d) or d > atol+rtol*max(abs(vx),abs(vy)):
                raise RuntimeError(f'{label}: line {i}, {sx} vs {sy}')
    print(f'PASS {label}: max_printed_absolute_difference={worst:.16e}')

def probes(lines):
    result=[]
    for line in lines:
        if line.startswith(('PLAQUETTE ', 'CLOVER_')):
            result.append(line)
        else:
            # Compare every printed scalar probe value, including unsummed sites.
            for token in line.split():
                if token.startswith(('S_raw=', 'S_RAW_', 'S_RAW=')):
                    result.append(token)
    if len(result)<10:
        raise RuntimeError('Too few probe records')
    return result

def combined(lines):
    return next(s for s in lines if s.startswith('COMBINED_DIAGNOSTIC '))

def main():
    root=Path(sys.argv[1])
    data={name:read(root/(name+'.out')) for name in ('baseline','off','zero','hyp40','hyp80')}
    compare(data['baseline'],data['off'],'baseline vs HYP-disabled')
    compare(data['baseline'],data['zero'],'baseline vs zero-alpha')
    compare(data['hyp40'],data['hyp80'],'HYP observables: piter40 vs piter80')
    for name in ('off','zero','hyp40','hyp80'):
        compare(probes(data['baseline']),probes(data[name]),'original probe vs '+name,atol=1e-14,rtol=0)
    for name in data:
        print(name.upper(),combined(data[name]))
    print('PASS temporal-HYP observable integration on tiny configuration')

if __name__=='__main__':
    main()
