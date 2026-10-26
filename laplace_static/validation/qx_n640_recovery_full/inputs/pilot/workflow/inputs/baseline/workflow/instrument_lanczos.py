#!/usr/bin/env python3
"""Add observation-only trace records to the captured qcd_lanczosCheby source."""
import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Instrumentation anchor is absent or ambiguous: ' + old[:80])
    return text.replace(old, new, 1)


def instrument(text):
    start = text.index('int qcd_lanczosCheby(')
    stop = text.index('int qcd_lanczosReortho(', start)
    before, body, after = text[:start], text[start:stop], text[stop:]
    body = replace_once(body,
        '   while(counter <= par->maxRestarts)\n   {',
        '   while(counter <= par->maxRestarts)\n   {\n'
        '      double qx_cycle_start = MPI_Wtime();\n'
        '      if(myid==root3d) {\n'
        '         printf("QX_CYCLE_BEGIN cycle=%d m=%d\\n", counter, m); fflush(stdout);\n'
        '      }')
    for name, call, code in (
        ('DSYEVR', 'useDSYEVR(T, evalT, evecT, m, par->redTol)', 3),
        ('DSTEMR', 'useDSTEMR(&(alpha[1]), &(beta[1]), m, evalT, evecT,par->redTol)', 4)):
        anchor = 'if(' + call + ' !=0)' if name == 'DSYEVR' else 'if( ' + call + ' != 0)'
        replacement = ('lapack_int qx_info = ' + call + ';\n'
            '         if(myid==root3d) {\n'
            '            printf("QX_REDUCED cycle=%d method=' + name + ' info=%d m=%d\\n", counter, (int)qx_info, m); fflush(stdout);\n'
            '         }\n         if(qx_info != 0)')
        body = replace_once(body, anchor, replacement)
    body = replace_once(body, '         resErr = normR/fabs(la);',
        '         resErr = normR/fabs(la);\n'
        '         if(myid==root3d)\n'
        '            printf("QX_RITZ cycle=%d index=%d lambda=%.17e residual_norm=%.17e relative_residual=%.17e\\n", counter, i, la, normR, resErr);')
    body = replace_once(body, '      counter += 1;\n      if (maxRes<par->res)',
        '      if(myid==root3d) {\n'
        '         int qx_checked = (i < 0) ? nvec : nvec-i;\n'
        '         printf("QX_CYCLE_END cycle=%d checked=%d max_checked_res=%.17e lambda_Nv=%.17e restarts=%d nRed=%d maxRed=%d reorthogonalizations=%d seconds=%.9f\\n",\n'
        '                counter, qx_checked, maxRes, eval[nvec-1], restarts, nRed, maxRed, reorthoCount, MPI_Wtime()-qx_cycle_start);\n'
        '         fflush(stdout);\n'
        '      }\n'
        '      counter += 1;\n      if (maxRes<par->res)')
    return before + body + after


def generate(output):
    report = json.loads((HERE/'INPUTS.json').read_text())
    source = HERE/'inputs/qcd_lanczos.c'
    expected = report['preserved_inputs']['inputs/qcd_lanczos.c']
    if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
        raise ValueError('Preserved solver source changed')
    Path(output).write_text(instrument(source.read_text()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output')
    generate(parser.parse_args().output)
