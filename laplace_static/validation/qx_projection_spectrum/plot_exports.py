#!/usr/bin/env python3
"""Combine matched numerical exports and plot all measured projection counts."""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
import campaign as c


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extra-raw',type=Path,action='append',required=True)
    parser.add_argument('--extra-windows',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();c.verify_sources()
    command=[sys.executable,str(c.ROOT/'laplace_static/analysis/qx_hyp_projection/make_plots.py'),'--output',str(args.output)]
    for p in args.extra_raw:command+=['--extra-raw',str(p)]
    for p in args.extra_windows:command+=['--extra-windows',str(p)]
    subprocess.run(command,check=True)
    summary=json.loads((args.output/'SUMMARY.json').read_text())
    configs=summary['configurations'];piters=summary['profile_piters']
    with (args.output/f'su3_central_windows_{len(configs)}configs.csv').open(newline='') as stream:
        windows=list(csv.DictReader(stream))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'qx-full-piter-windows'})
    fig,axes=plt.subplots(1,2,figsize=(13,5.4),layout='constrained')
    colors={4:'#cd5519',10:'#9c5f9b',20:'#278260',40:'#176b8e'}
    for ax,field,title in zip(axes,('unit_max','det_max'),('Frobenius unitarity defect','Determinant defect')):
        for p in piters:
            values=[max(float(r[field]) for r in windows if int(r['configuration'])==cfg and int(r['piter'])==p) for cfg in configs]
            ax.plot(configs,values,'o-',ms=4,color=colors.get(p),label=f'piter={p}')
        ax.set_yscale('log');ax.grid(axis='y',alpha=.2);ax.set_xlabel('Configuration n')
        ax.set_ylabel('Largest maximum over the three measured windows');ax.set_title(title);ax.legend()
    fig.suptitle(f'Central transporter SU(3) accuracy: {len(configs)} matched configurations\nOne temporal HYP step; midpoint98, τ/a=2,4,6',fontsize=14)
    fig.supxlabel('Measured data only. These Frobenius diagnostics differ from the entry norm used in the n100 iteration spectrum.',fontsize=10)
    name=f'su3_central_windows_all_piters_{len(configs)}configs'
    fig.savefig(args.output/(name+'.png'),dpi=180)
    fig.savefig(args.output/(name+'.svg'),metadata={'Date':None})
    path=args.output/(name+'.svg');path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    print('All-piter profiles and central-window plots:',args.output)


if __name__=='__main__':main()
