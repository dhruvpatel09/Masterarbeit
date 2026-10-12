#!/usr/bin/env python3
"""Plot paired systematic shifts relative to preliminary control errors."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
with (HERE/'results/ensemble_13_combined.csv').open(newline='') as handle:
    rows = list(csv.DictReader(handle))

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.4), sharey=True)
colors = {1: '#2378A0', 2: '#D88922', 4: '#8146A6'}
markers = {1: 'o', 2: 's', 4: '^'}
for ax, tau in zip(axes, (2, 4, 6)):
    selected = [r for r in rows if int(r['tau']) == tau]
    for separation in (1, 2, 4):
        subset = [r for r in selected if int(r['R']) == separation]
        ax.scatter([float(r['control_signal_over_se']) for r in subset],
                   [float(r['absolute_shift_over_control_se']) for r in subset],
                   s=24, alpha=.75, color=colors[separation], marker=markers[separation], label=f'R={separation}')
    ax.set_title(rf'$\tau={tau}$  (83 probe positions)')
    ax.set_xlabel(r'$|C_{40}|\,/\,\mathrm{SE}(C_{40})$')
    ax.set_yscale('log')
    ax.grid(axis='y', alpha=.18)
axes[0].set_ylabel(r'$|C_4-C_{40}|\,/\,\mathrm{SE}(C_{40})$')
fig.suptitle('qX n110–230: projection sensitivity of combined-axis profiles', fontsize=13, y=.98)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, loc='upper center', bbox_to_anchor=(.5, .93), ncol=3)
fig.text(.5, .025, '13 configurations · midpoint 98 · paired delete-one errors · autocorrelation not assessed', ha='center', fontsize=10)
fig.tight_layout(rect=(0, .06, 1, .84))
fig.savefig(HERE/'results/projection_sensitivity.png', dpi=180)
fig.savefig(HERE/'results/projection_sensitivity.svg')
