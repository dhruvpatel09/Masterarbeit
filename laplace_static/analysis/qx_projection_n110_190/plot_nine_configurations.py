"""Create standalone scientific figures from the paired nine-config analysis."""
import csv
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(directory):
    directory = Path(directory)
    with (directory/'ensemble_profiles_combined_axes.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(3,3,figsize=(12,10),layout='constrained')
    delta,delta_axes = plt.subplots(3,3,figsize=(12,10),layout='constrained')
    for i,tau in enumerate((2,4,6)):
        for j,R in enumerate((1,2,4)):
            selected = sorted([r for r in rows if int(r['tau'])==tau and int(r['R'])==R
                               and int(r['tr1'])==0 and int(r['tr2'])==0],key=lambda r:int(r['long']))
            x = [int(r['long'])-R/2 for r in selected]
            ax = axes[i,j]
            ax.errorbar(x,[1e3*float(r['C_p40']) for r in selected],
                        yerr=[1e3*float(r['C_p40_unblocked_jk_se']) for r in selected],
                        fmt='o-',color='#1d5d87',capsize=3,label='piter 40 + pilot jackknife')
            ax.plot(x,[1e3*float(r['C_p4']) for r in selected],'x--',color='#cb632d',
                    markersize=7,label='piter 4')
            ax.set_title(f'τ={tau}, R={R}')
            ax.axvline(-R/2,color='0.7',lw=.8,ls=':')
            ax.axvline(R/2,color='0.7',lw=.8,ls=':')
            ax.grid(alpha=.2)
            dx = delta_axes[i,j]
            dx.errorbar(x,[1e9*float(r['delta_p4_minus_p40']) for r in selected],
                        yerr=[1e9*float(r['paired_delta_unblocked_jk_se']) for r in selected],
                        fmt='o-',color='#6f3d82',capsize=3)
            dx.axhline(0,color='0.4',lw=.8)
            dx.set_title(f'τ={tau}, R={R}')
            dx.grid(alpha=.2)
            if i==2:
                ax.set_xlabel('Longitudinal offset from midpoint [lattice sites]')
                dx.set_xlabel('Longitudinal offset from midpoint [lattice sites]')
            if j==0:
                ax.set_ylabel('Connected C × 10³ [raw lattice units]')
                dx.set_ylabel('(C₄ − C₄₀) × 10⁹ [raw lattice units]')
    axes[0,0].legend(fontsize=8,loc='best')
    fig.suptitle('qX n110–190: longitudinal profiles at midpoint 98\nNine configurations; axes combined before ratio; unblocked jackknife is preliminary',fontsize=14)
    delta.suptitle('Paired projection sensitivity on the same configurations\nDelete the same configuration in both prescriptions; no autocorrelation correction',fontsize=14)
    fig.savefig(directory/'longitudinal_profiles.png',dpi=180)
    delta.savefig(directory/'paired_projection_difference.png',dpi=180)
    plt.close(fig);plt.close(delta)
    with (directory/'su3_windows_n110_190.csv').open(newline='') as stream:
        windows = [r for r in csv.DictReader(stream) if int(r['piter'])==4]
    fig,axes = plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for tau,color in zip((2,4,6),('#1d5d87','#46845b','#cb632d')):
        selected=sorted([r for r in windows if int(r['tau'])==tau],key=lambda r:int(r['configuration']))
        for ax,field,label in zip(axes,('unit_max','det_max'),('max ‖U†U − I‖F','max |det U − 1|')):
            ax.semilogy([int(r['configuration']) for r in selected],[float(r[field]) for r in selected],
                        'o-',color=color,label=f'τ={tau}')
            ax.set_xlabel('Configuration number');ax.set_ylabel(label)
            ax.grid(alpha=.25,which='both');ax.legend()
    fig.suptitle('piter 4 SU(3) defects on the central temporal paths\nOne temporal HYP step, α=(1,1,0.5); piter 40 controls below 8×10⁻¹⁶',fontsize=14)
    fig.savefig(directory/'su3_central_window.png',dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main(sys.argv[1])
