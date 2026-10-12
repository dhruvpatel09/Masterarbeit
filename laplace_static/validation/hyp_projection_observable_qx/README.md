# Four-iteration projection: observable comparison on qX_p2n100

Uses the unchanged, validated hyp_profile_scan_qx executable and the existing
n100 eigenvectors. HYP_PITER is already a supported runtime setting. No rebuild,
source patch, new eigenvectors, or change to the clover convention is required.

Measurements: one temporal HYP step, alpha=(1,1,0.5), piter=4;
midpoint t=98, tau=2,4,6, R=1,2,4, all spatial sources and three axes.
Probe offsets and record-count validation are unchanged from the earlier scan.
The original gauge field supplies the raw clover; only transport uses HYP.

The saved piter40 references are job 5904081. A fresh piter40 tau4 control must
reproduce that reference. Original clover probes and all S_RAW averages must
agree between piter4 and 40. L and LS differences are measurements, and are
never forced to agree.

Run from the Masterarbeit root:

    sha256sum -c laplace_static/validation/hyp_projection_observable_qx/package.sha256
    sbatch laplace_static/validation/hyp_projection_observable_qx/submit.sh

The allocation matches the validated scan: 8 nodes, 64 ranks, grid 1x4x4x4,
one thread per rank; test/test, 30 minutes. This differs from the eigenvector
benchmark, where the process grid is split in time as well.

Results are isolated under results_JOBID. Exports contain 1992 raw rows and
996 piter4-minus40 comparison rows. axis=-1 means the combined three-axis
average. Fractional changes divide by abs(piter40) and are blank for zero
references; values near zero can make fractional changes large. These are
single-configuration diagnostics, with no statistical errors or ensemble-level
connected-profile interpretation.

Job 5914979 established structural checks and repeatability of piter4 but found
max ||U^dagger U-I||_F=1.3828922260540608 and max |det(U)-1|=1.7457549354048014.
15678/50069504 physical temporal links exceeded a Frobenius defect of 0.01.
The paper arXiv:1012.3037 section 2.1 specifies a finite four-update projection;
the measurements here quantify its effect in this implementation/configuration.
PASS does not mean piter4 links are accurate SU(3) matrices, that all HYP staples
were independently verified, or that production accuracy has been established.

Local verification: shell/Python syntax; changed L/LS accepted as measurements;
changed clover or piter40 controls rejected; incorrect projection count rejected;
CSV schema, record counts, zero/negative reference handling and no-overwrite
behavior checked. The three structural validators are unchanged from the earlier
validated temporal-scan package. Cluster execution of this comparison is pending.
