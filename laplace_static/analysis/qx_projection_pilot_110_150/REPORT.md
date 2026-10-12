# qX central-window projection pilot: n110–n150

All 15 cases completed and passed structural/profile checks, unchanged original-field clover checks, and the piter40 SU(3) control. The n100 tau4 regression also reproduced its previous outputs exactly. These checks report piter4 defects without imposing a physics acceptance threshold.

## Setup and conventions

Five configurations (110,120,130,140,150), tau=2,4,6 with midpoint98, separations R=1,2,4, Nv=100, rho_i=1. One temporal HYP step with alpha=(1,1,0.5), compared at piter=4 and40. Projection iterations are separate from HYP smearing steps and from the spatial APE20,alpha0.5 Laplace-eigenvector operator.

The clover probe uses the original gauge field: F=(Q-Q†)/8 and S_RAW=sum over the six mu<nu planes of Tr(F†F), with no color-trace subtraction and no extra1/2. All reported connected quantities are dimensionless diagnostics in this raw convention, without a physical normalization or continuum interpretation.

Physical temporal links in the actual transport windows t0..t1-1 are audited. The open-boundary regions are outside this central-window claim. Observable jobs use grid1x4x4x4,64MPI ranks; this is distinct from the eigenvector-production grid6x4x4x4. FUNNELED initialization is retained.

## Projection diagnostics

Maxima over the three tau values for each configuration:

| Configuration | p4 max Frobenius unitarity defect | p4 max determinant defect | Max per-configuration combined-profile relative change |
|---|---:|---:|---:|
| 110 | 2.774e-06 | 1.389e-06 | 0.010585% |
| 120 | 6.681e-04 | 3.356e-04 | 0.268692% |
| 130 | 1.307e-05 | 6.552e-06 | 0.321693% |
| 140 | 2.563e-05 | 1.282e-05 | 0.013145% |
| 150 | 5.913e-05 | 2.964e-05 | 0.084902% |

The worst p4 window RMS unitarity defect is9.227e-7. No audited central-window links exceed1e-2 in either defect. The p40 control's largest unitarity/determinant defect is7.792e-16. A 1e-2 counter is a diagnostic threshold, not an established physics acceptance criterion.

The largest individual-axis relative connected change is5.0795%, at a small reference connected value (5.5601e-7); its absolute change is2.8243e-8. Small denominators can inflate these percentages. The largest per-configuration absolute connected change across all axes is8.3494e-8.

## Five-configuration raw-average comparison

For each matched profile point and prescription, compute

    C = mean_configuration(LS_RE)/mean_configuration(L_RE) - mean_configuration(S_RAW).

Each configuration has equal weight; spatial sites, axes and tau values are not extra independent configurations. This calculation does not average the per-configuration ratios.

The largest combined-axis absolute change after raw averaging is **1.675122e-08**, and the largest combined-axis relative change is **0.006449%**. Both maxima occur at tau6,R4,long5,tr1=2,tr2=0. There are996 matched mean-profile points, including249 combined-axis points. Largest all-axis absolute change is3.0020e-8; largest all-axis relative change is0.139968%.

`ensemble_profiles.csv` retains the raw means, both connected diagnostics, their paired difference, and exploratory delete-one-configuration jackknife scales. The same configuration is omitted from both prescriptions for each paired replicate. These five successive configurations have not been assessed for autocorrelation; the jackknife columns are not final errors or evidence of ensemble-wide significance. The plotted example cut is illustrative; the maxima above cover all exported points.

## Working decision and next work

The small measured profile sensitivity supports using piter4 as a working central-window pilot prescription, with piter40 retained as a control. This does not establish exact SU(3) projection or ensemble-wide acceptance. Recheck defects and paired observable sensitivity when adding configurations or changing windows; previous large open-boundary defects are not covered by this result.

Proceed with the next four eigenvector configurations n160,n170,n180,n190 using the already successful packed32-node launcher. Temporal HYP piter does not enter that eigenvector operator. After completion, audit saved residuals/normalization/Gram matrices on selected slices and extend the paired profile pilot. Existing saved-mode quality checks cover nine slices of n110–n150; tau2/6 endpoints are outside the independently audited slice set.

## Evidence and reproducibility

Source: `hyp_projection_configs_review_mixed_5918419.zip`; SHA256 `825765b3b5e25e7143feaa84126f408661076068757702776c2b938e273bade8`. The uploaded archive was CRC checked,26 archived source/checker hashes were verified, and the original15-case checkers were re-run on copied stdout. Re-exported individual raw/paired/SU(3) CSVs matched their uploaded counterparts byte for byte. This replay verifies archived numerical checks; it does not re-run MPI gauge computations or independently validate the absent libraries, gauge fields or eigenvectors.

`analyze.py` separately validates coverage, paired export arithmetic, unchanged S_RAW, and raw-before-ratio ensemble arithmetic using50-digit Decimal calculations checked against an independent double-precision mean. Original inputs are never modified. Reproduce into a new directory:

    python3 analyze.py hyp_projection_configs_review_mixed_5918419.zip --output reproduced --plot

The optional plot needs matplotlib. Other analysis uses Python's standard library. No HDF5/eigenvector data are included in this package.
