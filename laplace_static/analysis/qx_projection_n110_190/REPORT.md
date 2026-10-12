# Nine-configuration temporal-HYP projection pilot

All 27 configuration/temporal-separation cases passed revalidation: n110–190 in steps of 10, tau=2,4,6, midpoint 98. The new 12 cases completed at 01:10:50 CEST on 10 October 2026. Both n100 tau4 regressions reproduce the earlier measurements exactly. n100 is a regression reference and is excluded from the nine-configuration averages.

The two campaigns have identical numerical C, SU(3) helper, launcher and scientific checkers; six shared runtime-dependency hashes match. All 58 actual measurement launches returned zero. Locally available manifest entries and all copied run provenance were verified. Actual archived stdout was re-parsed, and every per-case and collected CSV matched byte-for-byte. Runtime binary, gauge and HDF5 bytes are outside these review archives.

## Estimator and sampling

`C = mean(LS_RE)/mean(L_RE) - mean(S_RAW)`, with equal configuration weights and the same configurations paired for piter4 and piter40. Combined-axis raw measurements are combined before the ratio. Spatial sites, axes and tau values are not additional ensemble samples. There are 17,928 raw rows, 8,964 per-configuration paired rows, 996 averaged diagnostic rows, and 249 combined-axis averaged profile points.

The mean and delete-configuration estimates use 50-digit Decimal arithmetic to protect the vacuum subtraction from cancellation. Error columns are paired configuration jackknife diagnostics. A matched block-three deletion check is also supplied (only three blocks). Nine ordered configurations cannot establish autocorrelation, final ensemble errors, or a physical temporal plateau. No UWerr/autocorrelation fit is claimed.

## Projection sensitivity of the averaged combined-axis profiles

| tau | maximum absolute C4−C40 | maximum relative difference | maximum shift / preliminary p40 jackknife error |
|---|---:|---:|---:|
| 2 | 2.860633e-09 | 0.000791% | 9.372697e-05 |
| 4 | 7.095291e-09 | 0.002229% | 1.358496e-04 |
| 6 | 1.172334e-08 | 0.004347% | 2.981235e-04 |

The largest combined-axis averaged relative difference is 0.0043465%; the largest absolute difference is 1.172334e-8 in raw lattice units. Individual-axis averaged records reach 0.0171007% and 3.088291e-8. Individual-configuration combined-axis ratios reach 0.321693%; these are different estimators, and ratios close to zero can exaggerate relative changes. All relevant absolute values and point locations are supplied in SUMMARY.json and the CSVs.

The largest combined-axis shift is less than 0.00030 of the unblocked p40 jackknife error at the same point. The exploratory block-three counterpart is below 0.00163. These are pilot sampling comparisons, not calibrated significance statements. In the latest four-configuration campaign, median full measurement times were 89.58 s for piter4 and 334.47 s for piter40 (3.73x ratio), including gauge/eigenvector reading and observable work; this is not a projection-kernel-only benchmark.

## SU(3) defects and decision scope

Worst central-path piter4 defects occur on n170, tau6, at t=95: max Frobenius unitarity defect 6.240982e-3 and max determinant defect 3.164495e-3. Over that six-link temporal window, the RMS values are 4.976312e-6 and 2.523243e-6. No inspected central-window link has a defect above 1e-2. That 1e-2 counter is a diagnostic threshold, not an agreed acceptance criterion. All piter40 controls remain below 8e-16.

These measurements support piter4 as a documented working prescription for this central-window observable pilot because the averaged projection sensitivity is small. They do not establish strict SU(3) membership, an ensemble-wide defect tolerance, or suitability near open boundaries. Earlier n100 boundary tests showed much larger defects. Retain piter40 as the precision control, preserve the paired results, and report the n170 tail when discussing the prescription with the advisors. A requirement for link-level accuracy near machine precision would select piter40.

## Operator and conventions

Global lattice 192x64^3; observable grid 1x4x4x4, 64 MPI ranks on eight nodes, one thread/rank, solver_mpiInit / MPI_THREAD_FUNNELED. One temporal HYP step, alpha=(1,1,0.5). Spatial APE20 alpha=0.5 belongs to the saved Laplace eigensystem. The clover probe is built from the original field: F=(Q−Q†)/8, sum over six planes Tr(F†F), without color-trace subtraction or an extra factor 1/2. Results remain in raw lattice units.

## Next work

The accompanying batch prepares n200,210,220,230 using the proven four independent eight-node MPI worlds in a 32-node general allocation. Budget approximately 900.42 GiB of new vectors and 20–30 minutes after allocation, subject to actual remaining scratch quota and I/O/convergence. Follow generation with the 13-timeslice saved-mode audit before extending profile statistics. Temporal HYP piter4/40 has no effect on eigenvector generation. Continue to retain physical core-hour and Slurm unit-hour capacity accounting separately.
