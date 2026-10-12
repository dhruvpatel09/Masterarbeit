# piter=4 versus 40: complete n100 observable scan

Data: observable job 5915129 and full-field projection job 5914979.
Configuration qX_p2n100; Nv=100; midpoint t=98; tau=2,4,6; R=1,2,4;
one temporal HYP step, alpha=(1,1,0.5). Original unsmeared field supplies the clover.

## Result and scope

All 996 paired profile records (1992 raw rows) were checked. The real loop and loop–probe averages
change by at most 1.8432e-7 relative to piter=40. Vacuum subtraction increases relative sensitivity:
the largest change in the single-configuration ratio diagnostic is 0.15851% across individual
axes and 0.022700% for the three-axis combined records. Absolute changes remain below 7.2675e-8.
These are numerical comparisons on one configuration, with no statistical error estimate.

The near agreement of these bulk averages does not establish near-SU(3) piter=4 links over
the whole field or ensemble accuracy. Keep piter=4 as the advisor-requested working candidate;
locate the already measured full-field tails, then compare the same observables on further configurations.

## Checks

- ZIP integrity and all expected keys passed: 332 pairs per tau, 83 per axis, four axis labels.
- Source counts are 262144 per individual axis and 786432 for axis=-1 (combined).
- Configuration, Nv, time endpoints, midpoint coordinates, finite values and piter/mode metadata passed.
- Every exported difference and fractional change matches an independent recomputation with the original float convention.
- S_RAW is exactly identical at the printed precision in every pair. The supplied job log also reports unchanged original probe samples.
- Combined records agree with the arithmetic mean of the three axes to floating-point reduction precision.
- The job log reports a fresh piter=40 tau=4 control matching the previous reference exactly.

The CSV consistency checks are independent of the existing validators. Halo/direct-sample and
fresh-control checks are evidence from the supplied job log; raw per-source data were not supplied.

## Vacuum-subtracted diagnostic

For each configuration-level, spatially averaged record, define

    C_p = LS_RE_p / L_RE_p - S_RAW_p

and Delta C = C_4 - C_40. This is a ratio of real parts, not the real part of a complex ratio.
The detailed CSV additionally reports both components of the alternative complex diagnostic

    Ccomplex_p = (LS_RE_p + i LS_IM_p) / (L_RE_p + i L_IM_p) - S_RAW_p.

The largest absolute change in Re(Ccomplex) is 7.2751e-8, consistent with the real-parts result.
Computation uses 50-digit Decimal arithmetic on the printed inputs to avoid introducing further
cancellation error; this does not increase the measurement precision.

For a later ensemble estimate, average the raw L, LS and vacuum quantities over configurations
before forming the ratio. An average of these single-configuration C values is a different estimator.
Ensemble uncertainties and autocorrelations have not been evaluated here.

## Complete scan, all axes and offsets

| tau | Max relative change L_RE | Max relative change LS_RE | Max abs Delta C | Max relative change C (%) |
|---:|---:|---:|---:|---:|
| 2 | 8.051176e-08 | 9.872636e-08 | 3.158362e-08 | 0.008408 |
| 4 | 6.848263e-08 | 8.665821e-08 | 5.674542e-08 | 0.158507 |
| 6 | 1.581224e-07 | 1.843114e-07 | 7.267465e-08 | 0.114235 |

## Three-axis combined records, all offsets

| tau | Max relative change L_RE | Max relative change LS_RE | Max abs Delta C | Max relative change C (%) |
|---:|---:|---:|---:|---:|
| 2 | 4.006382e-08 | 4.881785e-08 | 1.513415e-08 | 0.003080 |
| 4 | 3.447607e-08 | 4.322841e-08 | 2.205502e-08 | 0.011838 |
| 6 | 1.389759e-07 | 1.581676e-07 | 3.739231e-08 | 0.022700 |

## Absolute changes and small denominators

| Quantity | Max absolute change, all 996 records | Max fractional change |
|---|---:|---:|
| L_RE | 3.261021e-15 | 1.581224e-07 |
| LS_RE | 6.574020e-15 | 1.843114e-07 |
| L_IM | 7.446450e-16 | 6.610692e-04 |
| LS_IM | 1.719979e-15 | 4.103094e-02 |

The largest fractional change in LS_IM is about 4.1%, where its reference is very close to zero.
Its maximum absolute change across the scan is 1.7200e-15. Imaginary parts remain finite measured
fluctuations; their relative changes alone do not establish a bias or a statistical significance.

The maximum relative C change occurs at tau=4, axis=1, R=4, (long,tr1,tr2)=(5,1,0):
C40=9.0936863823e-6, C4=9.0792722388e-6, Delta C=-1.4414143473e-8.
The maximum absolute C change occurs at tau=6, axis=0, R=4, (-1,1,0):
C40=-8.9381166390e-4, Delta C=-7.2674654719e-8.
Each table column has its own maximizing record; the maxima do not necessarily coincide.

All relative changes use absolute piter=40 denominators. Values and RMS across coordinates
are descriptive diagnostics: overlapping sources and repeated loop values make these rows
correlated, and they are not independent configurations or statistical error samples.

## Existing full-field SU(3) evidence

| Diagnostic over 50069504 physical temporal links | piter=4 | piter=40 |
|---|---:|---:|
| Max Frobenius norm of Udag U - I | 1.3828922261 | 7.4985014987e-16 |
| RMS Frobenius defect | 4.4105126671e-3 | 2.6576002778e-16 |
| Max absolute determinant defect | 1.7457549354 | 7.7718709383e-16 |
| Links with Frobenius defect > 0.01 | 15678 (0.0313125%) | 0 |

Round-trip and zero-alpha controls, unchanged spatial/absent-boundary links, and the
independent piter=4 repeat passed. The existing global data do not locate the tails in time.
The attached follow-up measures them separately for each physical time slice and for
the actual transporter windows 97..98, 96..99 and 95..100. It does not assume a boundary explanation.

## Separate eigenvector performance comparison

The supplied latest queue snapshot still has 8-node job 5915130 pending for Priority.
Keep it queued. Previous n110 measurements give 877 s for the 16-node pilot and 559 s
for the 32-node full-write repeat, with a small timing-scope difference (the pilot includes its checker).
Using 48 allocated physical cores per node, these are approximately 187.093 and 238.507
physical-core-hours: about 1.569x speedup at 27.48% more allocated physical-core-hours.
These are capacity costs, not measured busy CPU time or confirmed project billing.
The 8-node timing is needed to choose 4x8, 2x16 or 1x32 allocation packing; a packed pilot
must also verify independent node placement and concurrent filesystem performance.

Observable job timings (roughly 89–90 s per piter=4 run and 333 s for the fresh piter=40 control)
are a different workload and must not be mixed with eigenvector scaling.

## Reproduce locally

Python 3.9 or later; standard library only:

```bash
python3 analyze_projection.py hyp_projection_results_5915129.zip reproduced
```

The destination must not exist. Outputs are vacuum_subtracted_diagnostics.csv (996 rows)
and summary.json (overall, combined, per-tau and per-tau/axis/R summaries and input hashes).
The original six-file result archive is included. No HDF5 data or eigenvectors are included.
