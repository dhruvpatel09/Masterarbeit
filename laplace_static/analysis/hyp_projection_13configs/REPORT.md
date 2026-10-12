# qX paired HYP profiles: 13 configurations n110–230

All 39 configuration/time-window cases were rechecked from the three original
review archives. The newest n200–230 campaign completed all 12 cases with
7,968 raw and 3,984 paired records. Across all batches there are 25,896 raw and
12,948 paired records; these produce 996 ensemble profile positions, of which
249 combine the three axes. n100 is an exact regression control and is excluded
from this ensemble. The source CSV files are preserved byte-for-byte in inputs/.

## Ensemble comparison

At each position and projection prescription, the estimator is

    C = mean_configuration(LS_RE) / mean_configuration(L_RE)
        - mean_configuration(S_RAW).

Configurations have equal weight. Raw measurements are averaged before taking
the ratio. Individual-configuration connected ratios are not averaged. Arithmetic
uses 50 decimal digits to protect the vacuum subtraction from cancellation.
The same configuration is omitted in both prescriptions for each delete-one
jackknife replicate; this preserves their covariance. The CSVs contain both
the control statistical error and the paired error of the projection difference.
Reported central values are the full-sample plug-in estimator, without a jackknife
bias correction. This definition is used consistently for the 5/9/13 subsets.

| Temporal separation tau | Maximum combined-axis relative shift | Maximum shift / preliminary control SE |
| --- | ---: | ---: |
| 2 | 0.0009325% | 0.0001335 |
| 4 | 0.0016373% | 0.0001210 |
| 6 | 0.0033033% | 0.0004194 |

Each maximum is taken over all 83 combined-axis probe positions for that tau;
the two maxima need not occur at the same position. The relative shift is
100*abs(C4-C40)/abs(C40). All ratios are defined for this dataset.
Across the individual axes as well, the largest ensemble relative shift is
0.0148668%.
The much larger single-configuration, individual-axis relative differences
(largest 16.654% in n220, tau=6)
must be interpreted with their connected-signal denominators; they are not an
ensemble bias estimate. The ensemble CSV also reports absolute differences,
control signal/error and matched difference errors.

These results support the use of piter=4 for further measurements in the tested
central windows while retaining piter=40 controls. They do not supply an automatic
physics acceptance gate or establish an ensemble-wide bound on projection bias.
The sample is only 13 ordered configurations. Autocorrelation has not been
assessed; the delete-one errors are preliminary, not final uncertainty estimates.
The data at different probe positions are correlated and are not independent
configurations. No effective sample size or significance of maxima is inferred.

## SU(3) diagnostics

Metrics are ||Udagger*U-I||_F and |det(U)-1| over physical temporal links in the
actual transport windows. The absent outgoing upper-boundary link is excluded.
For the latest four configurations, maximum piter4 defects are
0.00044664056 and
0.00023117719, on n230.
Across all 13 configurations the previously observed n170, tau6 maximum remains:
unitarity 0.0062409825, determinant
0.003164495. This is visibly above machine precision.
No tested central-window link exceeds the diagnostic 1e-2 count threshold;
that threshold is not an acceptance criterion. The piter40 control remains
below 8e-16 in both metrics. The original-field clover measurements remain
exactly identical between projection prescriptions in every archived case.

Global boundary-tail results from the earlier n100 study are not erased by
central-window insensitivity. This report covers midpoint98, tau2/4/6, R1/2/4,
Nv100, one temporal HYP step, alpha=(1,1,0.5). It does not demonstrate that
four projection iterations return all links to SU(3) at strict accuracy.

## Verification and reproduction

VERIFICATION.json records the 39 case identities, provenance and diagnostics.
The archived package/prepared manifests, controller ledger and case exports were
checked. Structural checkers were rerun on the actual MPI outputs. Their JSON and
CSV outputs matched the archived exports exactly, including the n100 regression.
The static science sources and regression tables match across all three batches.
Runtime libraries, gauge files and eigenvector HDF5 data are not in the reviews;
they were not independently reread for this archive analysis.

From Masterarbeit:

```bash
python3 laplace_static/analysis/hyp_projection_13configs/test_analysis.py
python3 laplace_static/analysis/hyp_projection_13configs/analyze.py
# Optional; requires matplotlib:
python3 laplace_static/analysis/hyp_projection_13configs/plot.py
```

Tests cover the ratio-of-means definition, analytical jackknife normalization,
shared-noise cancellation, and rejection of duplicate/missing records, changed
clover measurements and modified input bytes. An independent math.fsum calculation
checks the ensemble estimator on all 996 positions. plots use preliminary control
errors as a descriptive scale, with explicit autocorrelation limitations.

## Next production step

The accompanying mental/production/qx_evec_batch_240_270 package prepares four
new eigenvector runs on four independent eight-node groups in one 32-node general
allocation. It preserves the successful n200–230 implementation and settings.
Temporal HYP is separate from eigenvector generation. After generation, audit
saved modes before extending the same paired profile measurements.
