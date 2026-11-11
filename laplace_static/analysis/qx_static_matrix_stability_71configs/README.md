# qX matrix uncertainty and source stability

This is an **offline follow-up** to the full-ensemble mixed-source analysis at
commit `2b9830b3e5de49b8234b2257c24a0175f247cc61`. It uses the same validated
measurements. No new MPI campaign, eigenpair computation or clover measurement
is part of this package.

The questions are: do the complex matrix defects exceed their uncertainties;
does the held-out split determine the result; and do correlated fits to the
correlators confirm the effective-energy window behaviour?

## Data and statistical unit

- All **71 configurations**, n100 through n800 in steps of ten, including n640.
  Configuration independence is the working assumption from the ensemble
  information supplied by the advisors. These diagnostics do not independently
  establish thermalization or independence.
- Lattice `192 x 64^3`: open temporal and periodic spatial boundaries. The
  spatial half-box separation is `R/a=32`; half the temporal extent is 96.
- Original complex matrices at `R/a=4,12,24,32`, correlator times `tau/a=0..32`,
  source time `t0/a=90`; sink index is the row and source index the column.
- Source cutoffs `Nv=40,60,80,100`; the numerical projection iteration count
  remains `piter=10`. Cutoff and projection iteration count are different controls.
- The unit of resampling is a complete configuration, retaining its axes,
  matrix entries, times and source choices together. Neither axes, sites,
  folds nor repeated split choices increase the statistical sample size.

The previous package's primary cache was reconstructed from the validated
native outputs: 599,808 matrix rows and 872,448 rank-block rows, with 149,952
exact diagonal-reference matches. The original source archive remains in that
package's tracked chunks. This package pins eight published dependencies by
SHA-256 and byte count, rather than duplicating those chunks. Its data loader
also checks the previous manifest's published Git blob identity.

## 1. Complex components with uncertainties

For each configuration, form `A=(C-C^dagger)/2`. An anti-Hermitian 4x4 matrix
has 16 independent real components: four diagonal imaginary components and
the real and imaginary parts of its six upper off-diagonal entries. Each is
scaled by the fixed full-sample equal-time metric diagonal for readability.
This scale cancels in its studentized statistic.

The fixed test family includes times 4,8,12,16,24,32. At the periodic half box,
the matrix is analytically real; its ten imaginary anti-Hermitian components
are checked against a `64*machine_epsilon` bound and excluded from significance
tests. The six real antisymmetric components remain. There are **324 tests**.
All times and all matrix imaginary entries are also exported descriptively.
An imaginary off-diagonal entry alone does **not** imply non-Hermiticity.

The global-null check recenters the 71 component vectors, resamples whole
vectors 1,999 times, and recomputes each bootstrap studentizer. The maximum
absolute t over the entire family preserves its cross-component correlations.
There is no inversion of a 324-dimensional covariance from 71 observations.
Seeds and draw indices are archived; degenerate studentizers are errors,
not silently discarded draws. Adjusted tails and the global tail use the
`(exceedances+1)/(draws+1)` convention.

Results:

- Maximum observed component `|t| = 3.1696`.
- Bootstrap 95% global maximum threshold `3.9048`.
- Global empirical null tail **0.3475**; no component has a family-adjusted
  tail below 0.05.
- `SUMMARY.json` includes a binomial interval for finite-draw uncertainty,
  approximately `[0.3263,0.3685]`, conditional on the empirical centered distribution. It does not measure
  finite-ensemble bootstrap calibration bias.

Thus the defect plots alone do not establish a statistically resolved
anti-Hermitian mean in this prescribed family. **Non-rejection is not proof of
Hermiticity, spectral positivity or absence of open-boundary effects.** It
does not authorize replacing the measured matrices by their Hermitian parts.

## 2. Additional held-out source splits

The source-selection algorithm is unchanged: a Hermitian Rayleigh objective
at one early training time, normalized by the equal-time metric, chooses a
source vector. Evaluation always uses the **original complex matrix**.
This is a fixed-source control, not a physical positive-spectrum GEVP analysis.

The bases are pair `(40,100)`, triple `(40,60,100)` and full `(40,60,80,100)`,
with training times 4,6,8. Four balanced random splits (seeds 1701--1704) are
added to the original alternating 36/35 split. Each fold is evaluated using a
source trained on its complementary fold. The same-sample result remains a
control. Together these give **216 projected curves**.

Every held-out curve is recomputed under all 71 configuration deletions,
including metric estimation and source retraining. Original fold membership
is retained after deletion. The alternating central curves and all their
deletions reproduce the published cache. Same-sample deletion curves are
reused from that pinned, already verified cache.

For the full basis, training time 4a, and the 14--20 fit window, the five
held-out fitted-energy ranges are:

| R/a | Minimum aE | Maximum aE |
|---:|---:|---:|
| 4 | 0.1165710 | 0.1165731 |
| 12 | 0.1934787 | 0.1934888 |
| 24 | 0.2668795 | 0.2670834 |
| 32 | 0.3172129 | 0.3177925 |

These split changes are small compared with the individual energy errors.
Paired random-minus-alternating changes at this setting are at most 1.34
pointwise standard errors. The five choices share the same data: their spread
is a sensitivity diagnostic, not a new estimate based on five independent
ensembles. All bases and training times are retained, without choosing a
winner by its statistical tail.

## 3. Correlated fits to projected correlators

Fit the real ensemble mean to

`Re Q(tau) = A exp[-aE * (tau/a - tau_min/a)]`

in five seven-point windows: 6--12, 10--16, 14--20, 18--24 and 22--28. The
amplitude is signed and the energy search interval is `[-0.5,1.5]`; positivity
is not imposed on individual measurements or fit parameters. Fits to the
real part are diagnostic ansatz tests, not a demonstration that the imaginary
part vanishes.

Use the complete-retraining jackknife covariance of each mean curve, scaled
to a correlation matrix before Cholesky whitening. A relative covariance
eigenvalue threshold of `1e-12` rejects deficient systems. There is no
pseudoinverse, diagonal fallback or regularization. The full-sample covariance
is fixed during all 71 parameter refits. This propagates configuration and
source-training fluctuations, conditional on that covariance prescription.
It does not calibrate uncertainty from estimating the covariance itself.

The amplitude is eliminated analytically. Energy minima are bracketed by
profile-gradient changes on a 161-point grid, refined with Brent root solves
and bounded Newton updates, and compared with both interval endpoints.
Interior solutions must satisfy an energy stationarity bound of `1e-10`.
All **1,080 fits** are numerically defined, with no central or deletion energy
at a search boundary. Numerical resolution alone is not physical acceptance.

The exported Gaussian chi-square tails use estimated covariance matrices and
are **uncalibrated, exploratory** diagnostics. In particular, a large local
tail does not establish a plateau. The component max-t bootstrap above is a
separate null test and does not calibrate these fit tails.

Each seven-correlator fit is compared with the corresponding six-point
effective-energy constant: correlators t through t+6 correspond to log ratios
with lower endpoints t through t+5. Their difference is paired across identical
configuration deletions. Time-window, source-basis and training-time contrasts
are paired in the same way. A small error on a difference between highly
correlated estimators is not a small error on the energy itself.

### Representative results, not selected final potentials

Full basis, training time 4a, original alternating held-out split; parentheses
below are complete-retraining jackknife standard errors:

| R/a | Correlator window 6--12 | 14--20 | 22--28 |
|---:|---:|---:|---:|
| 4 | 0.117029 (0.000134) | 0.116571 (0.000194) | 0.116665 (0.000257) |
| 12 | 0.196742 (0.000674) | 0.193479 (0.001109) | 0.197949 (0.002826) |
| 24 | 0.271310 (0.002208) | 0.267016 (0.004831) | 0.292235 (0.018712) |
| 32 | 0.298572 (0.003800) | 0.317732 (0.015582) | 0.162659 (0.038882) |

Interpretation by radius:

- **R/a=4:** a strong signal and a promising stable region, but the paired
  14--20 minus 6--12 change is `-0.000458 +/- 0.000230` (about 2 pointwise SE).
  Do not assign zero time-window systematic uncertainty.
- **R/a=12:** the same change is `-0.003264 +/- 0.001245` (2.62 pointwise SE).
  High statistical precision still permits resolved time dependence.
- **R/a=24:** the 14--20 and early fits agree within their paired error. The
  22--28 window has a maximum central phase about 0.281 rad, a Gaussian fit
  tail about 0.010, and only about 4 SE minimum real signal. Its direct fit
  minus effective-energy constant differs by `0.02190 +/- 0.00666`.
  Late-window phase, shape and nonlinear-estimator sensitivity remain concerns.
- **R/a=32:** the late-minus-early change is `-0.13591 +/- 0.03956` (3.44
  pointwise SE). Nevertheless the late-window Gaussian tail is about 0.526.
  The late fitted energy remains near 0.163 in all five held-out splits.
  Thus neither a satisfactory local tail nor a different split removes the
  time-window dependence. Half-box reality does not establish a plateau.

These pointwise contrasts belong to a large diagnostic catalogue and are
not multiple-test-adjusted discoveries. Full numerical exports remain the
authoritative values; there is no automatic potential-point selection here.

Source-basis sensitivity should not be confused with split sensitivity. At
training time 4a in the early 6--12 window, the paired full-minus-pair energy
change is `+0.002512 +/- 0.000288` at R24, and `-0.009066 +/- 0.002453` at R32.
The R24 effect is particularly well resolved as a paired contrast even though
it is small compared with the energy itself. At R32 in the late window the
full-minus-pair change is `-0.07990 +/- 0.05298`, much less precisely resolved.
Thus stable fold choices do not establish stability under source design or
time-window changes. A larger basis or a lower fitted energy alone is not a
criterion for physical acceptance.

## Figure guide

All figures have PNG and SVG versions under `results/`. Pair sources use
blue hollow circles and solid lines; the full basis uses orange squares and
dashed lines. Triple-basis controls use magenta triangles and dotted lines.
Small paired changes are plotted directly rather than hidden under overlaid
central curves.

1. `antihermitian_family_71configs`: component excursions against the common
   family threshold; distinguishes central defects from resolved deviations.
2. `projected_phase_uncertainty_71configs`: phases with jackknife errors and
   pointwise imaginary-component diagnostics. Half-box roundoff is excluded.
3. `heldout_split_energy_shifts_71configs`: random-minus-alternating paired
   energies, demonstrating the scale of split sensitivity.
4. `correlator_fit_window_scan_71configs`: fitted energies as the time window
   moves. Phase/signal guides and approximate-tail flags remain visible.
5. `correlator_vs_effective_energy_71configs`: paired estimator differences,
   including the late R24 discrepancy.
6. `paired_basis_and_training_checks_71configs`: source-basis and training-time
   sensitivity, keeping the two controls separate.
7. `signed_correlator_fit_checks_71configs`: measured real correlators and
   model shapes, retaining signs on a symlog axis. Only a fixed display
   normalization is used; no absolute values or logarithms hide negative data.
8. `paired_time_window_changes_71configs`: early-to-middle and early-to-late
   changes with paired errors, showing why a plateau remains to be established.

The measured coverage remains through tau=32; effective energies end at lower
endpoint 31 (midpoint 31.5). The above fits end at 28. The full-time phase and
signed-signal plots retain the remaining measured times.

## Verification and regeneration

From the repository root, with NumPy/SciPy installed:

```bash
qx_analysis=laplace_static/analysis/qx_static_matrix_stability_71configs
python3 "$qx_analysis/verify_files.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 "$qx_analysis/test_analysis.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 "$qx_analysis/reproduce.py" --verify
```

`verify_files.py` uses only the standard library and is appropriate for the
SuperMUC receipt check. Numerical analysis belongs on Stromboli or another
analysis machine; no Slurm allocation is required. Python files use Python
3.6-compatible syntax; installed numerical-library versions must support the
chosen interpreter. The generation environment is archived in `inputs/`.

`reproduce.py --verify` actually recomputes the 1,999 bootstrap draws, five
held-out splits, all their deletions, and all 1,080 correlated fits with 71
parameter refits each. Allow several minutes, depending on the machine.
It checks labels, cohort, seeds, validity masks, CSVs, resample arrays and
summary. Floating comparisons allow narrow platform roundoff, reject
infinities, and require matching undefined-value masks. Covariance conditions
are checked in reciprocal form, while stationary roots are checked against
their accuracy bound. Same-sample Q(0) is a training normalization; its real
signal-to-error ratio is intentionally blank, avoiding significance from
roundoff in a theoretically zero error.

All source procedures, including held-out splits, also have a real equal-time
quadratic form because C(0) is Hermitian. Their tiny imaginary mean and
imaginary jackknife error are retained as raw numerical diagnostics, but the
ratio of those two roundoff values is undefined and therefore blank. Central
and every deletion projection must first satisfy a `64*machine_epsilon`
imaginary-reality bound relative to their equal-time real scale. A substantive
violation is an error, not a masked test. Held-out real equal-time uncertainty
is retained, and nonzero-time imaginary diagnostics remain unchanged. This
avoids a spurious replay mismatch between BLAS implementations at tau=0
without weakening tolerances for physical quantities. The corresponding
phase plot omits the undefined equal-time imaginary-significance point.

For intentional regeneration, run `run_analysis.py`, then `render_plots.py`
and `make_manifest.py`. Generation rewrites results; verification does not.
Both resampling NPZ files are explicitly retained despite broad repository
ignore patterns. Keep the previous raw-review chunks and primary cache.

## Remaining work

This completes the offline uncertainty/split/correlator-fit follow-up. It does
not produce final physical GEVP levels, a final static potential, a Cornell
fit, string tension, or new flux-tube profiles. A single fixed source time
cannot settle open-boundary systematics. The absence of a family-adjusted
anti-Hermitian rejection is insufficient to assume a positive spectral form.

The next choice is a measurement design review: preserve the usable short-
radius information, and decide whether R24/R32 need source-time averaging or
explicit source-time/boundary controls before further spectroscopy. Additional
matrix radii or source designs can be considered once that question is fixed.
No broad new SuperMUC submission is requested by this package. Any later
production measurement should follow a small validated control and retain
configuration-level pairing.

Method references: [Bruno and Sommer, arXiv:2209.14188](https://arxiv.org/abs/2209.14188)
and [Christ, Eranki and Kelly, arXiv:2409.11379](https://arxiv.org/abs/2409.11379).
The latter motivates empirical calibration; its fit-statistic bootstrap is
not claimed to have been performed for the one-exponential fits here.
