# qX mixed-source matrices: 71 configurations

This package analyzes the completed piter=10 matrix campaign, including all 71
configurations n100,...,n800 in steps of ten. It extends the technical pilot
n170/n640/n650 without counting those configurations twice. The source and sink
operators use nested Laplacian cutoffs Nv=40,60,80,100 from the same saved 100
eigenvectors. These cutoffs are operator choices; piter is a separate numerical
SU(3)-projection iteration count.

The matrix is indexed **sink row, source column**. Every ordered source/sink pair
is retained, including complex off-diagonal entries. The statistical sample is
71 configuration averages. Spatial sites, the three axes, source choices and
cross-validation folds are not additional independent configurations.

## Measurement and validation

The unchanged native executable measured R/a=4,12,24,32, at source time t0/a=90
and temporal separations tau/a=0,...,32. qX has geometry 192 x 64^3 with open
time and periodic space. R/a=32 is half the spatial box; the temporal half-extent
is 96, not 32. The last measured correlator is C(32), hence the last adjacent
effective energy uses C(31)/C(32), plotted at midpoint 31.5.

Cold replay verifies every archived member, native build receipt, preparation,
case marker, stdout/stderr hash and numerical CSV against the original pinned
checker code. All 71 saved-mode audits and 72 matrix scans are replayed (71 p10
scans plus the separate n650 p40 control). This includes:

- 599,808 p10 matrix rows and 872,448 spatial rank-block rows;
- 2,343 audited time slices and 234,300 saved eigenmodes;
- 149,952 exact p10 diagonal matches to the published temporal-extension data;
- 109,056 direct mixed-contraction checks, 27,264 legacy-diagonal checks,
  1,136 equal-time Gram checks and 213 periodic round trips;
- maximum saved-mode relative residual 1.5843e-11;
- maximum central-link unitarity/determinant defects below 8.0e-16;
- n650 p10-versus-p40 maximum matrix-component difference 6.618e-24.

The diagonal-only mean, independent standard error, effective energy and
jackknife error regressions reproduce all 528 available correlator points and
512 effective-energy entries for these four radii with the earlier 50-digit
Decimal arithmetic. Undefined points retain their status. No positive-sample
filter or absolute-value logarithm is used.

Native executables and two libraries are not present in the review archive;
their four hashes are retained and checked against the build manifest. Replay
does not execute them or reopen the multi-terabyte gauge/eigenvector inputs.
This is a replay of archived diagnostics, not a new independent MPI computation.

The replay permits platform roundoff only when comparing the independently
recomputed `scans[i].minimum_equal_time_relative_pivot` with its archived
summary value. The absolute comparison bound is 32 binary64 machine epsilons
plus 5e-12 times the larger magnitude; both values must be finite, nonnegative
and in the same zero/positive class. This diagnostic is derived through a
four-by-four Cholesky calculation and can differ in its last digits across
math-library/Python builds. The native Gram, positivity and all other checks
are executed unchanged before this comparison. Other summary fields, case
metadata, primary/rank-block CSV values and SHA-256 receipts remain exact.
A mismatch now identifies the configuration, nested field and compared values.

At the same-sample training normalization point Q(0)=1 by construction, and
its true sampling error is zero. Roundoff-sized jackknife errors give an
ill-conditioned signal/error ratio of order 1e14--1e15. For this one point per
candidate, replay checks both normalized complex means and both component
errors within 64 binary64 epsilons, verifies the internal signal/error ratio,
and compares reciprocal signal/error (relative uncertainty) with that roundoff
floor. The positive-infinite signal convention is allowed only if its error
is exactly zero. This allowance never applies at nonzero time or to held-out
normalizations. All physical-time signal checks and numerical estimates retain
their existing comparisons. The archived numerical exports are unchanged;
these normalization-point ratios do not measure physical precision.

## What is tested offline

The equal-time correlation metrics are positive definite for the full sample
and every delete-one sample in each of three bases:

| Basis | Laplacian cutoffs |
|---|---|
| pair | 40,100 |
| triple | 40,60,100 |
| full | 40,60,80,100 |

Diagonal scaling is applied before whitening. The full-basis condition numbers
are 1226.94, 455.70, 14.09 and 8.24 at R/a=4,12,24,32. The short-distance
operators are more redundant, but none meets the 1e-12 relative metric-rank
rejection threshold. No metric eigenmode is silently removed.

For a transparent source-selection test, at training time tau/a=4,6 or 8,
choose v to maximize

    v^dagger H(t_train) v / (v^dagger C(0) v),
    H(t_train) = [C(t_train) + C(t_train)^dagger]/2.

Normalize v^dagger C_training(0) v=1 and hold v fixed at every subsequent and
earlier evaluation time. This Hermitian training objective chooses a source;
it does not replace measured correlators with their Hermitian parts and does
not establish a positive spectral decomposition or variational energy bound.
Every measurement is Q_i(t)=v^dagger C_i(t) v using the original complex C_i.

Two procedures are compared for all 36 radius/basis/training-time candidates:

1. **Same sample:** train on all 71 configurations and evaluate on all 71.
2. **Cross fit:** train on alternating configuration-index fold A and evaluate
   only fold B, and vice versa; average all 71 held-out measurements. Fold A
   has 36 configurations and fold B has 35. Cross fit combines two sources;
   it is not itself a single fixed physical operator.

Every one of 71 outer delete-one configuration resamples repeats the metric,
source training, both holdout folds, averaging and nonlinear log ratio. Thus
the quoted cross-fit errors do not incorrectly treat overlapping trained
sources as independent observations. Contrasts share the same deletions.
The raw projection resamples are archived as numerical arrays.

The displayed energy is log(Re mean Q(t)/Re mean Q(t+1)), only for positive
endpoint means. The full complex means, phases and errors remain exported.
Phase magnitude 0.05 radians and real mean/uncertainty 3 are review flags, not
physics acceptance criteria. No failed sample is dropped. For this data all
2,304 projected energy points and all their 71 deletion sets are positive;
that algebraic fact alone does not establish a plateau.

Six-point constant fits are catalogued in windows 6--11, 10--15, 14--19,
18--23 and 22--27 for all candidates. Full-sample outer-jackknife covariance
defines GLS weights, held fixed in the deletion fits. The p values use an
approximate Gaussian chi-square reference with estimated covariance; many
choices were examined. They are not automatic acceptance rules or adjusted
discovery significances. Full covariance across 16 matrix entries x 33 times
is not inverted: 71 samples imply empirical covariance rank at most 70.

General complex pencils C(t)u=lambda C(t_ref)u, with reference times 0,2,4,
are also exported as central algebraic checks. They retain nonsymmetric,
complex matrices. Their eigenvalues are unordered roots, not tracked physical
states. All backward residuals must be below 1e-10. The replay matches whole
spectra and circular phases, rather than relying on LAPACK eigenvalue order.

## Main observations

The following numerical examples use the full four-source basis, training time
4 and the held-out procedure. Errors include full source retraining.

| R/a | Constant over 14--19 | Approximate local p | Minimum endpoint signal |
|---|---|---|---|
| 4 | 0.1165708 +/- 0.0001943 | 0.685 | 601.9 |
| 12 | 0.1934663 +/- 0.0011090 | 0.743 | 118.5 |
| 24 | 0.2664523 +/- 0.0049142 | 0.098 | 26.1 |
| 32 | 0.3150254 +/- 0.0148040 | 0.662 | 7.91 |

These are examples from a catalogue, not four selected ground-state energies.

- **R/a=4:** the pair and full source bases give very similar effective
  energies. The 14--19 fitted full-minus-pair contrast is
  -0.00000345 +/- 0.00000609. Expanding a nearly redundant basis offers little
  clear improvement in this comparison. The early-to-middle window shift is
  about two paired errors, so local flatness still needs a window check.
- **R/a=12:** the signal remains strong, but window dependence is visible.
  The 14--19 minus 6--11 fitted difference is
  -0.0032736 +/- 0.0012454 (a pointwise diagnostic, not a global significance).
- **R/a=24:** at correlator time 24 the held-out full-basis signal/uncertainty
  is 15.84, versus 0.504 for the original Nv=100 diagonal and 14.78 for Nv=40.
  Thus source choice removes much of the poor Nv=100 overlap, but the existing
  low-cutoff operator already supplies most of this gain. At 14--19 the
  optimized-minus-Nv40 constant is only 0.0008170 +/- 0.0012790. At time 32
  the projected signal falls to 2.28 and its phase reaches about -0.70 radians.
  A full-range extrapolation from the earlier region is not justified.
- **R/a=32:** the full basis improves some signal measures, but the late
  apparent plateau is not stable. The 22--27 constant is
  0.1555700 +/- 0.0386091, compared with 0.2983989 +/- 0.0037700 in 6--11.
  Their paired shift is -0.1428289 +/- 0.0392684, despite a local late-window
  p of 0.616. This illustrates why a flat noisy patch is insufficient.
- Held-out and same-sample results are close for the early training choices
  examined, with nonzero paired uncertainty retained. This single alternating
  split is a useful overfitting check, not proof that every optimized result
  is unbiased. Training-time and basis comparisons are exported in full.
- Nonzero-time matrices become less Hermitian as signal decreases. At the
  half spatial box the measured matrix is real to roundoff, yet can still be
  nonsymmetric and can have complex generalized roots. Algebraic diagnostics
  alone cannot distinguish statistical fluctuations from boundary/operator
  systematics. A standard positive Hermitian GEVP interpretation is not assumed.

No final static potential, string tension, Cornell fit or excited-state level
is claimed by this package. Before that step, the source, time-window and
open-boundary controls must support the intended interpretation. R=4 and 12
are useful high-signal controls, R=24 tests source overlap and phase stability,
and R=32 remains a large-distance/window-sensitivity control.

## Figures and numerical exports

The v2 distribution revises figure styling. Pair-basis curves use blue hollow
circles and solid lines; full-basis curves use orange filled squares and dashed
lines. Triple-basis curves use magenta hollow downward triangles and dotted
lines. Nv40 and Nv100 reference curves use grey upward triangles/dotted lines
and purple diamonds/dash-dot lines, respectively. Reference-time and
training-time comparisons also use distinct colours and shapes. Colours,
markers, fills and line patterns provide redundant distinctions when curves
overlap or figures are printed in greyscale. All numerical CSV, JSON and NPZ
results and all statistical methods are unchanged from the original distribution.

Every figure has PNG and SVG versions under `results/`:

| Figure | Reading guide |
|---|---|
| [matrix_conditioning_71configs.png](results/matrix_conditioning_71configs.png) | Equal-time metric spectra; delete-one ranges are not confidence intervals. |
| [matrix_complexity_71configs.png](results/matrix_complexity_71configs.png) | Relative non-Hermitian and imaginary norms; central diagnostics only. |
| [fixed_source_effective_energies_71configs.png](results/fixed_source_effective_energies_71configs.png) | Nv40/Nv100 diagonals and pair/full held-out sources, all times. |
| [heldout_effective_energies_detail.png](results/heldout_effective_energies_detail.png) | Pair/full sources only, times 6--21; reference diagonals remain in companion figure. |
| [heldout_signal_and_phase_71configs.png](results/heldout_signal_and_phase_71configs.png) | Signal and phase expose unreliable long-time regions. |
| [paired_source_sensitivity_71configs.png](results/paired_source_sensitivity_71configs.png) | Paired anchor/basis shifts, with complete-retraining errors. |
| [holdout_comparison_71configs.png](results/holdout_comparison_71configs.png) | Held-out minus same-sample effective energies. |
| [constant_window_stability_71configs.png](results/constant_window_stability_71configs.png) | Window sensitivity, phase/signal flags and local-fit diagnostics. |
| [complex_pencil_spectra_71configs.png](results/complex_pencil_spectra_71configs.png) | Unordered complex roots; not a physical energy spectrum. |

The CSVs include all entries, gaps, review flags, paired contrasts and every
attempted constant window. `projection_resamples.npz` contains the complete
projection means and 71 outer resamples for further paired analyses.

## Installation and reproduction

The original review ZIP is 124,921,131 bytes, above GitHub's ordinary 100 MiB
per-file limit. `import_review.py` checks its full SHA-256 and imports four
chunks (three 32 MiB plus one 24,257,835-byte chunk). These preserve the exact
original ZIP bytes and stay below the per-file limit. The small distribution
ZIP contains the scripts, plots and tables; source chunks must be imported from
the original SuperMUC review. Do not stage the oversized original ZIP itself.

From the repository root on Stromboli, after copying the original review ZIP
and unpacking this result package:

```bash
qx_analysis=laplace_static/analysis/qx_static_matrix_71configs
python3 "$qx_analysis/import_review.py" --archive qx_static_matrix_full71_review.zip
python3 "$qx_analysis/verify_files.py"
OPENBLAS_NUM_THREADS=1 python3 "$qx_analysis/test_analysis.py"
OPENBLAS_NUM_THREADS=1 python3 "$qx_analysis/reproduce.py" --verify
```

`verify_files.py` needs only the standard library and supports Python 3.6.
Scientific replay requires NumPy and SciPy; figure regeneration also requires
Matplotlib. Existing scientific environments may be used. No pip installation
or environment replacement is performed by these scripts.

The cold replay uses a temporary folder under `work/` and removes it on exit.
Allow several minutes for complete hashing, extraction and native stdout replay;
no Slurm allocation is needed. Requirements and generation versions are archived.
The package depends on the already-published temporal-extension raw CSV and two
result CSVs, and on the tracked original matrix-pilot review. All four are pinned.

For deliberate regeneration, `replay.py` rebuilds the primary cache,
`run_analysis.py` rebuilds tables and resamples, `render_plots.py` rebuilds
figures, and `make_manifest.py` refreshes checksums. Do not refresh receipts to
silence an unexplained verification mismatch.

When committing, include all four chunks plus both NPZ files. The package's
local ignore exceptions retain them; `git ls-files --error-unmatch` can confirm
that all six binary inputs/results are staged. Keep the original review until
both clusters validate the same published commit. Transfer bundles can be
removed after their destination validates successfully; native campaigns,
builds and saved eigenpairs remain the underlying evidence.

The next step is scientific review of source/window/boundary stability using
these offline exports. No new production measurement or MPI job is submitted
by this package.
