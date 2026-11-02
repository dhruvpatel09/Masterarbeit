# Correlated exponential fits of the qX static correlators

All **71 configurations n100–800 in steps of 10**, including n640, enter every
fit. This is post-processing of the temporal-extension measurements published
at `b4a9d6a1c76ccfe41f3b36fb6a255a4fc43485b3`; no new gauge-field measurements
are required. The original measurements and earlier analyses remain the source
of the numerical inputs.

## Setup and scope

- qX: `192 x 64^3`, open time boundaries, periodic spatial directions.
- Fixed source time `t0/a = 90`, measured separations `tau/a = 0,...,32`.
- One temporal HYP step, `alpha = (1,1,0.5)`, `piter = 10`.
- Spatial source modes: APE20, `alpha = 0.5`; cutoffs `Nv = 40,60,80,100`.
- Distances `R/a = 1,2,4,8,12,16,20,23,24,25,28,32`.
- Configuration measurements are already averaged over the specified sources
  and axes. The statistical sample size is **71**, not the number of sites,
  source positions, axes, cutoffs or time slices.
- Thermalization and approximate independence are taken from the ensemble
  assessment. No additional thermalization cut or binning is imposed here.

The primary table is the signed real static correlator, not the connected
clover profile. No logarithm, absolute value, positive-configuration selection,
or signal-based configuration exclusion is used in the fits. The source cutoff
changes the operator; `Nv=100` is not an exact-reference extrapolation limit.

There are **384 one-exponential fits** over all twelve distances and **160
two-exponential fits** at `R/a = 4,12,16,24,32`. Each cutoff/distance uses starts
`tau_min/a = 4,8,12,16` and ends `tau_max/a = 24,32`. All endpoints are inclusive.
Every attempted fit is exported; no preferred potential is selected.

## Models and covariance

With `d = (tau-tau_min)/a`, the dimensionless models are

```
C1(d) = s A0 exp(-e0 d)
C2(d) = s [A0 exp(-e0 d) + A1 exp(-(e0 + delta) d)]
e0 = a E0,  delta = a (E1-E0).
```

`s` is the maximum absolute full-sample mean over the fitted interval. It is a
fixed numerical normalization; amplitudes in the exports refer to `tau_min`.
Amplitudes, `e0` and `delta` are nonnegative. The energy and gap ceiling of 10
is a numerical guard, not a physics prior accepted as a measured bound. Near
boundary solutions are explicitly flagged. The forward, positive-amplitude
models are a bulk approximation with open temporal boundaries; no periodic-time
cosh or backward image is assumed. Residual boundary contributions are not
separately measured by this analysis.

The covariance of the configuration mean is `cov(C_cfg)/71`, including all
fitted times. Cholesky whitening is applied through the correlation matrix for
numerical scaling. There is **no eigenmode removal, diagonal approximation or
shrinkage**. A non-positive covariance is an error. A single fit contains at
most 29 time points; this does not invert the rank-deficient joint covariance of
all four cutoffs and all times.

SciPy bounded least squares uses an analytic Jacobian, five initial guesses
for one exponential, and thirteen for two exponentials. The latter includes
the nearly nested one-exponential solution. The best converged candidate is
retained, and its objective must not exceed the one-state objective beyond the
specified numerical tolerance. All 4,000 candidate records are archived in
`results/solutions.json`.

## Errors and parameter diagnostics

Naive errors propagate the full covariance through the **implicit derivative
of the fitted optimum**. The Hessian includes the residual-dependent second
derivative of the model, rather than only a Gauss–Newton approximation. This
is conditional on the full-sample covariance used as the fit weight.

The field `identified` means only that the listed **local parameter checks**
pass: no boundary hit, positive nonsingular local curvature, resolvable
amplitudes and gap, and no sufficiently different lower energies among nearby
successful starts. Exact thresholds appear in `fit_models.py`. In particular,
this flag does not establish a good fit, an asymptotic plateau, correct model
order, or a physical ground-state energy. Parameter solutions and formal local
errors for flagged fits remain in the exports.

A leave-one-configuration-out audit is requested for 32 selected models:
`Nv=40,100`, `R/a=4,12,24,32`, starts `8,12`, end `32`, both model orders.
Twenty-two full-sample fits pass the local checks and are audited, preserving
all 71 omissions and **fixed full-sample covariance and scale**. The other ten
are explicitly marked as unperformed; no configuration is removed. The
jackknife uses a single warm start from the full optimum and therefore checks
the local branch, not an exhaustive global optimization search in every
replicate. Twenty-one audits retain resolved parameters in every replicate.
For these, jackknife/naive energy-error ratios range from **0.99821 to 1.00288**.

At `Nv=100, R/a=24, tau/a=8...32`, the two-state full-sample local checks pass,
but five of the 71 jackknife solutions fail local parameter checks. Its formal
jackknife error is retained with `jackknife_reliable=False`.

Cutoff differences and two-minus-one-state energy differences use **paired**
configuration fluctuations and, where available, matched jackknife omissions.
Their small errors can resolve operator or model sensitivity even when
individual estimates overlap. They do not include model-truncation or
window-selection uncertainty. No multiple-comparison discovery claim is made.

## Goodness-of-fit calibration

The known-covariance Gaussian chi-square tail is only a reference, because the
time covariance is estimated from the same finite ensemble. A Hotelling
reference based on the local interior-model residual dimension is also
exported. It is an approximation for these nonlinear fits and is suppressed
when the local checks fail.

Eight one-state fits, `Nv=40,100`, `R/a=4,12,24,32`, `tau/a=8...32`, additionally
use a recentered whole-configuration bootstrap. Each empirical configuration
is shifted by fitted mean minus observed mean; 399 matched resamples then
re-estimate the covariance and refit under that null distribution. The common
resampling indices and seed are archived. Every draw succeeds. The upper-tail
estimate is `(1 + exceedances)/400`. Reported Monte Carlo errors describe
finite-draw uncertainty only; finite-sample bootstrap bias is not quantified.
These tests are exploratory and are not an automatic fit acceptance rule.

The Gaussian, linearized Hotelling and empirical bootstrap tails differ
substantially. For example, `Nv=40, R/a=24, tau/a=8...32` gives approximately
`1.16e-7`, `0.0080`, and `0.0825 +/- 0.0138`, respectively. The last uncertainty
is only the sampling uncertainty of 399 bootstrap draws. No single number in
this comparison establishes ground-state isolation. Increasing model order is
not assigned a chi-square likelihood-ratio p-value: the one-state limit lies
on a boundary with an unidentified extra-state energy.

## Main numerical observations

The following are one-exponential **fit estimates**, with conditional naive
errors; they are not a selected final static potential. The endpoint is 32.

| R/a | Nv | start 8: aE | start 12: aE |
|---:|---:|---:|---:|
| 4 | 40 | 0.116784 +/- 0.000092 | 0.116732 +/- 0.000124 |
| 4 | 100 | 0.116884 +/- 0.000090 | 0.116742 +/- 0.000128 |
| 12 | 40 | 0.194350 +/- 0.000499 | 0.193203 +/- 0.000607 |
| 12 | 100 | 0.195682 +/- 0.000472 | 0.194453 +/- 0.000641 |
| 24 | 40 | 0.265004 +/- 0.001774 | 0.263586 +/- 0.002867 |
| 24 | 100 | 0.519349 +/- 0.019130 | 0.541720 +/- 0.082956 |
| 32 | 40 | 0.341065 +/- 0.018493 | 0.251794 +/- 0.025682 |
| 32 | 100 | 0.296888 +/- 0.003421 | 0.295407 +/- 0.007937 |

At `R/a=4`, later one-state energies approach a common value across cutoffs.
At `R/a=12`, the paired `Nv40-Nv100` difference over `tau/a=12...32` is
`-0.001250 +/- 0.000329`; residual source/model dependence matters at the
current precision. The much larger cutoff difference at `R/a=24` and the
window dependence at `R/a=32` prevent a uniform ground-state interpretation
across the spatial scan. Small statistical errors or acceptable fit curves
alone do not remove that issue.

Only **50 of 160 two-state fits pass the local parameter checks**. Later starts
frequently remove enough information that the excited amplitude or gap cannot
be separated. This is a constraint of the present data, not evidence that the
missing state has zero physical contribution. The exported curves can remain
stable when the two fitted energies are unresolved.

## Figures

1. [One-state window and cutoff scan](results/single_exponential_window_scan_71configs.png):
   all four cutoffs at six representative distances, endpoint 32.
2. [One versus two states](results/one_vs_two_state_71configs.png):
   energy bars only for locally resolved two-state fits; unresolved cases have
   crosses at the foot of the panel, not at a measured zero energy.
3. [Signed correlators and whitened residuals](results/signed_correlator_model_checks_71configs.png):
   direct fits at R/a=12,24. A symmetric-log axis retains negative means and
   error bars. Cholesky residuals mix time slices and are not independent
   point-by-point pulls.
4. [Covariance calibration](results/goodness_covariance_calibration_71configs.png):
   the same eight fits under three reference/calibration choices.
5. [Two-state parameter resolution](results/two_state_parameter_resolution_71configs.png):
   all 160 models; colours describe parameter checks, not fit acceptance.

Every figure also has an SVG version. Tables include all distances and windows,
not just the displayed cuts. `effective_energy_reference_comparison.csv` links
each one-state fit to the earlier six-effective-energy constant fit with the
same starting time. Their time support differs, so no significance is assigned
to the difference between those estimates.

## Reproduction and integrity

From the repository root, Python with NumPy and SciPy suffices for numerical
verification. Matplotlib is required only to redraw the figures.

```bash
python3 laplace_static/analysis/qx_correlated_exponential_71configs/test_analysis.py
python3 laplace_static/analysis/qx_correlated_exponential_71configs/reproduce.py --verify
python3 laplace_static/analysis/qx_correlated_exponential_71configs/verify_files.py
```

The numerical verifier reloads all 112,464 combined primaries from the published
raw table, re-evaluates every archived optimization candidate and covariance,
and replays 1,562 jackknife and 3,192 bootstrap solutions and the paired exports.
It does **not** rerun the optimizer search or claim bit-identical optimization
across SciPy versions. The nine tests check the model derivatives, whitening,
synthetic recovery, implicit error propagation, normalization, signed inputs,
singular covariance rejection, boundary flags and null recentering.

`verify_files.py` needs only the Python standard library. It checks every package
file and both published dependencies by SHA-256. It is suitable for the
SuperMUC Python 3.6 installation without adding scientific Python packages.

To rerun the complete optimization and redraw:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 laplace_static/analysis/qx_correlated_exponential_71configs/run_analysis.py
python3 laplace_static/analysis/qx_correlated_exponential_71configs/derived_results.py
python3 laplace_static/analysis/qx_correlated_exponential_71configs/render_plots.py
```

Results can vary slightly across numerical-library versions or in poorly
identified models. After an intentional regeneration, review the new exports
and run `make_manifest.py` in this analysis directory to regenerate the package
manifest; the archived manifest protects the published results. Python 3.6 syntax is supported by the core scripts; optional
rendering depends on the installed compatible Matplotlib version.

## Next physics checks

Use time-window and source stability together before choosing an energy for
each distance. A focused common-energy fit to several source operators would
be a useful next test, retaining their cross covariance and respecting the
71-sample covariance rank. Long-distance energies require particular care:
the additional-state parameters can be unresolved even when a curve fits the
means. Only after a defensible energy range is established should a correlated
potential/Cornell analysis be promoted beyond diagnostics. This package does
not quote a final string tension or a physical-unit scale.

## References

- [SciPy `least_squares`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html):
  bounded trust-region least squares with callable derivatives.
- N. H. Christ, C. Eranki, C. Kelly,
  [Bootstrap-determined p-values in Lattice QCD](https://arxiv.org/abs/2409.11379),
  Phys. Rev. D 111, 074514 (2025): recentered bootstrap calibration of fit statistics.
- M. Bruno, R. Sommer,
  [On fits to correlated and auto-correlated data](https://arxiv.org/abs/2209.14188),
  Comput. Phys. Commun. 285, 108643 (2023): fit and covariance diagnostics.
