# Shared-energy and source-operator diagnostics for qX

Every fit uses all **71 configurations n100–800 in steps of 10**, including
n640. The measurements are the published temporal-extension correlators, with
`piter=10`, fixed source time `t0/a=90`, `tau/a=0,...,32`, and mode cutoffs
`Nv=40,60,80,100`. No further gauge-field measurements are required for this
analysis. Dependencies are pinned to the results available at
`62011935cbe1612cb55e2ad9fe0b318fd46a19a5`.

The question is whether different source operators support a common energy,
with independent amplitudes, over the measured time intervals. A mode cutoff
changes the source operator; these cutoffs are not an extrapolation sequence
with `Nv=100` as an exact reference. Different finite-window energy estimates
can reflect different mixtures of states, even though the physical spectrum
should be common in the bulk approximation.

## Data and fit designs

The qX lattice is `192 x 64^3`, with open time boundaries and periodic spatial
directions. Spatial modes use APE20 at alpha=0.5; temporal transport uses one
HYP step at alpha=(1,1,0.5). The statistical unit is one configuration after
the specified source and axis averages. The configurations are treated as
approximately independent. No spatial position, axis, cutoff or time slice is
counted as an additional independent configuration.

The focused joint fits cover `R/a=4,8,12,16,24,32`, starts `8,12,16`, and end
`32`, with inclusive endpoints. They complement the previous full-distance
analysis rather than replacing its twelve measured distances.

| Design | Cutoffs | Times |
|---|---|---|
| pair_dense | 40,100 | every 2 units from the start through 32 |
| pair_matched | 40,100 | every 4 units from the start through 32 |
| all_four | 40,60,80,100 | the same every-4 grid as pair_matched |

These are **54 joint cases**, each with three model choices: **162 fits**.
The complete covariance across the selected sources and times is used. Its
largest dimension is 28, below the maximum sample-covariance rank of 70.
The explicit time-grid control avoids silently regularizing a covariance with
too many components. There is no SVD cut, shrinkage, or diagonal replacement.
The maximum correlation-matrix condition number in this scan is about 37,681.

All signed configuration measurements are retained. No logarithm, absolute
value, positive-configuration filter or signal-based configuration exclusion
enters the fits. A per-source numerical normalization is fixed from the maximum
absolute ensemble mean over its fitted times. This changes the amplitude
coordinates, not the energy or objective.

## Models

With `d=(tau-tau_min)/a`, the normalized model for source n is one of

```
shared_one:       C_n(d) = A_n exp(-e d)
independent_one:  C_n(d) = A_n exp(-e_n d)
shared_two:       C_n(d) = A_n exp(-e0 d) + B_n exp(-(e0+gap) d).
```

Here `e=aE`. The second model relaxes only the common-energy constraint while
retaining the same source/time covariance. The third assumes two common
energy levels with independent amplitudes for each source. Its gap is common,
not a separate effective excited energy for every operator. This is an
explicit two-level approximation, not a complete spectral representation.

Amplitudes, energies and gap are nonnegative. An energy/gap ceiling of 10 is a
numerical guard; boundary solutions are flagged. Forward exponentials are a
bulk approximation with open time boundaries. Periodic-time cosh terms are
not assumed, and residual temporal-boundary contributions are not separately
quantified here.

The optimizer uses analytic first derivatives and multiple starts: five for
shared_one, six for independent_one, and fourteen for shared_two. The two
larger models include a nearly nested shared_one start; their best successful
objectives cannot exceed the shared_one objective beyond numerical tolerance.
All 1,350 full-sample candidate solutions are archived.

## Errors and local parameter checks

Naive errors propagate the independent-configuration covariance through the
implicit derivative of the actual optimum. The curvature contains the
residual-dependent model Hessian. These errors are **conditional on the
full-sample covariance used as the fit weight**. They do not include model
order, time-window selection, or covariance-weight uncertainty beyond that
conditional propagation.

The field `local_checks_pass` tests boundaries, curvature, resolvability of
amplitudes and gap, and sufficiently different near-optimal starts. It is not
a test of goodness of fit or physical ground-state isolation. Some source
amplitudes can be poorly measured while an energy is more stable; failure of
the complete parameter check does not by itself prove the energy is wrong.
Formal estimates and flags are retained for every model.

Only **5/54 shared_two models pass all full-sample local checks**. All five
start at 8. The parameter audit includes the original selected start-12
shared_one/shared_two models, and additionally every shared_two model passing
the full-sample checks. This gives **29 requested audits**, of which **17**
have locally resolved full-sample parameters and are performed. The other
twelve remain explicitly unperformed, without removing any configuration.

Each audit omits the same one of 71 configurations, uses fixed full-sample
covariance and source normalization, and starts on the full-sample local
branch. All 1,207 refits converge. Fifteen audits preserve all local parameter
checks in every replicate. Among those, jackknife/naive error ratios range
from **0.99949 to 1.01436**. Two extra-state audits have failures of local
checks: pair_dense R16/start8 has 65/71 resolved replicates; all_four R4/start8
has 70/71. Their formal errors remain exported with the audit flag unset.
Warm-start audits do not exhaustively search all branches in every replicate.

Source-energy differences within independent_one and shared_two-minus-
shared_one contrasts preserve the covariance of the 71 matched configuration
fluctuations. The exported `paired_z_diagnostic` is a conditional linearized
diagnostic, not a separately bootstrap-calibrated probability or a
multiple-testing discovery claim.

## Two different bootstrap questions

Every design at start 12 is calibrated: **18 cases, 399 whole-configuration
resamples each**, using archived common resampling indices. Each empirical
configuration is shifted by fitted shared_one mean minus observed mean.
Every resample re-estimates the full covariance and refits shared_one and
independent_one, retaining matched source/time measurements.

The resulting null distribution answers two separate questions:

1. Does the common one-state model describe the data overall? The statistic is
   its full correlated objective.
2. Is the common-energy constraint compatible with the data within the
   one-state family? The statistic is the improvement when source energies
   are allowed to differ, using the same covariance.

All **7,182 resamples and 14,364 null refits succeed**. Upper-tail estimates
are `(exceedances+1)/400`. The tables include exceedance counts and exact
95% binomial intervals for finite-draw uncertainty conditional on the empirical
ensemble. These intervals do not quantify finite-sample calibration bias or
model uncertainty. The minimum numerical resolution is 1/400; zero
exceedances would not establish an exact zero probability.

Known-covariance Gaussian references are also exported, clearly labelled.
The extra-state model is not assigned a chi-square likelihood-ratio p-value:
the nested one-state limit has boundary amplitudes and an unidentified gap.
No bootstrap threshold automatically selects a plateau or a preferred model.

## Numerical observations

The following common one-state estimates use all_four at start 12 and end 32.
They are **conditional fit estimates**, not selected final potential points.
The last column tests the common-energy constraint in that same four-source
one-state family.

| R/a | Common aE +/- naive SE | Bootstrap tail for common-energy constraint |
|---:|---:|---:|
| 4 | 0.116589 +/- 0.000123 | 0.9375 |
| 8 | 0.165665 +/- 0.000362 | 0.5425 |
| 12 | 0.193719 +/- 0.000664 | 0.0950 |
| 16 | 0.216283 +/- 0.001182 | 0.0775 |
| 24 | 0.261222 +/- 0.002832 | 0.0175 |
| 32 | 0.301945 +/- 0.010313 | 0.0550 |

At R4 and R8 the later source/time-grid comparisons support a common-energy
description at the present precision. A high tail probability does not prove
ground-state isolation. At R12, the focused Nv40/Nv100 test gives 0.0125 on
the dense grid and 0.0100 on the matched coarser grid. The effect therefore
persists under that time-grid change. The four-source test imposes three
energy equalities instead of one and has a different covariance dimension;
its omnibus tail of 0.095 does not replace the focused comparison.

Within the four-source unrestricted fit at R12, the paired Nv40-Nv100 energy
difference is **-0.0009148 +/- 0.0002393**. Overlapping individual energy bars
can conceal this sensitivity. At R24 the full common one-state goodness tail
is 0.035 and the common-energy-constraint tail is 0.0175, both exploratory
calibrations. A constrained estimate near 0.26 alone does not settle the
large-distance ground-state problem. R32 also retains source and window
uncertainty and weak signal for some cutoffs.

## Equal-time source shape and long-time signal

The additional source-shape export uses all twelve measured distances and
all 71 configurations. It contains ratios of ensemble means, with paired
naive errors, for

```
<C_N(R,0)> / <C_N(1,0)>
<C_N(R,24)> / <C_N(R,0)>.
```

Among the measured distances, the Nv100 equal-time shape has its minimum at
R/a=25, while the Nv80 sampled minimum is at 28. For Nv100 the first ratio is
0.006583 +/- 0.000011 at R25 and rises to 0.010298 +/- 0.000043 at R32.
The minimum changes with the source cutoff; it is not a feature of the
physical static potential inferred from this plot.

At R24 and tau24 the second ratio is 0.001324 +/- 0.000090 for Nv40 but only
0.000040 +/- 0.000078 for Nv100. This directly illustrates the differing
long-time signal. Equal-time shape is an operator diagnostic, **not a direct
ground-state overlap measurement**. The pattern motivates further study of
source filtering and state overlap; it does not by itself establish the
mechanism or a physical spectral effect. Negative means remain visible.

## Figures and tables

- [Shared energy versus fit start](results/shared_energy_window_scan_71configs.png)
- [Bootstrap source compatibility and absolute fit checks](results/source_compatibility_bootstrap_71configs.png)
- [Shared one versus two states](results/shared_one_vs_two_state_71configs.png)
- [Source-dependent energies in the joint fit](results/source_dependent_energies_71configs.png)
- [Equal-time source shape and signed long-time signal](results/source_shape_and_long_time_signal_71configs.png)

All five have SVG versions. The fit, contrast and resampling tables retain
every case, including unresolved models. Crosses at the foot of model-order
panels mark unresolved parameters, not measured zero energies. Distinct plots
use distinct normalization or fit quantities; their ordinates should not be
equated.

## Reproduction

From the repository root, with NumPy and SciPy available:

```bash
python3 laplace_static/analysis/qx_shared_energy_71configs/test_analysis.py
python3 laplace_static/analysis/qx_shared_energy_71configs/reproduce.py --verify
python3 laplace_static/analysis/qx_shared_energy_71configs/verify_files.py
```

The fifteen tests check derivatives, exact model recovery including four
sources, cross-source covariance, implicit propagation, normalization/order
invariance, nesting, signed data, singular covariance rejection, matching of
configuration/source/time indices, bootstrap recentering, finite-draw tails,
paired ratio cancellation, and summary roundoff, structure and nonfinite-value
checks. The numerical verifier reloads the complete
112,464-primary input table and replays all archived fits, audits, null refits
and derived exports. Optimizer searches are not repeated by `--verify`.

Computed floating-point fields in the bootstrap summary are compared with
relative tolerance `1e-10` and absolute tolerance `1e-13`, allowing numerical
library roundoff in probability functions. Configuration labels, case order,
field names, counts and other discrete values must match exactly. Nonfinite
values are rejected. The file-integrity checks remain exact; numerical
verification does not rewrite archived results or plots.

`verify_files.py` needs only the Python standard library and checks package
bytes and four published reference/dependency files. It is suitable for the
SuperMUC Python 3.6 environment. Core scripts use Python 3.6 syntax; numerical
and plotting libraries must be installed in compatible versions.

To rerun the full optimization and redraw, use

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 laplace_static/analysis/qx_shared_energy_71configs/analysis.py
python3 laplace_static/analysis/qx_shared_energy_71configs/render_plots.py
```

After intentionally regenerating and reviewing results, run `make_manifest.py`
in this directory. Numerical minima in poorly identified models can vary
across library versions; archived parameter solutions and file hashes remain
the reference for the published package.

## Next measurement decision

The remaining issue is source/state control before final potential extraction.
Use both the source-shape and fitted-state diagnostics when designing a
focused source-operator pilot. Correlations between existing measurements are
already retained here; an operator cross-correlator matrix for a variational
analysis is a different measurement. Its construction and validation need to
be checked in the measurement code before a production campaign. The existing
71-configuration mode set can be reused. This package does not quote a final
Cornell potential, string tension or physical-unit scale.

## References

- [SciPy bounded least squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).
- N. H. Christ, R. Eranki, C. Kelly,
  [Bootstrap-determined p-values in Lattice QCD](https://arxiv.org/abs/2409.11379),
  Phys. Rev. D 111, 074514 (2025): recentered bootstrap calibration.
- M. Bruno, R. Sommer,
  [On fits to correlated and auto-correlated data](https://arxiv.org/abs/2209.14188),
  Comput. Phys. Commun. 285, 108643 (2023): covariance and fit diagnostics.
