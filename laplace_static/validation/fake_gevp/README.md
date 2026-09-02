# Fake-data GEVP validation

## Purpose

This validation exercises the numerical GEVP implementation used by
`laplace_static/scripts/analyze_laplace_gaussian_gevp.py` on an externally
supplied synthetic correlation matrix.  It is separate from the Gaussian
profile physics analysis.  A successful result demonstrates that the solver can
recover well-resolved levels from a suitable operator basis; it does not by
itself establish the hidden generator energies.

The validator deliberately calls the production functions
`hermitian`, `metric_diagnostics`, `metric_status`, `solve_gevp`,
`central_and_jackknife_matrices`, `jackknife_error`, and
`effective_from_principal`.  Generalized eigenvalues are independently checked
against `scipy.linalg.eigh`.

## Input provenance and layout

The raw input is kept locally inside the ignored validation-data directory:

```text
laplace_static/validation/fake_gevp/data/fakedat.mat
```

It is not committed because it is an externally supplied 12 MB binary input.
The checksum below identifies the exact file required for reproduction.

The validated SHA-256 is

```text
0ad8a271954fa47559e7d307a2bb40983a5d97a63fee4829bcd34f27ddbff6be
```

The MATLAB-v5 file contains one `4 x 4` cell array named `corr`.  Every cell is
a real `1000 x 101` array: 1000 measurement-history entries and times
`t = 0,...,100`.  MATLAB column `t+1` corresponds to Python time index `t`.

For each measurement and time, the validator forms

```text
C(t) <- (C(t) + C(t)^dagger) / 2
```

before solving the GEVP.  This is necessary because `corr{i,j}` and
`corr{j,i}` have the same expectation value but independently generated noise.

## Statistical treatment

The primary errors use a contiguous delete-one-block jackknife with block size
10, giving 100 jackknife blocks.  A block-size scan over
`1, 2, 5, 10, 20` diagnoses the autocorrelation correction.  The validator also
writes the median and 5--95% range of the autocorrelation functions of all
`4 x 4 x 101 = 1616` primary input streams.

The primary reference time is `t_ref=5`.  Reference times
`0, 1, 2, 3, 5, 8, 10` are all analyzed.  The choice `t_ref=5` is a diagnostic
working point that leaves an early-time interval for the higher levels while
retaining a well-conditioned metric; the final interpretation must use the
full reference-time scan.

No expected energy is supplied to the program.  Its mechanical correlated-fit
scan writes all complete windows and marks those with `p >= 0.05`; those rows
are candidates for inspection, not automatic final plateau selections.

## Numerical acceptance gates

The run is successful only if:

- all selected central metrics pass the production positive-definiteness and
  `rcond=1e-10` rule;
- all selected block-jackknife metrics are valid;
- no central or block-jackknife GEVP solve fails;
- production and SciPy generalized eigenvalues agree within the recorded
  floating-point tolerance;
- generalized eigenpair residuals and metric-orthonormality residuals pass
  their thresholds;
- both the existing production self-test and the validation-driver analytic
  self-test pass.

Non-positive principal correlators at late times are not repaired.  Their
effective energies remain invalid and are shown as diagnostic/incomplete
points where a central result still exists.

## Run on Stromboli

From the repository root:

```bash
cd /home/m2130292/Masterarbeit

python3 -B laplace_static/scripts/validate_gevp_fakedat.py --self-test

sbatch laplace_static/scripts/submit_validate_gevp_fakedat.sh
```

The launcher obtains the repository root from
`SLURM_SUBMIT_DIR`.  Slurm executes a copied script below its spool directory,
so resolving the root from `BASH_SOURCE[0]` would incorrectly select the spool
tree.  Submit the job from the repository root as shown above.

No Slurm array or explicit wall-clock limit is needed.  The job uses one task,
one CPU, and 2 GB of memory on `compute2011`.

## Generated evidence

The job writes the following under
`laplace_static/validation/fake_gevp/results/`:

- `fake_gevp_validation.txt`: human-readable audit and numerical gates;
- `fake_gevp_metadata.json`: machine-readable provenance and parameters;
- `fake_gevp_effective.csv`: principal effective energies, errors, and validity;
- `fake_gevp_conditioning.csv`: metric spectra and jackknife support;
- `fake_gevp_autocorrelation.csv`: primary-stream autocorrelation summary;
- `fake_gevp_fit_scan.csv`: correlated constant-fit window scan;
- `fake_gevp_effective_tref5.png`: four principal effective-energy channels;
- `fake_gevp_reference_time_scan.png`: reference-time stability;
- `fake_gevp_autocorrelation.png`: autocorrelation diagnostic;
- `fake_gevp_block_error_ratios.png`: blocked/unblocked error comparison;
- `fake_gevp_metric_conditioning.png`: reference-metric conditioning.

## Validated numerical result

The canonical run completed as Stromboli job `156982` on 2026-09-02.  It used
Python 3.8.10, NumPy 1.24.4, and SciPy 1.10.1.  Every selected central metric
and every block-jackknife metric was accepted, and no GEVP solve failed.  At
the primary choice `t_ref=5`, block size 10, the normalized metric condition
number is `31.8827`, with `100/100` valid jackknife metrics.

The production solver agrees with `scipy.linalg.eigh` to a maximum absolute
difference of `7.77e-15`.  The largest generalized-eigenpair residual is
`5.95e-16`; metric-orthonormality residuals are of order `1e-14`, far
below the `5e-10` acceptance threshold.

The longest mechanically accepted correlated constant-fit windows are:

| state | time window | energy | chi2/dof | p-value |
|---:|---:|---:|---:|---:|
| 1 | `10..29` | `0.10004281(2405)` | `1.390` | `0.119` |
| 2 | `14..30` | `0.19989901(21289)` | `1.465` | `0.102` |
| 3 | `6..19` | `0.30038316(11739)` | `1.689` | `0.056` |
| 4 | `6..22` | `0.40438741(18499)` | `1.497` | `0.091` |

The central reference-time curves are mutually compatible over their common
early-time ranges.  The unit-diagonal metric condition number remains moderate,
from `30.5` to `68.6`, over `t_ref=0,1,2,3,5,8,10`.  The primary-stream
autocorrelation has median integrated time `tau_int=0.971`; blocking increases
the median effective-energy errors by about `1.32--1.38` at block size 10,
with little further change at block size 20.  This supports block size 10 as
the primary error estimate.

States 1 and 2 remain usable far into the measured time range.  States 3 and
especially 4 become noise dominated later; non-positive principal correlators
then produce intentionally invalid effective energies.  Their early-time
plateaux remain resolved.

These checks establish that the production implementation recovers four
well-separated levels from this suitable four-operator matrix.  Confirmation
of the numerical levels against the hidden generator values remains an
independent blind check.

After the Python validation is interpreted, selected effective-energy points can
be cross-checked independently with the existing Octave-compatible `UWerr.m`
Gamma-method implementation.  That statistical cross-check should be performed
after choosing the diagnostic time points without using the hidden generator
values.

## Interpretation boundary

Passing the numerical gates validates the GEVP implementation and its use on
this external dataset.  Extracting four stable energy levels additionally
requires compatible plateaux under changes of time window, reference time, and
block size.  The extracted values should then be sent to the data generator for
blind confirmation.
