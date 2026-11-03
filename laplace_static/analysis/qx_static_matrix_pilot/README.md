# qX mixed-source pilot: replay, conditioning and exploratory spectra

This package analyzes job array **5961383**, using configurations **170, 640
and 650**. It completes the offline review of the mixed-source MPI validation
pilot at commit `646175fd325c07ec5aafcb8cbb0322be4c3ba9dc`.
It contains no production submission commands.

## Validated measurements

The original review ZIP is retained byte for byte in
`inputs/job5961383_review.zip`, with its SHA-256 in `inputs/SOURCES.json`.
`replay.py` securely extracts it temporarily, verifies preparation and build
bindings, case markers, launch inventories, stdout/stderr hashes and exported
CSV hashes, and replays the original strict saved-mode and contraction checks.
All archived numerical tables are also compared directly with the parsed logs.

| Check | Result |
|---|---:|
| Configurations / audits / scans | 3 / 3 / 4 |
| Audited saved modes | 9,900 across 33 times per configuration |
| Complex matrix records | 33,792 |
| Complex rank-block records | 49,152 |
| p10 diagonal reference matches | 6,336; maximum difference **0.0** |
| n650 p40 diagonal difference from p10 reference | 3.31e-24 maximum |
| n650 full p10/p40 component difference | 6.62e-24 maximum |
| Independent rectangular / legacy checks | 6,144 / 1,536 |
| Equal-time Gram checks / periodic payload roundtrips | 64 / 12 |
| Largest saved-mode relative residual | 8.24e-12; unchanged limit 1e-9 |
| Largest central-link unitarity / determinant defects | 7.68e-16 / 7.80e-16 |

The replay does not open the physical gauge/eigenvector paths or re-execute
native binaries. Build artifact hashes are retained; files present in the
review are verified, and omitted binaries/libraries are explicitly listed in
`results/VALIDATION.json`. Their original SuperMUC validation remains part of
the recorded build provenance.

All seven MPI launches returned zero and their stderr files are empty. The
batch-script stderr records the already observed sticky-module unload messages.

| Configuration | Audit | p10 scan | p40 control | Total recorded MPI time |
|---|---:|---:|---:|---:|
| 170 | 3.94 min | 3.57 min | — | 7.51 min |
| 640 | 4.11 min | 3.57 min | — | 7.68 min |
| 650 | 3.96 min | 3.55 min | 6.99 min | 14.50 min |

These totals exclude scheduler delays and Python preparation/replay overhead.

## Matrix definition and conditioning

The source/sink cutoffs are independently `40,60,80,100`; separations are
`R/a=4,12,24,32` and `tau/a=0..32`, with source time `t_src/a=90`.
All original complex ordered entries are retained. Analysis uses the recorded
three-axis average, with **sink rows and source columns**. The ensemble mean
is an unweighted mean over the three configurations.

Equal-time conditioning is measured after normalization
`A(tau) = D^-1 C(tau) D^-1`, where `D_ii = sqrt(C_ii(0))`.
Thus operator normalization does not dominate the quoted condition numbers.

| R/a | Four cutoffs 40/60/80/100 | Three cutoffs 40/60/100 | Two cutoffs 40/100 |
|---|---:|---:|---:|
| 4 | 1192.97 | 410.05 | 125.57 |
| 12 | 442.21 | 175.23 | 48.55 |
| 24 | 13.64 | 8.46 | 1.18 |
| 32 | 8.53 | 3.34 | 2.61 |

Every equal-time metric is positive definite for the examined configurations,
their mean and leave-one-out means. The nested cutoffs have strong correlations
at short distance. The quoted condition numbers do not by themselves invalidate
double-precision eigensolutions. They motivate explicit basis-size and
reference-time comparisons. No metric modes are removed automatically.

Finite-time matrices are assessed with singular values and the relative
anti-Hermitian norm `||A-A^dagger||_F / (2 ||A||_F)`. The latter grows strongly
for some long-distance, late-time matrices. A real half-box correlator at
`R=32` can still be nonsymmetric in its source/sink labels.

## Exploratory generalized pencils

The numerical problem is

\[
C(\tau)v_n = \lambda_n(\tau,\tau_{ref}) C(\tau_{ref})v_n.
\]

The equal-time metric is diagonally scaled and whitened to improve numerical
coordinates. Reference separations `0,1,2,4` are examined. All finite-time
pencils are solved using the **general complex solver**, preserving their
nonsymmetric entries. The Hermitian-part eigenvalues are exported separately
as `hermitian_part_diagnostics.csv`; they do not replace the measured matrices.

At the first time beyond the reference, modes are ordered by decreasing real
eigenvalue. Later labels follow maximum right-eigenvector overlap in the
equal-time metric. Leave-one-out modes are matched to the central modes in the
central equal-time metric. Vector phases are fixed canonically for portable
replay. The largest generalized-eigenpair relative residual is below 7e-16 in
the recorded run; the numerical acceptance limit is 1e-10.

Some late-time conjugate modes have equal or nearly equal overlaps. Different
LAPACK column order or roundoff can then give the same complete spectrum
different tracked labels, particularly at R=32. Portable replay uses the
archived mode labels as a naming convention: every pencil is solved again,
then the freshly computed complex eigenvalues are assigned to the archived
labels by minimum spectral distance. The recomputed values, eigenvectors,
residuals, overlaps, projections and logarithmic ratios remain independent
calculations. Material spectral differences still fail verification.

Principal phases and imaginary parts of complex logarithmic ratios are
compared using circular distance, so equivalent +pi and -pi representations
agree. Both values must remain in the principal range. All other numerical
tolerances, source hashes and saved-mode/projection checks remain unchanged.
This portability correction does not alter any archived input, numerical table
or figure. Tests explicitly exercise reordered eigenpairs, arbitrary vector
phases, roundoff-sized imaginary sign changes and materially incorrect spectra.

The leading four-source result is sensitive to later reference times. Smaller
bases can reduce this sensitivity. The reference-time plot retains these
failures and gaps; selecting a favorable reference does not establish a
ground-state plateau.

The standard Hermitian GEVP spectral interpretation assumes an appropriate
positive transfer-matrix representation and control of finite-time effects
[Blossier et al., arXiv:0902.1265]. Fixed source time and open temporal boundaries
also require boundary-state control. The present small-sample nonsymmetric
pencils are numerical feasibility diagnostics; their complex branches are not
identified as physical energy levels.

## Fixed-source and held-out checks

A leading vector determined at `tau/a=4`, with reference separation zero, is
held fixed in the contraction `q(tau)=v^dagger C(tau)v` at every later time.
The result is compared with individual Nv=40 and Nv=100 correlators and with
two- and three-source bases. This avoids independently retuning the displayed
source at every plotted time.

For each omitted configuration, the vector is retrained on the other two and
evaluated both on that training pair and on the omitted configuration. The
held-out configuration never enters its training mean. These curves test source
sensitivity with minimal data; they do not provide three independent precision
estimates or confidence intervals.

Complex logarithmic ratios are exported in full. The real effective-energy
figures show a point only when both correlators have positive real part and
absolute phase at most **0.05 radians**. This is an explicitly labelled display
gate, not a statistical acceptance test. Negative and complex results remain
in the CSVs. Plots show the first 16 adjacent-time ratios; all 32 are exported.

The equal-time matrix rank and the statistical covariance rank answer separate
questions: the four-source metric has rank four, while three configuration
samples imply empirical covariance rank at most **two**. Consequently this
package performs no correlated energy fits, confidence-interval claims, final
potential extraction or string-tension fit.

## Results and figures

- `matrix_conditioning_pilot.(png,svg)`: basis-size conditioning, metric
  eigenvalues, finite-time non-Hermiticity and singular-value conditioning.
- `complex_principal_correlators_pilot.(png,svg)`: all four tracked complex
  branches at reference zero; signed real and imaginary parts use symlog axes.
- `fixed_source_effective_energies_pilot.(png,svg)`: individual sources and
  fixed-vector projections; gaps follow the declared phase/positive-real gate.
- `reference_time_sensitivity_pilot.(png,svg)`: four-source versus two-source
  leading branches for all four reference times.
- `fixed_source_holdout_pilot.(png,svg)`: vectors trained on two configurations
  and tested on the third; no confidence bands.

`SUMMARY.json` records the principal findings and numerical counts. The eight
CSV exports retain conditioning, metric spectra, original complex generalized
modes, fixed-vector coefficients/correlators, projection controls and separately
labelled Hermitian-part checks. `ANALYSIS_PLAN.json` records the exploratory
choices. Basis and display choices are diagnostic choices, not a blinded
production-selection procedure.

## Reproduce on Stromboli

Use the existing scientific Python environment with NumPy, SciPy and Matplotlib
(`requirements-analysis.txt`). The verification scripts create no Slurm jobs.

```bash
cd /home/m2130292/Masterarbeit
qx_analysis=laplace_static/analysis/qx_static_matrix_pilot
python3 "$qx_analysis/test_analysis.py"
python3 "$qx_analysis/reproduce.py" --verify
python3 "$qx_analysis/verify_files.py"
```

`verify_files.py` is the lightweight standard-library integrity check, suitable
also for SuperMUC. Full numerical replay additionally needs NumPy and SciPy.
Numerical replay uses explicit tolerances for eigensolver roundoff and verifies
residual bounds, preserves the archived mode-label convention and compares
principal phases modulo 2*pi. It does not require floating-point summaries to
match exactly between library versions. File SHA-256 checks remain exact.

To intentionally regenerate numerical tables and plots:

```bash
python3 "$qx_analysis/run_analysis.py"
python3 "$qx_analysis/render_plots.py"
python3 "$qx_analysis/make_manifest.py"
```

## Remaining work

1. Decide production geometry and retain the full matrix as an offline basis
   selection resource. The technical pilot covers four radii. Additional radii
   require corresponding measurements on the pilot configurations as well.
2. Implement full71 preparation, validated pilot reuse, provenance and array
   scheduling. For unchanged geometry, 68 configurations require new measurements;
   the three validated pilot configurations can be reused.
3. Submit that production campaign when ready; no such submission is made here.
4. Assemble all71 matrices and use the agreed independent-configuration treatment
   with paired covariance. Test basis size, reference time, open-boundary effects,
   projection accuracy, mode tracking and time windows before energy fitting.
5. Extract a potential only where the full-ensemble stability checks support it;
   then assess physical-distance coverage and comparisons with the earlier
   diagonal and EM1 analyses.

Retain the original SuperMUC pilot build, campaign directories, mode storage and
recovery paths for later reuse. Transfer ZIPs and bundles can be removed after
their destination has verified the corresponding files.
