# qX source-time pilot: completed review

Offline analysis of the completed native source-time pilot. All three cases
(n170, n640, n650) pass the archived native checks. This package finishes the
technical pilot review and stops before another measurement campaign.

The analysis has no native build, MPI launch or scheduler submission step.
The VS Code spell-checker edit can be published in its own commit using
[PUBLISH.md](PUBLISH.md); it is independent of the numerical analysis.

## Question and fixed scope

Does moving the temporal interval change the signed complex matrices and
the projections of exactly the same frozen sources?

- qX geometry: `192 x 64^3`, open temporal and periodic spatial boundaries.
- Source positions: `t0/a=58,90,122`; separations `tau/a=0..32`.
  The largest intervals are 58 to 90, 90 to 122, and 122 to 154. No temporal wrap.
- Radii: `R/a=4,12` are controls; `R/a=24,32` are the main long-distance checks.
  The value 32 is the spatial half-box, not the temporal half-extent (96).
- Every ordered pair of cutoffs `Nv=40,60,80,100` is retained. Sink is the
  matrix row and source is the column; no nonzero-time Hermitian replacement.
- Sources: pair `(40,100)`, triple `(40,60,100)`, full `(40,60,80,100)`,
  each trained at separations 4,6,8, at each radius: 36 frozen candidates.
- Their weights were chosen at `t0/a=90` using the other 68 configurations,
  excluding all three pilot cases, and are never retuned at shifted positions.
  This is a fixed independent-training control, not the earlier full71
  alternating cross-fit estimator. Values from those procedures should not
  be spliced together as if they were identical operators.

The plotted main comparison uses the pair and full bases trained at separation
4a, fixed before viewing shifted outputs. All candidates remain in the CSV.

## Validation result

Cold replay reconstructs the CSV exports from the original stdout with the
unchanged published checker helpers. It verifies source bytes, preparation and
build receipts, physical endpoints, launch order, completion markers, stdout
and stderr hashes, all ordered entries, rank/axis reductions and mode checks.
It also checks the baseline against the published primary matrix cache and
replays the excluded-pilot source training using phase-invariant invariants.
An independent dense contraction checks every archived frozen-source value.

| Check | Result |
|---|---:|
| Validated configurations | 3/3 |
| Native launches | 3 audits + 8 scans; all exit codes zero |
| Audited time slices / eigenmodes | 229 / 22,900 |
| Maximum saved-mode relative residual | 1.0656e-11 (criterion 1e-9) |
| Native matrix / rank-block rows | 67,584 / 147,456 |
| Independent direct / legacy diagonal checks | 12,288 / 3,072 |
| Equal-time Gram checks / periodic roundtrips | 128 / 24 |
| Fresh n650 t0=90 axis-averaged matches | 2,112 exact complex entries |
| Maximum shifted p10/p40 component difference | 3.3087e-24 |
| Maximum axis-averaged p10/p40 relative matrix-norm difference | 1.9819e-15 |
| Archived fixed-source projections / paired contrasts | 13,068 / 9,504 |
| Baseline primary-cache matches / frozen training candidates | 6,336 exact / 36 |

The p10/p40 control compares n650 at identical `t0/a=122` endpoints. It tests
projection-iteration accuracy, not temporal-position independence.

The original review and ten replay dependencies are included. The primary
cache has SHA-256 `6dd357bf702aed84f5373aea634c81ddbf363cf0b8c3b0835ce106e7de0131aa`.
All dependency hashes are pinned to repository snapshot
`69368a64e4f9a4982641c255ae11a2b713e7145f`. Original binaries and some build
inputs are outside the review: their recorded hashes remain in the build
receipt, and their absence is explicitly listed in `VALIDATION.json`.
Replay does not execute MPI or reopen the gauge/eigenvector files.

## What the position comparisons show

Short-distance R4 is a useful control: its position changes are much smaller
than the large-distance changes at R24/R32. At R12 the late-time dependence
and phase fluctuations are more visible. At R24/R32, substantial within-case
changes persist even though the same weights are used at each source position.

For a matrix comparison, define a fixed diagonal scale from the same case at
`t0=90,tau=0`: `D_ii=sqrt(Re C_90,ii(0))`, and
`Ctilde_t0=D^-1 C_t0 D^-1`. The following ranges are over the six paired
observations (three cases times two shifts) of
`||Ctilde_t0(tau)-Ctilde_90(tau)||_F / ||Ctilde_90(tau)||_F`.
They are **ranges of observed changes, not uncertainty intervals**.

| R/a | At tau/a=16 | At tau/a=32 |
|---:|---:|---:|
| 4 | 0.54% to 3.51% | 1.91% to 7.11% |
| 12 | 2.52% to 15.15% | 21.82% to 67.24% |
| 24 | 53.34% to 75.50% | 116.34% to 158.48% |
| 32 | 63.92% to 126.75% | 118.47% to 180.17% |

Small late-time denominator norms amplify relative changes. These numbers
do not measure a percentage change in the physical potential. They do not
distinguish ordinary configuration/position fluctuations from an ensemble
boundary effect.

Signed frozen projections show late negative real measurements and large
phases in some individual cases, especially at R24/R32. A negative
single-configuration correlator is allowed and is not a rejection criterion.
It does not establish a negative ensemble mean. All signs remain in the
plots and exports; no absolute values or positive-case filtering are used.
At R32 the imaginary part is at roundoff, but a negative real measurement has
principal phase +pi or -pi. Such a phase indicates sign, not a resolved
imaginary contribution. The two branch endpoints represent the same direction.

No source position is selected as better from this small pilot. In particular,
moving the source does not provide a demonstrated universal cure for the
late-time signal/window issues. It would take an ensemble position comparison
to assess the mean effect and its uncertainty.

## Statistical interpretation and normalization

The statistical cases number **three**, not nine positions or 99 time points.
All positions, axes, matrix entries, sources and times within a configuration
remain correlated. The technical cases were chosen for validation; no pilot
standard errors, significance tests, seven-point covariance fits, effective-
energy logarithms, physical GEVP levels or final potential are reported.

For the main projection plots, the common normalization is
`Re Q_cfg,t0(tau) / Re Q_cfg,90(0)`, separately for each frozen source and case.
The paired change is `(Re Q_t0-Re Q_90)/Re Q_90(0)`. This retains changes in
the equal-time source amplitude and never divides by a weak late-time Q.
A companion figure uses each position's own equal-time normalization to
separate shape from amplitude. Neither normalization creates new independent
measurements or an energy estimate. Raw real and imaginary values are retained.

## Plots

All eight figures have PNG and SVG versions. Blue/orange/purple denote source
positions 58/90/122 except where the legend explicitly names a control.
Each configuration is shown separately; these figures have no ensemble error bars.

1. [Signed pair/full projections](results/source_position_signed_correlators.png):
   three case rows and four radius columns. Common baseline equal-time scale;
   dashed pair and solid full curves; all positions and signed data retained.
2. [Own-position normalized projections](results/own_time_normalized_correlators.png):
   same measurements after normalization by each position's own equal-time Q.
   Compare with figure 1 to separate amplitude and time-shape changes.
3. [Phases and signs](results/source_position_phases.png): full basis, training
   separation 4a. Principal phases include negative-real sign flips at R32.
   The +/-0.05-radian lines are descriptive guides, not significance thresholds.
4. [Paired position changes](results/source_position_paired_changes.png): shifts
   58-minus-90 and 122-minus-90 for the same case and frozen source. Zero means
   no change; signed symlog axes retain positive and negative differences.
5. [Full complex matrix changes](results/matrix_position_changes.png): all
   sixteen ordered entries after one fixed diagonal scaling. Identical Y ranges
   aid radius comparisons; no matrix symmetrization or eigenmode deletion.
6. [All frozen basis/training controls](results/frozen_basis_training_controls.png):
   pair/triple/full and training separations 4/6/8 at tau=16 and 24. Each
   configuration/shift is retained. These are source-design controls, not a
   selection of a winning candidate from the pilot.
7. [Shifted projection control](results/shifted_projection_control.png): the
   n650 p10/p40 matrix comparison at t0=122. Differences are at floating-point
   scale and far smaller than the observed source-position changes.
8. [Measured launch timings](results/native_launch_timings.png): sum of MPI
   durations, distinguishing the audit, baseline, two shifts and p40 control.

## Tables

- `position_projections.csv`: all 10,692 primary position/source/time values,
  complex phases, real signs, two normalizations and paired changes.
- `full_anchor4_checkpoints.csv`: 216 display checkpoints, retaining every case.
- `matrix_position_diagnostics.csv`: 1,188 case/radius/position/time matrix
  norms, differences and complex/non-Hermitian descriptive diagnostics.
- `matrix_checkpoint_ranges.csv`: 24 range summaries; six paired changes each.
- `projection_control.csv`: all 132 axis-averaged R/time comparisons.
- `launch_timings.csv`, `launch_totals.csv`: eleven launches and three totals.
- `descriptive_flags.csv`: 108 source/position summaries of sign/phase counts.
  Counts are correlated configuration/time entries, not trials or probabilities.
- `VALIDATION.json`, `SUMMARY.json`: complete validation and scope.

## Measured timings

| Case | Audit | Scans | Sum of recorded MPI launches |
|---|---:|---:|---:|
| n170 | 4.76 min | two p10 | 12.06 min |
| n640 | 4.72 min | two p10 | 12.06 min |
| n650 | 5.23 min | fresh baseline, two p10 shifts, p40 control | 23.39 min |

Individual p10 scans took about 3.65 to 3.70 minutes; the p40 scan took
7.09 minutes. These sums exclude queue waiting, setup/shell overhead and gaps
between launches; the tasks need not have started together. Native receipts
preserve the actual Slurm task job IDs and array indices.

## Reproduction

Use the existing NumPy/Matplotlib analysis environment on Stromboli. No
installation or environment change is required by this package.

```bash
qx_analysis=laplace_static/analysis/qx_static_source_time_pilot
python3 "$qx_analysis/verify_files.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 "$qx_analysis/test_analysis.py"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 "$qx_analysis/reproduce.py" --verify
```

The integrity check is standard-library only and appropriate on SuperMUC.
Replay uses an isolated temporary extraction and removes it afterwards. It
reconstructs the native exports before recomputing every analysis table.
No optimizer is rerun on the shifted data and no plot is regenerated by replay.
The frozen-training check compares invariant normalization, leading-vector
overlap and residual bounds, avoiding platform-dependent eigenvector phases.
Only the original equal-time pivot diagnostic gets its established narrow
roundoff tolerance; primary native CSV exports remain exact.

For an intentional new empty output directory:

```bash
OPENBLAS_NUM_THREADS=1 python3 "$qx_analysis/run_analysis.py" \
    --output temp/qx_source_time_replay
python3 "$qx_analysis/render_plots.py" --output temp/qx_source_time_replay
```

The original 17 MB review ZIP and the 2.7 MB published primary cache must be
tracked with the code and results. Do not regenerate integrity receipts to
silence an unexplained mismatch. Exact source references are in
`inputs/SOURCES.json`. All package Python uses Python 3.6-compatible syntax.

## Stop point

Technical pilot validation, position-resolved review and reproducible export
are complete. The earlier full71 source/window analysis remains valid as its
own dataset; this pilot adds a temporal-position control without replacing it.

Any ensemble source-position extension, within-configuration position averaging,
further source designs, physical energy extraction or optimized flux profiles
remain future decisions. No next campaign is prepared or submitted here.
