# qX source-time pilot: design and input preflight

This package fixes the next small measurement design after the matrix
uncertainty study at commit `8b8b69863cd802fde5ea0d2b044aea5c98f91aa8`.
It contains an executable plan, frozen source vectors, existing reference data,
and a filesystem preflight. **The new native scan, its build, submission and
collection commands are not implemented by this design package.** Running its
commands cannot launch MPI or submit a Slurm job.

## What the existing results say

The full71 analysis used the same source time, t0/a=90, throughout. Its joint
anti-Hermitian component bootstrap did not reject the zero-mean hypothesis
(tail 0.3475 over the prescribed family), and additional held-out splits
changed the fitted energies little. Neither result proves that the remaining
time-window dependence is a boundary effect or that the measured matrix can
be replaced by its Hermitian part.

For the full source basis, training time 4a, original alternating held-out
split, the existing fitted energies are:

| R/a | 6--12 correlator window | 14--20 | 22--28 |
|---:|---:|---:|---:|
| 4 | 0.117029 +/- 0.000134 | 0.116571 +/- 0.000194 | 0.116665 +/- 0.000257 |
| 12 | 0.196742 +/- 0.000674 | 0.193479 +/- 0.001109 | 0.197949 +/- 0.002826 |
| 24 | 0.271310 +/- 0.002208 | 0.267016 +/- 0.004831 | 0.292235 +/- 0.018712 |
| 32 | 0.298572 +/- 0.003800 | 0.317732 +/- 0.015582 | 0.162659 +/- 0.038882 |

The paired late-minus-early R32 change is -0.13591 +/- 0.03956, even though
the late fit has an uncalibrated Gaussian tail near 0.526. R24's late window
has a phase reaching about 0.281 radians and a tail near 0.010. Source choice
also matters: in the early window the full-minus-pair change is
+0.002512 +/- 0.000288 at R24 and -0.009066 +/- 0.002453 at R32.
These are diagnostic, pointwise contrasts, not selected final potential values.
Exact source values are copied into `inputs/DIAGNOSTIC_REVIEW.json` and pinned
to the published fit and paired-contrast tables.

Thus the next measurement should vary the temporal position while keeping
the operator definitions and comparisons fixed. It should distinguish
source-position sensitivity from source retuning and weak late-time signal.

## Proposed technical pilot

| Parameter | Fixed choice |
|---|---|
| Configurations | n170, n640, n650; the already validated technical pilot set |
| Source times | 58,90,122; reuse the existing 90 measurements |
| Separations | Every integer tau/a=0..32 |
| Spatial targets | R/a=24,32 |
| Short-distance controls | R/a=4,12, retaining the existing kernel geometry |
| Cutoffs | 40,60,80,100 independently at source and sink; all 16 ordered pairs |
| Spatial averaging | All 64^3 sites on each axis and the combined axis average |
| Temporal links | Same one-step HYP parameters; piter=10 baseline |
| Fresh reference | Repeat t0=90,piter=10 on n650 |
| Projection control | t0=122,piter=40 on n650 |
| Proposed allocation | Three tasks, each eight nodes and 64 ranks, one thread/rank |

There are six shifted p10 scans, one fresh p10 reference and one p40 control:
**eight scan launches plus three saved-mode audits, eleven MPI launches**.
Each scan has 8,448 matrix records; 67,584 new records altogether. The existing
baseline at t0=90 is reused on all three configurations. A fresh reference
regression on n650 must succeed before interpreting source-time differences.
No new eigenpairs are generated and no full71 source-time array is authorized
or implemented by this package.

The most distant sink and nearest boundary distance at tau=32 are:

| t0/a | Interval | Distance to left boundary | Distance to right boundary |
|---:|---:|---:|---:|
| 58 | 58..90 | 58 | 101 |
| 90 | 90..122 | 90 | 69 |
| 122 | 122..154 | 122 | 37 |

Coordinates are 0..191. None of these paths wraps in time. This is a
source-position comparison, **not** a pair of exactly reflected intervals;
in particular fixed sources 58 and 122 are not exact time-reflection partners
for all tau. The spatial half box is R/a=32, while half the temporal extent is
96. Temporal measurements here still stop at tau=32.

## Keep source training separate from the pilot

`inputs/FIXED_SOURCES.json` contains 36 source vectors: three previously
defined bases, training times 4,6,8, and four radii. They are trained only on
the **68 existing configurations outside n170,n640,n650**, at reference source
time 90. Metric normalization and the source-selection algorithm are unchanged
from the published matrix analysis. The primary controls are pair/full bases
at training 4a; the other fixed choices remain secondary diagnostics.

Exactly the same coefficients must be used at source times 58,90,122. Do not
optimize separately on each of the three pilot samples or on each new time
window. This makes a temporal-position change interpretable without an
additional change in source training. The frozen-source pilot checks are
conditional on that separate training sample; they are not final GEVP levels
or a full propagation of its training uncertainty.

`inputs/reference_t90_axis_average.csv` contains 6,336 **already measured**
baseline entries for the three pilot configurations. They come directly from
the published primary cache, with sink rows and source columns. The CSV
retains explicit source/sink labels, complex values and the count 3*64^3.
It supplies the combined-axis reference; future native validation must also
retain the per-axis and rank-block checks rather than treating this CSV as a
replacement for them.

## Statistical scope

The three configurations are for technical validation, orientation and
projection controls, input audit and timing. A seven-time covariance from
three configurations has rank at most two. **Do not run the full71 correlated
fit machinery on this three-configuration pilot or manufacture a nonsingular
covariance by counting axes, spatial blocks or source times as configurations.**

For a later ensemble analysis:

1. First retain position-resolved matrices and correlators, and compare
   shifted-minus-reference quantities on identical configurations.
2. Apply the same frozen or correctly cross-fitted source definition to every
   position in each comparison. A future production jackknife must carry all
   source times from a configuration together and redo any source training.
3. Compare signed signals, phases, effective energies and consistent fit
   windows. Preserve negative or unresolved measurements and their validity
   masks.
4. Only then assess an equal-weight source-time average within each
   configuration. Three source times do not triple the ensemble size, and an
   apparent statistical gain does not establish absence of boundary bias.

The target for a later physics measurement remains all 71 configurations.
If the shifted-time design is retained, the six validated shifted pilot scans
can be reused; the other 68 configurations would require 136 shifted p10
scans. That production step is a separate measurement decision.

## Input preflight

All 100 modes on all required slices must be audited before the future scan.
The two shifted sources require 66 unique slices on n170 and n640. The fresh
reference on n650 fills the gap between them and requires 97 slices there.
The total planned inventory is 229 slices, or 22,900 saved eigenvectors.

`preflight.py` uses the already published RUNS registry, checks its full71
labels, verifies original completion markers and the recovered n640 eigenvalue
checksum, then checks nonempty gauge, eigenvalue and required mode files.
It records file size, resolved path and nanosecond modification time without
reading or hashing the large mode payloads. The eigenvalue files and admission
markers are small and are hashed. Its report is saved under this repository's
ignored `work/` directory.

A successful preflight means **the input inventory is available**. It does not
check HDF5 attributes, eigenvector residuals, normalization or Gram matrices;
those require the later native saved-mode audit. The preflight does not load
the gauge field, execute the old binary, compile code or contact Slurm.

On Stromboli, validate the installed design and regenerate its existing-data
references without writing them:

```bash
qx_plan=laplace_static/validation/qx_static_source_time
python3 "$qx_plan/verify_files.py"
python3 "$qx_plan/test_plan.py"
OPENBLAS_NUM_THREADS=1 python3 "$qx_plan/generate_references.py" --verify
```

Only reference regeneration needs NumPy. Integrity checks, geometry tests and
the SuperMUC preflight use the Python standard library and Python 3.6 syntax.

After synchronization, run on SuperMUC:

```bash
qx_plan=laplace_static/validation/qx_static_source_time
python3 "$qx_plan/verify_files.py"
python3 "$qx_plan/test_plan.py"
python3 "$qx_plan/preflight.py" --output "$qx_plan/work/PREFLIGHT.json"
```

Missing/empty files, duplicate registry entries and changed admission receipts
are failures. A blocked case remains explicitly blocked. Physical inputs are
not changed, moved, generated or deleted.

## Native implementation checklist after preflight

The existing driver is intentionally left unchanged. It hardcodes t0=90 in
the generated geometry, emitted matrix/direct/block metadata, begin header and
SU(3) window. A later isolated extension must change these consistently; merely
changing the start used in the temporal product would mislabel the measurement.

The extension needs:

- Explicit source-time selection, endpoint caches and non-wrapping temporal
  products; every output and receipt includes actual t0 and t1.
- A strict native 100-mode audit over each case's declared 66/97 slices.
- All previous independent rectangular contractions, diagonal-kernel controls,
  complete periodic payload shifts, axis/rank reductions, equal-time Gram
  checks and R32 imaginary cancellation.
- SU(3) diagnostics for the actual shifted transporter windows, using the
  existing thresholds. A failing p10 projection must not be relabelled as a
  boundary effect or admitted automatically by switching piter.
- Fresh n650 t0=90 regression against published references, plus matched
  p10/p40 comparison at t0=122. Independent direct tests alone are insufficient
  to establish that source-time indexing is correct.
- Distinct immutable attempts and receipts, strict missing/duplicate output
  detection, reference replay, and prevention of duplicate submissions.

## Timing basis

The published full71 native p10 scans took **213.1--231.5 seconds**, with mean
219.7 seconds, for the same matrix geometry at source time 90. This is a useful
scan benchmark, not a runtime measurement of the extension. The new audit
reads two to almost three times as many slices per case; the n650 task has four
scans. A provisional planning allowance is **20--45 minutes after allocation**
for the slowest task, with a one-hour limit to be checked against actual pilot
timings. Queue delay and staggered task starts are additional. No scheduler
start time or new benchmark is claimed here.

## Method context

Open temporal boundaries make temporal position a relevant control. Luescher
and Schaefer discuss the open-boundary formulation in
[Lattice QCD without topology barriers](https://arxiv.org/abs/1105.4749).
Section 4.1 of
[The determination of potential scales in 2+1 flavor QCD](https://inspirehep.net/files/9bc1018085876caa2ba424f0ac2c5ce8)
explicitly examines Wilson-loop dependence on source position and the temporal
region retained for averaging. This motivates a direct position check; that
paper's numerical boundary cuts and its symmetrized Wilson-loop GEVP procedure
are not automatically prescriptions for this Laplacian-source qX observable.

This design does not diagnose the existing late-R32 change as a boundary
effect, identify a final plateau, or request further clover production.
