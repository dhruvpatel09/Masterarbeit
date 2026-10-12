# Projection sensitivity across qX configurations

This additive package reuses the profile calculation validated in jobs 5903817,
5904081 and 5915129. It reads existing gauge fields and eigenvectors for n110,
n120, n130, n140 and n150; it creates no new eigenvectors and changes no existing
solver library or measurement source. `adaptation.diff` records changes from the
hash-pinned n100 source.

## Measurements

Fifteen array tasks cover five configurations times tau=2,4,6, at midpoint98.
Each task measures piter4 and piter40 independently from the original field:
one temporal HYP step, alpha=(1,1,0.5), R=1,2,4, the existing 332 probe offsets
and axis combinations. The clover probe uses the original field and the existing
RAW convention (no traceless subtraction or extra factor1/2). Spatial APE20,
alpha0.5 defines the existing eigenvector operator; it is separate from temporal HYP.

Task0 also repeats n100 tau4 at both projection iteration counts and compares
all 664 raw profile records to the actual export from job5915129. Complete
collection success requires this regression. The reference provenance and
original source hashes are included.

Each task uses the previously validated observable grid1x4x4x4 (64 ranks,
8 ranks per node, one thread per rank), not the eigenvector generation grid.
`solver_mpiInit` provides MPI_THREAD_FUNNELED. Two tasks may run concurrently,
each on eight nodes in the test partition. The requested time is30 minutes;
previous n100 paired measurements took approximately seven minutes, and task0
has an additional pair. New configurations and shared I/O can change runtime.

## Quality and scope

The saved-mode audit5917325 passed on nine selected slices in each configuration.
All eigenvalue time coverage was checked. Tau4 endpoints96/100 were included in
the residual/Gram audit; tau2/6 endpoints95/97/99/101 were not. Their inventory is
checked before submission and their BC, spatial APE, configuration, time, Nv and
run-name metadata are checked when read. The new measurement is finite/structural
validation, not a new independent eigenmode residual audit.

SU(3) diagnostics cover every physical temporal link used by each central
transporter, t0 through t1-1. They report maximum and RMS ||U-dagger U-I||_F,
|det U-1|, and counts exceeding1e-2, both by time and by window. They do not
include the open-boundary t0/t190 tails outside these transporter windows.
piter40 is required to meet a1e-12 final-field control tolerance. piter4 has no
automatic physics acceptance threshold. Its defects and profile differences
must be assessed together across configurations. SU(3) membership alone does
not establish projection convergence or absence of observable bias.

The inherited checks enforce profile coverage, translations, source counts,
MPI reductions, direct samples, halo handling, transporter factorization,
unchanged spatial/boundary links and unchanged original clover probes.

Exports retain raw complex L, LS and real S quantities. `C=Re(LS)/Re(L)-S` and
its relative differences are per-configuration diagnostics, evaluated using
50-digit Decimal arithmetic on the printed inputs. They are not ensemble
estimates. Form ensemble raw averages before division/subtraction and estimate
correlated uncertainties; sites, axes and tau values are not independent
configurations. Relative diagnostics are undefined when the reference
denominator is zero. No group averaging or statistical significance is claimed.

## Run

From `/hppfs/work/pn29se/di24gus/Masterarbeit`, after unpacking:

```bash
bash laplace_static/validation/hyp_projection_configs_qx/build.sh
python3 laplace_static/validation/hyp_projection_configs_qx/prepare.py
sbatch laplace_static/validation/hyp_projection_configs_qx/submit.sh
```

The build/preparation checks stop on changed provenance, missing inputs or
existing preparation/build paths. A checksum failure should be investigated;
do not regenerate hashes merely to bypass it. No existing results are overwritten.

After all array tasks leave the queue, replace ARRAY_JOB_ID below with the
submission's numeric parent ID:

```bash
python3 laplace_static/validation/hyp_projection_configs_qx/collect.py ARRAY_JOB_ID
```

The collector rechecks all outputs, requires all15 cases and the n100 regression,
and creates `hyp_projection_configs_review_ARRAY_JOB_ID.zip` inside Masterarbeit.
It includes small logs, sources, manifests and CSVs, without HDF5 data, binaries
or libraries. On missing/failed tasks it still creates a partial review ZIP and
exits nonzero; partial output is not labelled a complete success. Download and
share the ZIP before deleting that exact archive.

Local validation of this package checks C syntax using exact recovered qcd/reader
headers and MPI/HDF5/solver declaration stubs, the SU(3) metric against
NumPy, checker/export logic with positive and negative fixtures, checksum and
ZIP integrity, and shell syntax. Compilation/linking and the fresh n100 physics
regression must run on SuperMUC; they have not been executed locally.
