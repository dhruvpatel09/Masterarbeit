# qX mixed-source static-correlator pilot

This workflow measures the complete source/sink matrix for four Laplacian
mode cutoffs. It extends the published diagonal measurements at commit
`43fbfc21413a81c5536285f53e96de3de67f13d1` without changing their parent code.
The purpose is to validate the contractions needed for a later variational
analysis of source overlap and excited-state contamination.

## Definition and orientation

Let

\[
K_N(x,y;t)=\sum_{n=0}^{N-1}v_n(x,t)v_n(y,t)^\dagger,
\qquad
T_x(t_0,t_1)=\prod_{t=t_0}^{t_1-1}U_0^{\rm HYP}(x,t).
\]

For source cutoff \(N_s\), sink cutoff \(N_k\), and \(y=x+R\hat e_a\),
the measured complex number is

\[
C_{N_s,N_k}(R,\tau)=\frac{1}{N_{\rm src}}
\sum_x\operatorname{Tr}\left[
T_y K_{N_k}(y,x;t_1)T_x^\dagger K_{N_s}(x,y;t_0)
\right],\qquad t_1=t_0+\tau.
\]

This is also the independent rectangular sum
\(\sum_{i<N_s,j<N_k} A_y(i,j)A_x(i,j)^*\), where
\(A_x(i,j)=v_i(x,t_0)^\dagger T_xv_j(x,t_1)\).
The diagonal \(N_s=N_k\) has exactly the existing contraction and normalization.

The C helper returns `out[source_slot][sink_slot]`. The CSV records explicit
`Nv_src` and `Nv_sink` labels. A conventional later GEVP matrix with **sink
rows and source columns** must be assembled as
`C[sink_slot, source_slot] = record(Nv_src, Nv_sink)`.
This indexing convention matters for the off-diagonal entries.

All 16 ordered pairs and both real and imaginary parts are retained. At nonzero
time separation a fixed-source, single-configuration matrix is not assumed
Hermitian. The workflow neither symmetrizes it nor reconstructs a missing
triangle. At zero separation it checks the Hermitian positive-semidefinite
Gram identity within roundoff. This check is not a condition-number cut for a
future GEVP. Swapping the spatial endpoints conjugates each ordered-pair
correlator; the global imaginary part therefore cancels at the half-box
separation `R=32`.

Cross-source correlators are new measurements. They cannot be reconstructed
from the covariance of the existing diagonal correlators.

## Technical pilot

| Setting | Value |
|---|---|
| Lattice and boundary conditions | `192 x 64^3`; open time, periodic space |
| Configurations | `170, 640, 650`, including recovered n640 |
| Source time | `t0=90` |
| Time separations | Every integer `tau=0..32`; 33 saved-mode slices |
| Spatial separations | `R=4,12,24,32` |
| Source and sink cutoffs | Independently `40,60,80,100`; all 16 combinations |
| Spatial averaging | All `64^3` sites for each of three axes, plus their average |
| Temporal HYP | One step `(1,1,0.5)`; baseline `piter=10` |
| Projection control | Full `piter=40` matrix on n650 |
| MPI allocation | 64 ranks, grid `1x4x4x4`, eight nodes, one thread per rank |
| Submission | Three array tasks, at most three concurrent, one-hour limit each |

Each task first audits all 100 saved modes on all 33 required times, then
measures its p10 matrix. The n650 task additionally measures the p40 matrix.
There are seven MPI launches altogether: three audits and four scans.
No eigenpairs are generated.

This is a validation pilot. The intended physics production cohort remains
**all 71 configurations, n100 through n800 in steps of 10**. A full-cohort
submission should follow review of the pilot's complex matrices, controls and
actual timings. This package does not implement a production-stage command.

## Validation and output inventory

The original 100-mode HDF5 buffers and complete MPI endpoint payload are reused.
The mixed kernel forms four prefix sums at each endpoint once, then contracts
all source/sink combinations. No partial-mode MPI buffer is transmitted.

Per scan, collection independently checks:

| Output | Expected coverage |
|---|---:|
| `MATRIX_RT` | 8,448 complex records: 16 pairs x 33 times x 4 radii x 4 axis choices |
| Published diagonal matches | 2,112 records, including every time and every cutoff |
| `MATRIX_DIRECT` | 1,536 independent quadratic checks at `tau=0,4,16,32`, all radii and axes, local sites 0 and 4095 |
| `MATRIX_LEGACY` | 384 comparisons to the unchanged diagonal kernel |
| `MATRIX_BLOCK` | 12,288 complex records at `tau=0,32`, `R=4,32`, all 16 pairs, three axes and all 64 ranks |
| Full-payload roundtrips | Three, each over all 64 periodic spatial cells |
| Equal-time Gram checks | 16: four radii x four axis choices |
| Central-link diagnostics | 32 times and their combined time window |

The pilot contains 33,792 complex matrix rows and 49,152 rank-block rows across
the four scans. Block coordinates are checked for completeness and consistent
rank mapping; their averages must reproduce the global sums. Blocks diagnose
spatial averaging and do not become independent ensemble configurations.

Saved-mode relative residual, normalization and Gram thresholds remain
`1e-9`. The central-link unitarity/determinant threshold remains `1e-12`.
The replay recomputes direct-test errors and tolerance bounds from the logged
complex numbers. Published p10 diagonal comparisons use relative tolerance
`5e-12`, absolute tolerance `1e-20`; the p40 diagonal and full-matrix projection
controls use relative tolerance `1e-10`, absolute tolerance `1e-20`.
No tolerance is enlarged automatically after a failure.

`INPUTS.json` pins the new sources and their published dependencies, including
the existing full71 temporal-extension raw CSV. Native build provenance also
records actual libraries, headers, generated helpers and compiler/module logs.
Preparation binds the binaries and physical inputs to a checksum of the plan.
Audit, stdout, stderr and launch inventories are replayed before a successful
case marker is written. Failed attempts are preserved; a retry gets a new
attempt directory. Active or ambiguous submission receipts prevent a duplicate
submission. An incomplete collection creates an explicitly incomplete review
archive and exits unsuccessfully.

## Local tests and scope

Run `python3 laplace_static/validation/qx_static_matrix/test_workflow.py`.
It requires only Python's standard library, Bash and GCC; Python code is
compatible with 3.6. Tests exercise:

- Full rectangular contractions at both 4 and 100 saved modes, against
  independent sums, ordered-pair one-hot controls, non-diagonal SU(3) gauge
  changes, endpoint conjugation and the equal-time Gram identity.
- The existing full-payload periodic-shift test and native driver syntax
  against the captured interfaces; generated I/O, shift and audit helpers
  must match their parent byte for byte.
- Complete synthetic complex output, including asymmetric finite-time matrices
  and negative real values; missing/duplicate records, wrong reductions,
  inflated tolerances and altered references must fail replay.
- Real campaign preparation, mocked MPI launch receipts, failed attempts,
  retries, partial/final review exports, changed inputs and changed logs.

Synthetic fixtures are labelled as such and are never presented as physical
measurements. These tests do not replace compilation against the actual MPI,
HDF5 and solver libraries or the subsequent multi-node pilot.

## Install and publish on Stromboli

Place `qx_static_matrix_pilot_workflow.zip` in the repository root, then:

```bash
(
    set -e
    cd /home/m2130292/Masterarbeit
    git pull --ff-only origin main
    unzip -n qx_static_matrix_pilot_workflow.zip

    qx_matrix=laplace_static/validation/qx_static_matrix
    python3 "$qx_matrix/test_workflow.py"
    python3 "$qx_matrix/campaign.py" check-inputs

    git add -- "$qx_matrix"
    git commit -m "Add mixed-source qX static-correlator validation pilot"
    git push origin main

    git bundle create qx_static_matrix_sync.bundle main \
      ^43fbfc21413a81c5536285f53e96de3de67f13d1
    scp qx_static_matrix_sync.bundle \
      di24gus@skx.supermuc.lrz.de:/hppfs/work/pn29se/di24gus/Masterarbeit/
)
```

## Import, build and submit on SuperMUC

The expected starting commit is `43fbfc21413a81c5536285f53e96de3de67f13d1`.
The compare-and-swap update below follows a successful Stromboli push.

```bash
(
    set -e
    cd /hppfs/work/pn29se/di24gus/Masterarbeit

    git bundle verify qx_static_matrix_sync.bundle
    git fetch ./qx_static_matrix_sync.bundle main
    git merge --ff-only FETCH_HEAD

    qx_matrix=laplace_static/validation/qx_static_matrix
    python3 "$qx_matrix/campaign.py" check-inputs
    git update-ref refs/remotes/origin/main "$(git rev-parse HEAD)" \
      43fbfc21413a81c5536285f53e96de3de67f13d1
    rm -- qx_static_matrix_sync.bundle

    bash "$qx_matrix/build.sh"
    python3 "$qx_matrix/campaign.py" prepare --name pilot
    python3 "$qx_matrix/campaign.py" submit --name pilot --concurrent 3
    squeue -u "$USER" -o '%.18i %.20j %.2t %.10M %.6D %R'
)
```

Do not rerun `build.sh` over an existing build: it deliberately stops. A source
or dependency change needs an inspected new build and preparation, not an
overwrite of a running campaign.

The previous temporal-extension pilot's slowest audit plus p10 scan took
581.17 seconds (9.7 minutes). For this new kernel, a provisional budget is
**15–30 minutes once all three pilot tasks have allocations**. The extra n650
projection scan is included in that estimate; it is not a new benchmark or a
scheduler prediction. Queue delays, delayed allocation of individual tasks and
filesystem variability add to elapsed time. Each task has a one-hour limit.
`LAUNCHES.json` records actual seconds separately for every audit and scan.

## Collect and review

After the array has finished:

```bash
cd /hppfs/work/pn29se/di24gus/Masterarbeit
python3 laplace_static/validation/qx_static_matrix/campaign.py collect --name pilot
ls -lh qx_static_matrix_pilot_review.zip
```

Upload `qx_static_matrix_pilot_review.zip` for replay and inspection. It contains
the plan, build provenance, launch timings, all stdout/stderr logs, full complex
matrix CSVs and rank-block CSVs. Physical gauge files and saved eigenpairs are
not copied into it. Retain the original campaign and mode storage. Remove
transfer ZIPs and bundles after the destination has verified them.

The next decision is whether the validated basis has adequate numerical rank
and a reproducible signal for a full71 mixed-matrix measurement. GEVP choices,
covariance handling and source/time-window comparisons belong to that analysis;
this pilot alone does not establish a ground-state potential or string tension.
