# qX full-ensemble mixed-source matrix

Measure the existing four-source matrix on the complete 71-configuration qX ensemble, using the unchanged validated matrix pilot binary and output checker. This directory adds preparation, scheduling, exact pilot reuse, resumable case execution, and full-cohort collection. It does not change the native contractions, cutoffs, projection iterations, saved modes or numerical tolerances.

The workflow is based on commit 5dd70f242ab81ca825b8c7c9cc221725fe9438a8. The existing three-configuration matrix pilot is job 5961383.

## Scientific purpose

The radius-by-radius analysis shows that a locally flat effective-energy window does not by itself establish a common ground-state energy. At R/a=4, source cutoffs agree closely, while at R/a=12 the paired Nv=40 minus Nv=100 difference in window 14–19 is approximately -0.0016367 +/- 0.0004814. Near R/a=24, late-time signal is much stronger for Nv=40 than Nv=100. These are source-overlap, time-window and signal questions.

The mixed entries supply information absent from the existing diagonal correlators. They allow later comparisons of two-, three- and four-operator bases, fixed optimized projections, and statistically assessed matrix properties. This measurement campaign does not automatically apply a Hermitian projection, accept a GEVP level, choose a potential or determine a string tension.

## Fixed measurement plan

| Parameter | Value |
|---|---|
| Configuration labels | 100 through 800 in steps of 10, all 71 including n640 |
| Validated pilot reuse | n170, n640, n650 |
| New configurations | Remaining 68 |
| Spatial separations R/a | 4, 12, 24, 32 |
| Correlator separations tau/a | Every integer 0 through 32 |
| Absolute source time | 90 |
| Independent source and sink cutoffs | 40, 60, 80, 100; all 16 ordered pairs |
| Primary SU(3) projection iterations | piter=10 |
| Saved eigenvectors | Existing lowest 100 at absolute times 90 through 122 |
| Allocation per array task | 8 nodes, 64 MPI ranks, one thread per rank |
| Array organization | 34 tasks, two new configurations per task, at most four concurrent |
| Native work per new configuration | Strict 33-slice saved-mode audit, then one p10 matrix scan |

R/a=32 is the spatial half-extent on 192 x 64^3, not the temporal half-extent. The latter is 96. Times through 32 are the already validated measurement range.

The four radii retain the pilot geometry: R4 is a strong-signal control; R12 tests the observed cutoff dependence; R24 addresses the sensitive source region; R32 covers the spatial half-box. Adding neighboring radii or changing to exact smooth eigenvalue profiles requires corresponding new pilot measurements and a new native validation stage.

The production combined tables contain 599,808 p10 matrix rows and 872,448 rank-block rows. The n650 p40 control remains in the archived pilot evidence and is replayed, but does not enter the p10 production table a second time.

## Validation and provenance

- Check all 53 pinned native source/reference/archive dependencies and the new orchestration files.
- Require the original build manifest and the same binary/library checksums used by the pilot.
- Require exact pilot geometry, three completed cases, correct p10/p40 controls, unchanged raw stdout/stderr, and preparation/case/launch receipts matching the published job5961383 archive.
- Reparse the pilot outputs and require exact agreement with their CSV exports before preparing production.
- Preserve the native saved-mode residual/norm bounds, SU(3) checks, direct rectangular contractions, legacy diagonal regression, periodic shifts and equal-time Gram checks.
- Retain complex ordered entries at nonzero time. No per-configuration Hermitian symmetrization is introduced.
- Bind every new attempt to its immutable prepared plan. A failed attempt never acquires a successful latest marker.
- Verify all 71 configurations before full collection, reparse each admitted attempt, and stream the combined tables without retaining the entire production dataset in memory.

The new tests cover orchestration and provenance, using explicit mocks for MPI and physical storage. They also reparse the real published pilot archive with the unchanged native output checker and require exact CSV agreement. That offline replay substitutes only the unavailable physical-build/storage lookup; it does not execute native binaries. The native numerical implementation is not modified; its previous native validation and current source/build pins remain the numerical foundation. The actual SuperMUC preparation and subsequent native runs still have to pass their physical-input and measurement checks.

## Installation and source checks on Stromboli

From the repository root, after extracting the workflow ZIP:

    qx_matrix=laplace_static/validation/qx_static_matrix_full71
    python3 "$qx_matrix/verify_files.py"
    python3 "$qx_matrix/test_workflow.py"
    python3 "$qx_matrix/campaign.py" check-inputs
    python3 "$qx_matrix/campaign.py" plan

These checks submit no jobs. They require only Python's standard library and Bash for the shell syntax check; they do not need physical gauges, eigenvectors or the SuperMUC binaries. Source syntax is compatible with Python 3.6.

Commit and push the new directory from Stromboli, then transfer a Git bundle to SuperMUC using base commit 5dd70f242ab81ca825b8c7c9cc221725fe9438a8.

## Prepare on SuperMUC

Use the original retained native build and the original pilot campaign. No rebuild is needed when their pinned files remain intact.

    cd /hppfs/work/pn29se/di24gus/Masterarbeit
    qx_matrix=laplace_static/validation/qx_static_matrix_full71
    python3 "$qx_matrix/campaign.py" check-inputs
    python3 "$qx_matrix/campaign.py" prepare --name full71 --pilot pilot

Preparation validates the physical records for all configurations and strictly admits the three pilot cases. It creates a new immutable full71 plan and does not call sbatch. A changed native build or incomplete pilot stops preparation. Retain the original pilot directories, gauge/eigenvector storage and n640 recovery paths.

## Submit and monitor

This command starts the 34-task production array:

    python3 "$qx_matrix/campaign.py" submit \
      --name full71 --start 0 --count 34 --concurrent 4

    squeue -u "$USER" -o '%.18i %.20j %.2t %.10M %.6D %R'
    python3 "$qx_matrix/campaign.py" status --name full71

Prior live or ambiguous submissions block another submission. Completed tasks are skipped. A task interrupted after its first configuration resumes by checking and reusing that completed case.

The pilot recorded approximately 4 minutes for an audit plus 3.6 minutes for a p10 scan. A two-configuration task therefore has about 15 minutes of native work. Allow roughly 2–3 hours for the remaining 68 configurations with four tasks continuously active, plus scheduler delays, filesystem variation and Python validation overhead. The one-hour task limit leaves headroom; this is not a guaranteed completion time.

Each running case has an exclusive RUNNING.json receipt. A hard-killed task may leave that receipt behind. Inspect the recorded job and confirm it has ended before clearing that one stale receipt. There is no automatic expiry or broad deletion of campaign outputs.

## Collect

Once the jobs finish:

    python3 "$qx_matrix/campaign.py" collect --name full71

Successful collection requires exactly all 71 configurations and produces:

- campaigns/full71/matrix_correlators.csv: sorted-by-configuration p10 production entries.
- campaigns/full71/rank_blocks.csv: p10 rank-reduction evidence.
- campaigns/full71/SUMMARY.json: full cohort, reuse, checksums, native checks and launch timing.
- campaigns/full71/ARCHIVE_FILES.json: byte counts and SHA-256 inventory.
- qx_static_matrix_full71_review.zip in the repository root.

The archive includes the selected pilot evidence, new attempts, preparation/submission receipts, and the native source/headers/text build provenance listed by the original build manifest. Native executable/library bytes, the already published 37 MB diagonal reference, and the original published pilot ZIP are not duplicated; their hashes and reference paths are retained, with omitted native inputs listed explicitly. The reused pilot's actual attempt files are included directly. Incomplete temporary files are listed explicitly as omitted. A final ZIP CRC check precedes publication of the archive.

Transfer the review archive for offline analysis. Remove transfer ZIPs/bundles after the receiving side verifies them; retain the native build, campaigns, saved modes and recovery evidence.

## Analysis after collection

Use configurations as the independent statistical units, retaining cross-source/time covariance. Its empirical rank is at most 70: the complete source/time feature vector must not be inverted as one unrestricted covariance matrix. Diagnose appropriately sized fit blocks and projections instead. Compare the full ensemble's finite-time Hermiticity/phase behavior, basis conditioning and smaller bases before interpreting variational spectra. Train fixed source projections consistently and assess fitting/selection uncertainty. Compare resulting energy candidates with the existing diagonal window scans, especially at R12 and R24. Exact smooth-profile measurements and optimized clover profiles remain subsequent stages.
