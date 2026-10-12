# Physical temporal-link projection by time

This adds time-resolved statistics to the successful full-field test in
`hyp_projection_qx` (job 5914979). It locates the measured piter=4 SU(3) tails
and checks the link intervals traversed by the midpoint-98 observable scans.
The previous observable comparison is included in `analysis/REPORT.md`.

## Run

From the Masterarbeit root, after extracting this package:

```bash
bash laplace_static/validation/hyp_projection_timeslices_qx/build.sh
sbatch laplace_static/validation/hyp_projection_timeslices_qx/submit.sh
```

The job uses 8 nodes, 384 MPI ranks, grid 6x4x4x4, and one compiled OpenMP
thread per rank. Initialization remains `solver_mpiInit` (FUNNELED).
It uses the existing qX build of libsolver and reads qX_p2n100. No eigenvectors
are generated or read by this diagnostic.

Build and submission verify the new package, relevant solver source and
geometry hashes, and the original test source. A changed input requires
inspection before updating its checksum. Existing build/output directories
are preserved by failure guards.

## Measurement

The C program keeps the existing round-trip, zero-alpha, spatial-link,
absent-boundary-link, independent piter=4 repeat, and piter=40 SU(3) controls.
It uses the same one-step temporal HYP calls and alpha=(1,1,0.5).
Added arrays accumulate the already computed final-link defects by global
time; MPI reductions use maxima, sums of squares, and integer counts.

For ORIGINAL, PITER4 and PITER40, it reports each physical temporal-link
slice t=0..190: 262144 links per slice, maximum entry and Frobenius unitarity
defects, RMS Frobenius defect, maximum and RMS determinant defects, and
exceedance counts at 1e-12, 1e-8, 1e-4 and 1e-2. It also reports each slice's
maximum real/imaginary component difference between piter=4 and 40.
Absent outgoing temporal links at t=191 are excluded.

The Python checker requires all 573 SU(3) records and all 191 difference
records and reconstructs the separate global maxima, RMS and counts. It
exports into `results_JOBID/`:

- `diagnostics.out`: application output;
- `su3_by_time.csv`: 573 SU(3) rows;
- `piter_difference_by_time.csv`: 191 difference rows;
- `SUMMARY.json` and `SUMMARY.md`: full-field and measured-window summaries;
- `manifest.sha256`: source/binary provenance.

Relevant physical link windows are t=97..98 (tau=2), t=96..99 (tau=4), and
t=95..100 (tau=6). The upper endpoint itself has no additional outgoing
link in the transporter.

No piter=4 closeness threshold is imposed. This measures where the existing
tails occur; it does not assume they are boundary-localized. It checks final
temporal links, rather than intermediate decorated HYP links, and does not
establish ensemble accuracy.

## Local validation

C syntax was checked against the reviewed qcd/solver headers and an MPI
declaration shim; shell and Python syntax also passed. Per-time matrix
defects, exceedance counts and component differences were checked against
independent NumPy calculations, including exclusion of a deliberately invalid
upper-boundary link. Output inventory, global reconstruction and malformed
time/count rejection were checked. These local checks do not exercise
distributed MPI or substitute for execution on SuperMUC.

## Independent scaling task

Leave eigenvector benchmark 5915130 queued. Once it completes, collect its
application output, `TIMING.json` and `COMPLETED.json`. The 8/16/32-node
comparison will use application-plus-writing elapsed time and allocated
physical-core-hours, with the previously noted timing-scope difference for
the 16-node pilot. The observable timing is a separate workload.
