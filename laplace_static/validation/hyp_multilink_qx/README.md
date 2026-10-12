# qX four-link temporal HYP and clover integration

Additive extension of the qX one-link driver validated in job 5903317.
Build from Masterarbeit with `bash laplace_static/validation/hyp_multilink_qx/build.sh`.
Submit with `sbatch laplace_static/validation/hyp_multilink_qx/submit.sh`.
The executable and reader object live in this directory's ignored `build/`.
The existing solver library for grid 1x4x4x4 is reused unchanged.

## Fixed validation case

- qX_p2n100, lattice 192x64^3, 64 MPI ranks, eight nodes, one thread/rank.
- solver_mpiInit, MPI_THREAD_FUNNELED, -fopenmp.
- Existing Nv=100 eigenvectors at t=96 and 100; uniform mode weights.
- P(x) = U0(96,x) U0(97,x) U0(98,x) U0(99,x), in this order.
- Products are constructed locally and exchanged into positive spatial halos.
- No temporal wrapping. Even temporal extent and interior endpoints required.
- R=1 pairs: (14,0,0)->(15,0,0) and (15,0,0)->(16,0,0).
- Probe at t=98, spatial position (14,1,2) or (15,1,2), respectively.
  These are integer-site offsets (0,1,2) from the adjointed endpoint;
  they are not the half-integer spatial midpoint of the R=1 pair.
- Clover uses original gauge links: F=(Q-Qdag)/8, no trace subtraction,
  S_RAW=sum over six planes Tr(Fdag F), no factor 1/2.
- Also prints the original-field spatial average of each plane at t=98.
- Clover storage is released before HYP allocation.

## Controls and interpretation

Five separate MPI runs use the same driver:

1. baseline: mode 3 copies the original field, bypassing solver conversion/HYP.
2. off: mode 0 imports/exports through solver, with no smearing.
3. zero: mode 1 uses alpha=(0,1,0.5), piter=40.
4. hyp40: mode 2, one temporal HYP step, alpha=(1,1,0.5), piter=40.
5. hyp80: same with piter=80.

Here `baseline` denotes this driver's original-field control, not an independent
implementation of the multi-link observable. Direct versus factorized
contractions provide an algebraic cross-check. The split-product check also
checks temporal multiplication grouping; neither proves the HYP algorithm
independently. All-rank field isolation checks run before links are replaced
by transporters. The transport buffer is never subsequently smeared or used
for clover. The probe is computed before HYP from the original field.

Comparison requires baseline/off/zero agreement, hyp40/hyp80 stability,
unchanged probes, complete output counts, finite data, local/cross-boundary
factorization PASS, and valid L*S products. It compares printed physical values
with absolute tolerance 1e-18 plus relative tolerance 1e-11. Factorization
requires relative difference <1e-10; split products require absolute 1e-12.

The output contains raw complex L and LS for two pairs, and real S_RAW.
There is no source average or ensemble connected subtraction in this test.
Do not interpret LS/L-S on one pair/configuration as a flux profile.
The larger production source/axis/R/temporal scans follow this validation.

## Local verification and remote gate

The complete C file passed GCC C11 -Wall -Wextra -Werror syntax checking
using the uploaded project headers and minimal MPI/HDF5 declaration stubs.
That is not a link or MPI runtime test. The actual ordered-product helper
also passed an independent explicit four-factor contraction check on
noncommuting complex matrices (max difference 1.24e-16), and the clover
helper passed a normalization/trace-inclusion check. The comparator was
checked with a synthetic fixture and rejected bad LS, wrong time endpoints,
and failed factorization. Shell syntax checks passed.

SuperMUC compilation/linking and the full job are still required.
Submission checks the compiled manifest, preserves it in results_JOBID,
and runs the comparator before printing END. The test partition limit is
30 minutes. The prior five-run one-link job took 15m40s; this extension adds
clover work, so that is not a runtime guarantee.
