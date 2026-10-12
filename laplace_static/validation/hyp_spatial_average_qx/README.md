# qX spatial-source and three-axis averaging

This additive validation extends the four-link case passed in job 5903809.
The earlier sources/results are retained and used for pointwise regressions.

## Run

From Masterarbeit:

    bash laplace_static/validation/hyp_spatial_average_qx/build.sh
    sbatch laplace_static/validation/hyp_spatial_average_qx/submit.sh

Requires the existing solver/build_supermuc_qx_p1x4x4x4/libsolver.a,
qcd/lib/libqcd.a, and all five .out files in
laplace_static/validation/hyp_multilink_qx/results_5903809/.
The source/header checksums must match the files used to prepare this package.

The job uses 8 nodes, 64 MPI ranks, one thread/rank, test partition/QoS,
30-minute limit, solver_mpiInit (MPI_THREAD_FUNNELED), and the existing
qX_p2n100 eigenvectors at t=96 and 100. It runs baseline, off, zero, hyp40,
and hyp80 sequentially. HYP is one step, alpha=(1,1,0.5).

## Measurement

- Global lattice 192x64^3; process grid 1x4x4x4.
- R=1 and tau=4; ordered temporal product t=96..99, midpoint t=98.
- Nv=100, constant mode weights; existing spatial eigenvectors are reused.
- Every spatial source contributes once per axis: 262144 pairs per axis,
  786432 combined. Sources are correlated samples, not independent configs.
- Three positive spatial directions; no separate negative-direction average.
- Axis a probe offset: 0 along a, +1 along (a+1)%3, +2 along (a+2)%3,
  modulo 64. This rotates the established offset (0,1,2).
- Original gauge field supplies F=(Q-Qdag)/8 and S_RAW=sum of six
  Tr(Fdag F), no colour-trace subtraction and no factor 1/2.
- HYP field supplies transporters only.
- The 2 MiB midpoint scalar-density slice is replicated by MPI_Allreduce;
  arbitrary rotated probe positions need no scalar corner halos.
- Pair L is computed by the tested colour-kernel factorization, O(Nv).
  Complex LS is formed at each source BEFORE summation.
- Local sums use compensated summation; global sums use MPI_Reduce.
- Per-axis and combined outputs contain mean L, mean S_RAW, and mean LS.
  No ratio, connected subtraction, ensemble error, or physical flux claim.

## Gates

1. Existing baseline/off/zero and piter40/80 comparisons remain.
2. All existing point outputs match job 5903809 within the comparator tolerance.
3. Coordinate sentinel fields independently check qcd plus3d halo indexing
   along x, y, z, including periodic global wraps.
4. Direct Nv^2 and factorized results are checked at three sources per axis:
   interior coordinate 14, MPI-face coordinate 15, global-wrap coordinate 63,
   with the two other coordinates zero.
5. Per-axis counts: interior 245760, halo 16384, periodic wraps 4096
   (wraps are a subset of halo pairs). Combined count is 786432.
6. Root prints all 192 rank-axis partials. Python independently re-sums them
   with math.fsum and checks MPI reduction and normalization.
7. Probe means for each axis equal the spatial clover mean, since each
   fixed-offset translation is a bijection of the periodic spatial lattice.
8. Original probes remain unchanged for all five modes.

The coordinate test validates the spatial spinor halo maps. The real MPI
execution is still required; this is not an independent HYP implementation
or a decomposition-independence study.

## Local verification performed

- GCC C11 -Wall -Wextra -Werror syntax check using uploaded project headers
  and minimal MPI/HDF5 declaration stubs (not a production link test).
- Actual averaging helper tested on 64 simulated spatial blocks with Nv=3
  synthetic eigenvectors and periodic face buffers; its axis averages agree
  with an independent NumPy full-lattice calculation of L, S and LS.
  This exercises rotated offsets, normalization, and periodic wrap handling,
  but does not execute the MPI library or qcd communication implementation.
- Output parser independently re-summed synthetic rank partials and rejected
  corrupt source, halo and combined counts.
- Bash syntax and Python compilation checks passed.

The remote build and Slurm job are the remaining validation gate. Inspect
qx_hypavg_JOBID.out and .err, requiring the final PASS and END. Retain
results_JOBID and its manifests. The previous multi-link five-run job took
16m44s; new averaging costs additional time, with no guaranteed runtime.
