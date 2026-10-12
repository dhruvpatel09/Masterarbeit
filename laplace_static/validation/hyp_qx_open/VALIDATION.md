# Temporal HYP validation on qX_p2n100

Configuration: /hppfs/scratch/00/di98qes3/qX/qX_p2n100
Global lattice: 192x64x64x64
MPI grid: 6x4x4x4; 384 ranks; one OpenMP thread per rank
Initialization: solver_mpiInit (MPI_THREAD_FUNNELED)
Smearing: one temporal HYP step; alpha=(1,1,0.5)

## Results

- Job 5902645: piter=10 failed final SU(3) checks.
  Maximum temporal component difference, 10 vs 20: 0.72123965145790625.
- Job 5903127: piter=20 and 40 passed final SU(3) diagnostics,
  but differed by 0.026137617585367933.
- Job 5903146: piter=40 vs 80 passed.
  Independent repeat at piter=40: zero component difference.
  Maximum temporal component difference: 5.8841820305133297e-15.
  Unitarity and determinant defects: below 8e-16.
  Round-trip and zero-alpha differences: zero.
  Spatial and absent upper-boundary temporal link changes: zero.
  PASS and END printed; script duration 214 seconds.

## Scope

Use piter=40 as the working choice for integration on this configuration.
This establishes structural checks, repeatability, and 40-vs-80 stability
for this configuration and MPI decomposition. It is not an independent
verification of HYP or an ensemble-wide convergence study.

One HYP smearing step is distinct from the projection iteration count.
The advisor's intended smearing-step count remains to be clarified.
