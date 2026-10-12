# qX temporal-HYP one-link observable integration

Additive package based on the uploaded test_tau_qX_one_point.c.
Run build.sh then sbatch submit.sh from the Masterarbeit root.
The baseline source SHA256 must match. No existing source is overwritten.
A new solver/build_supermuc_qx_p1x4x4x4 is built; build.sh refuses to overwrite it.

## Geometry and resources

Global 192x64^3, process grid 1x4x4x4, local 192x16^3, one thread.
64 ranks across 8 nodes (8 ranks/node) leave memory headroom for HYP work fields.
This is a diagnostic allocation, not a throughput-optimized production layout.
The 384-rank solver library cannot be used with this driver's geometry.
The test partition request is 30 minutes; the runtime is not yet measured.

The existing script's eigenvector path is retained. Preflight checks require
t96 and t97 HDF5 files. Nv=100, rho=1, configuration qX_p2n100.
The driver reads the two existing eigenvector slices without modifying them.

## Coverage

Original baseline plus disabled, zero-alpha, one-HYP-step piter40, piter80.
The new executable uses solver_mpiInit. The original executable keeps MPI_Init
as the control. The original gauge field is retained separately. Only the new
transport field is passed into the tau and factorized pair calculations, and
its positive spatial halo is refreshed before measurements.

Compare all 10000 tau entries at (0,0,0), selected entries, R=0 factorization,
and R=1 local (14->15) and cross-rank (15->16) pairs. Factorization PASS statuses
are enforced by the comparator. Numeric tolerance: 1e-18 + 1e-11*max(|a|,|b|).
Link-control/isolation tolerance: 1e-13. Nonfinite outputs are rejected.
Results are written to results_JOBID. The parent log reports progress and PASS.

This is a qX transporter and pair-correlator smoke test. It does NOT yet compute
qX clover insertion, arbitrary temporal extents, spatial averages, or production
flux-tube profiles. The tiny clover integration is already separately validated.
No independently known HYP observable value is claimed. A pass checks baseline
recovery, finite results, algebraic factorization, halo use, and 40/80 agreement
for these samples and this decomposition.

## Preparation checks

Source diff, shell syntax, Python syntax, and comparator tolerance rejection
were checked locally. MPI/HDF5/MKL compilation and numerical checks must run
on SuperMUC. No production parameter approval is inferred from this test.
