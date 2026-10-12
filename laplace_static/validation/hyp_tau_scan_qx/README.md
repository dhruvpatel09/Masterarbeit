# Fixed-midpoint temporal HYP scan on qX_p2n100

Reuses the existing validated hyp_profile_scan_qx executable without rebuilding or modifying C. Its solver_mpiInit initialization is retained. Runs HYP_MODE=2, one HYP step, alpha=(1,1,0.5), piter=40, Nv=100 and R=1,2,4.

Endpoint pairs (97,99), (96,100), (95,101) give tau=2,4,6 and fixed midpoint 98. Uses 64 ranks on eight nodes, the same 1x4x4x4 process grid. Three serial invocations avoid increasing memory demand. Approximately 17 minutes is a planning estimate from the previous HYP40 run, not a guarantee; allocation limit is 30 minutes.

Checks source hashes and existing build manifest. Checks complete point records, path endpoints, factorization, halo counts, rank reductions, probe translations and all 332 profile records per tau. Tau4 is compared with job5903817. The original probe at the fixed midpoint is compared across tau. Produces 996 raw CSV records, retaining complex components. No baseline/HYP80 reruns or new eigenvectors required.

This is a single-configuration temporal-dependence measurement, not a statistical plateau determination. Original clover conventions remain unchanged. The wrapper does not independently validate the original observable implementation.

Run from Masterarbeit: sha256sum -c laplace_static/validation/hyp_tau_scan_qx/inputs.sha256, then sbatch laplace_static/validation/hyp_tau_scan_qx/submit.sh. Existing build and reference must be present. If a checksum fails, inspect the changed files rather than overriding the guard.

Review of supplied mental.c: all local timeslices are processed; no selected-timeslice input was found. A nonzero qcd_lanczosCheby return aborts MPI. Achieved residuals are not printed here. mental uses MPI_Init; this package does not modify or submit mental. Source inspection does not establish that the existing mental executable was built from exactly this source. Old n100 job timestamps span 11m39s.
