# Temporal-HYP observable integration: tiny open lattice

This is an additive test package for the existing Masterarbeit repository.
It does not modify the existing drivers, solver sources, or builds.

## Run on SuperMUC

From the Masterarbeit root:

    bash laplace_static/validation/hyp_observable_tiny/build.sh
    sbatch laplace_static/validation/hyp_observable_tiny/submit.sh

The build checks the SHA256 of the existing baseline source, then creates
solver/build_supermuc_tiny_serial (one rank, local 8x4x4x4, one thread).
It stops if that build directory already exists. A failed build can be resumed
manually after diagnosing its error; do not blindly delete existing builds.

The submission uses partition=test, qos=test, one node, one rank, 10 minutes.
It reads the existing tiny gauge configuration and t=2,6 HDF5 eigenvectors.
Five executions: original baseline, solver roundtrip with HYP disabled,
zero-alpha HYP, one temporal HYP step with piter=40, and the same with piter=80.
Outputs are in results_JOBID and hyp_obs_JOBID.out/.err alongside these files.

## Conventions

Nv=16, t_src=2, t_sink=6, R=2, probe=(0,1,2), all 192 axis/source samples.
Only the temporal transporter receives the separate transport gauge field.
Clover and plaquette use the original gauge field. Existing eigenvectors are
read without modification. The current baseline uses unprojected colour
clover F=(Q-Qdagger)/8 and S_raw=sum_{mu<nu} Tr(Fdagger F), without a factor 1/2.
The unsmeared L reference check stays active in the disabled/zero-alpha runs.
It is inapplicable in the nonzero-HYP runs. The clover reference stays active.

The comparator rejects nonfinite values, compares all printed baseline/control
records and 40/80 records, and checks original-probe values across all modes.
Observable comparison uses atol=1e-13 plus rtol=1e-12; probe isolation uses
atol=1e-14. Link isolation/control checks use 1e-13. Existing baseline regression
tolerances are retained. No known HYP observable reference is asserted.

This is an integration check, not an independent HYP implementation validation,
a production qX driver, or an ensemble-wide convergence determination.
The baseline executable retains its original MPI_Init intentionally as the
control; the new integration executable uses solver_mpiInit.

## Local preparation checks

Source replacements were checked against the uploaded baseline. Shell syntax,
Python syntax and comparator rejection behaviour were checked locally.
The MPI/HDF5/MKL compiler environment is not available locally; compilation
and numerical execution must be verified on SuperMUC.
