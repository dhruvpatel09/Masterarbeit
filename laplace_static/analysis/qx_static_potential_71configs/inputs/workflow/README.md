# Fixed-source qX static-correlator scan

The production target is the complete supplied ensemble: all 71 configurations, n100–800 in steps of 10, including recovered n640. Existing saved eigenpairs are reused. No new eigenpair generation is requested.

The measured object is the existing Laplacian static correlator, with Nv=100 and rho=1, all spatial source sites and the three positive spatial axes. It is not a rectangular Wilson loop with thin spatial strings. The contraction and HDF5 readers are extracted unchanged from the validated qX profile source. This program computes no clover insertion.

## Measurement geometry

| Parameter | New scan | Preserved regression |
|---|---|---|
| Lattice | 192 x 64^3; open time, periodic space | Same |
| Spatial separation R/a | Every integer 1–32 | 1, 2, 4 |
| Temporal separation tau/a | Every integer 0–16 | 2, 4, 6 |
| Source time t0/a | Fixed at 90 | 98 − tau/2 |
| Eigenpairs | Existing 100 spatial APE20, alpha=0.5 modes | Same |
| Temporal HYP | One step, alpha=(1,1,0.5) | Same |
| Projection | p10; matched p40 in the pilot | p10/p40 archived correlators |

R/a=32 is half the periodic spatial box. The earlier values 2,4,6 were temporal separations. Half the spatial size does not prescribe a temporal limit. The source interval 90–106 stays near the temporal centre; extending tau beyond 16 requires a separate geometry and signal review.

After collection, effective potentials can be formed from `log(mean(L_tau)/mean(L_tau+1))` in lattice units, with covariance propagated across tau. A tau=0 point is a normalization diagnostic, not a preferred plateau-fit point. Negative or statistically unresolved means must be reported; do not replace them with absolute values. Fit ranges require signal and plateau checks. Periodic images and finite-volume effects remain relevant near R/a=32. Physical comparison with Em1 also needs a scale and an additive-potential convention.

## Pilot and complete-ensemble run

The pilot uses n170 (previous central-link outlier), recovered n640, and n650 (production control). Each of three array tasks runs a selected-time saved-mode audit followed by p10 and p40 scans. The audit covers every time slice from 90 through 106 and retains the relative residual limit 1e-9. It does not regenerate or modify eigenpairs.

Each scan exports 2,212 correlators (2,176 fixed-source and 36 original-geometry controls). All three axes and their pooled average are retained. It checks every shifted endpoint coordinate, 120 direct-versus-factorized contractions, and three full periodic payload roundtrips. The p10 and p40 final central-link defects must each be at most 1e-12. Their matched correlators must agree within `1e-20 + 1e-10*max(abs(x),abs(y))` in each real/imaginary component. These are explicit numerical validation tolerances for this workflow.

At selected tau and R, 64 regional means from 16^3 source cubes are also exported. These are correlated spatial-block diagnostics for investigating source averaging; they are not 64 extra independent configurations. Their pooled mean must reproduce the full spatial mean.

On SuperMUC, from the repository root:

```bash
python3 laplace_static/validation/qx_static_potential/test_workflow.py
bash laplace_static/validation/qx_static_potential/build.sh
python3 laplace_static/validation/qx_static_potential/campaign.py prepare --name pilot
python3 laplace_static/validation/qx_static_potential/campaign.py submit --name pilot --concurrent 3
```

When its array finishes:

```bash
python3 laplace_static/validation/qx_static_potential/campaign.py collect --name pilot
```

The review ZIP is written in the repository root as `qx_static_potential_pilot_review.zip`. Collection reparses the logs and checks their hashes. An incomplete run exports an explicitly incomplete review and exits with code 2.

After the pilot has passed and its results and timings have been inspected, prepare the full stage with `prepare --name full71 --stage full --pilot pilot`. The three validated p10 pilot measurements are reused. The remaining 68 configurations are packed two per array task only if the measured audit-plus-p10 timings leave at least 50% margin within 50 minutes; otherwise the controller uses one per task, or stops if even one would not fit. At most 39 array tasks may be submitted per batch. A previous submission must have left the queue before another batch or retry is submitted. Use `--start` and `--count` for a second batch if the prepared plan has more than 39 tasks. Full collection enforces complete coverage of all 71 configurations. There is no automatic acceptance of a partial production cohort.

## Reproducibility and local checks

`INPUTS.json` pins the parent kernels, saved-mode audit source, input inventory, original-geometry reference table and this workflow. The build is isolated and verifies the prior profile build's library/environment manifest. Its own manifest records the generated code, executables, libraries and headers. Preparation records physical-input paths, symlink targets, sizes and timestamps, plus hashes for the small eigenvalue table and plan. Large gauge/eigenvector files are not copied into Git.

The standard-library test suite checks complete geometry and data coverage, original-correlator regressions, strict residual rejection, projection checks, and immutable case provenance. Native C tests execute the extracted complex contractions against an independent Nv^2 sum and execute the exact face-shift implementation on a simulated 4^3 process grid through 64 periodic shifts, checking complete payloads at both process-grid corners. A syntax check uses the preserved qcd interfaces and test-only MPI/HDF5-reader declarations. The production build uses the actual MPI and HDF5 libraries; execution on SuperMUC is the final distributed validation.

GitHub synchronization continues through Stromboli. Do not fetch or push GitHub directly from SuperMUC. Keep the original recovered-mode directories: assembled n640 files contain symlinks to those data.
