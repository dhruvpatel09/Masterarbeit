# qX static correlators: temporal extension and mode cutoffs

The published 71-configuration scan reaches temporal separation `tau/a=16`.
Several effective-energy fits change when the fit window moves to later times.
The weak signal near `R/a=24` also warrants an operator-overlap check before
interpreting a fitted value as a static potential. This workflow extends the
time range and compares projectors formed from subsets of the existing modes.
It does not generate eigenpairs or change the clover measurements.

## Measurement

- All 71 configurations, n100 through n800 in steps of 10, remain the production
  target. The technical pilot uses n170, recovered n640, and n650.
- Source time `t0/a=90`; `tau/a=0,...,32`; sink times `90,...,122`.
- Spatial distances `R/a=1,2,4,8,12,16,20,23,24,25,28,32`, along all three
  periodic axes. Every spatial source is included; the combined-axis result is
  their average.
- One temporal HYP step with alpha `(1,1,0.5)`, `piter=10`, `rho=1`.
- `Nv=40,60,80,100` are summed in the same p10 launch, after reading the same
  saved 100-mode endpoint buffers. `piter` counts projection iterations;
  `Nv` counts eigenvectors used in the projector.
- Each pilot configuration additionally has a `piter=40`, `Nv=100` launch
  over the same temporal and spatial geometry.
- Geometry `192 x 64^3`, rank grid `1 x 4 x 4 x 4`, 64 ranks on eight nodes,
  one thread per rank, MPI FUNNELED. The existing gauge and saved-mode paths,
  including the recovered n640 symlinks, are taken from the admitted inventory.

For a cutoff N, the endpoint projector is
`K_N(x,y;t) = sum_(n=0)^(N-1) v_n(x,t) v_n(y,t)^dagger`.
Changing N changes the source operator and its state overlaps. The amplitudes
can therefore differ substantially without indicating a different energy.
Compare normalized correlators and effective energies with paired errors before
comparing fitted plateaux. The sequence is not an extrapolation to an exact
answer at `Nv=100`, and smaller N is not guaranteed to improve every distance.
At half the spatial box, periodic-image effects also remain part of the physics.

## Controls

Before either scan, all 100 modes on each of the 33 selected time slices are
checked against their saved eigenvalues. Relative residual, norm and Gram
tolerances remain `1e-9`; eigenvalue ordering is checked before admitting a
cutoff interpretation. Temporal transporter SU(3) diagnostics cover links at
`t/a=90,...,121`, including every new time segment. Their numerical maximum
tolerance remains `1e-12`; spatial and top-boundary links must be preserved.

Every Nv=100 fixed-source record through tau=16 is compared with the published
matched configuration and piter. The original midpoint-98 tau=2,4,6 controls
are also retained at Nv=100. The reference CSV contains 63,048 archived records:
p10 for all 71 configurations and p40 for the three pilot controls. Source
paths, hashes and the publication commit are in `inputs/REFERENCE_SOURCES.json`.
Contraction comparisons use absolute tolerance `1e-20` plus relative `5e-12`;
the pilot p10/p40 comparison uses relative `1e-10`. These tolerances allow
roundoff from compilation and reductions without changing the measurement.

The direct O(Nv^2) contraction is checked independently of its factorized
implementation at four fixed-source times and four distances, for each cutoff.
The driver verifies every shifted coordinate and performs a complete 64-cell
payload roundtrip along each periodic axis. It reads and transmits the original
100-mode payload for every cutoff. Regional outputs are retained at Nv=100
for tau=0,8,16,24,32 and R=4,8,16,24,32.

Per configuration, p10 yields 6,372 correlator records, 396 direct checks and
4,800 regional records. The p40 control yields 1,620, 108 and 4,800 respectively.
The audit contains 3,300 individual mode records. Collection replays stdout,
checksums and coverage; a zero process exit code alone does not admit a case.

The statistical unit remains one configuration after source and axis averaging.
Use the advisor-confirmed independent-configuration baseline, with paired
covariances across times, cutoffs and controls. Regional outputs are correlated
diagnostics. There are at most 71 configuration samples, including at different
cutoffs and temporal separations. Negative individual correlators are retained;
only statistically supported positive ensemble means admit logarithms.

## Pilot on SuperMUC

After synchronizing the committed workflow:

```bash
cd /hppfs/work/pn29se/di24gus/Masterarbeit
python3 laplace_static/validation/qx_static_extension/test_workflow.py
bash laplace_static/validation/qx_static_extension/build.sh
python3 laplace_static/validation/qx_static_extension/campaign.py prepare \
  --name late_modes_pilot
python3 laplace_static/validation/qx_static_extension/campaign.py submit \
  --name late_modes_pilot --count 3 --concurrent 3
```

This is three array tasks, each with a strict audit and two scan launches.
The p10 launch evaluates all four cutoffs. The allocation has a one-hour limit;
actual audit and scan durations are recorded for production batching. Queue
waiting time cannot be estimated from this workflow.

After the array finishes:

```bash
python3 laplace_static/validation/qx_static_extension/campaign.py collect \
  --name late_modes_pilot
```

The resulting `qx_static_extension_late_modes_pilot_review.zip` is written
atomically in the repository root. Review its matched baseline, projection
control and timing results before preparing production. Three pilot
configurations establish technical validity and runtime; they do not replace
the full-ensemble physics analysis.

## Full-ensemble continuation

After pilot review, `prepare --name late_modes_full71 --stage full
--pilot late_modes_pilot` validates and reuses the three p10 pilot cases. It
plans only the remaining 68 configurations. A 50% timing margin determines
whether an array task contains one or two configurations. If even one will not
fit safely, preparation stops for a batching review.

Submission accepts at most 39 tasks per batch and concurrency at most four.
If the prepared plan has 68 tasks, submit tasks 0..38, wait for that batch to
finish, then submit tasks 39..67. If it has 34 tasks, one batch is sufficient.
Completed validated tasks are reused; active submissions block duplicate work.
Keep the case directories and saved modes until the results are verified and
published. Collection includes the pilot records and their admission markers.

Numerical analysis, plotting and Git publication can run on Stromboli after
collection; gauge-field MPI measurements run on SuperMUC. Transfer bundles
inside the repository root and remove them after the destination import and
verification succeed. The workflow supports Python 3.6 and requires no NumPy
for preparation, collection or its Python coverage tests. Local native tests
use simulated MPI face exchanges and preserved interfaces; the SuperMUC pilot
provides the actual MPI/HDF5 execution check.
