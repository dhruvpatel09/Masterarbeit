# Continue the qX projection comparison across separate job IDs

Job 5918419_0 has passed: n110 at tau2, with the fresh n100 tau4
regression matching all 664 reference records exactly. This helper continues
indices 1..14 of the already prepared comparison, without repeating index0.
It changes no measurement source, executable, library, prepared input or
scientific acceptance rule and creates no new eigenvectors.

## Submission and collection

The controller submits **one single-element array at a time**, using the original
8-node, 64-rank, test-partition script. It submits the next case only after the
previous job has left the queue and its actual output checks have passed.
Each accepted parent ID is recorded in RUNS.json. The final collector handles
these distinct parents and requires all 15 cases plus the successful n100
regression. The old collector assumes one parent and should not be used here.

The original `--array=0-14%2` could hit the QoS submitted-job limit: `%2` limits
running tasks, not the number of submitted array elements. This helper leaves
at most one of its cases pending/running. It does not require the accounting
database or call sacct. It queries the user inventory with squeue, and treats
a query error as a stop rather than as proof of completion.

Indices 1 and 2 cover n110 at tau4/6. Indices 3..14 cover n120, n130, n140 and
n150, each at tau2/4/6. All use midpoint98, R1/2/4, Nv100, one temporal HYP step,
alpha=(1,1,0.5), independent piter4/40 measurements and original-field clover.
The original profile checks and piter40 SU3 control are reused unchanged.
No piter4 physics acceptance threshold is added.

## Start from Masterarbeit

Unpack this additive package once. Do not rebuild or prepare the measurement
again, and do not resubmit the full original array.

```bash
cd /hppfs/work/pn29se/di24gus/Masterarbeit
test ! -e laplace_static/validation/hyp_projection_configs_qx_queue
unzip hyp_projection_configs_qx_queue_integration.zip
sha256sum -c laplace_static/validation/hyp_projection_configs_qx_queue/package.sha256
sha256sum -c laplace_static/validation/hyp_projection_configs_qx_queue/inputs.sha256

nohup python3 -u \
    laplace_static/validation/hyp_projection_configs_qx_queue/run_remaining.py \
    --first-job 5918419 \
    > laplace_static/validation/hyp_projection_configs_qx_queue/controller.out \
    2>&1 < /dev/null &
```

The controller is lightweight submission/bookkeeping work; every MPI
measurement still runs as a Slurm compute job. It polls every 30 seconds.
nohup allows the controller to continue after a normal terminal disconnect;
local login-session policies can still terminate it, so the workflow is resumable.

Monitor with:

```bash
tail -n 25 laplace_static/validation/hyp_projection_configs_qx_queue/controller.out
squeue -u "$USER" -o '%.18i %.14j %.2t %.10M %.6D %R'
```

Success ends with COMPLETE and creates this review archive **inside Masterarbeit**:

`hyp_projection_configs_review_mixed_5918419.zip`

It contains the aggregate raw/paired/SU3 CSVs, report, case-to-job mapping,
individual logs, measurement sources and provenance. It excludes eigenvectors,
HDF5 data, executables and libraries. Download and share that exact ZIP before
deleting it. Its final SUMMARY requires 9960 raw rows, 4980 paired rows,
120 SU3 time records and 30 SU3 window records.

## Resume and failure handling

The controller uses a lock to prevent two copies from submitting the same
campaign. It flushes the ledger before and after sbatch. After an interrupted
controller, or a temporary squeue failure, keep RUNS.json and resume:

```bash
nohup python3 -u \
    laplace_static/validation/hyp_projection_configs_qx_queue/run_remaining.py \
    --resume \
    >> laplace_static/validation/hyp_projection_configs_qx_queue/controller.out \
    2>&1 < /dev/null &
```

Known pending/running jobs are monitored rather than submitted again. If sbatch
explicitly rejects a case with QOSMaxSubmitJobPerUserLimit, that case remains
unsubmitted; resume after other submissions finish. Any MPI failure, failed
scientific/structural check or changed provenance stops further allocation.
Failed measurements are preserved and are not automatically retried.

If interruption occurred after Slurm may have accepted a submission but before
its reply was durably recorded, the controller stops on an ambiguous submission
instead of guessing. Keep the ledger and logs and resolve that case before
continuing. Stopping the controller does not cancel already accepted jobs.

After the controller stops, a partial review can be generated explicitly:

```bash
python3 laplace_static/validation/hyp_projection_configs_qx_queue/collect_mixed.py --partial
```

This is labelled partial and exits nonzero. Running collect_mixed.py without --partial
requires complete coverage and can safely reuse an already completed, verified
archive. Neither a queue disappearance nor a completion marker alone is taken
as proof of successful output: the collector rechecks the original logs.

## Timing and interpretation

The first n110 piter4/40 pair took 424.8 seconds, approximately 7m05.
Fourteen similar pairs would take about 99 minutes of execution, plus queue
waits, polling and checking; configuration-dependent work and shared I/O can
change that. The 30-minute limit remains per case. This estimate does not
predict queue wait time.

The n110 tau2 piter4 central-window maximum Frobenius unitarity defect was
1.4431382465295746e-6, determinant defect 7.2528849170972100e-7, with no links
above 1e-2 in that window. The largest relative combined-axis connected
diagnostic change was 7.898347177072421e-6 (0.0007898%). These are promising
single-configuration diagnostics, not an ensemble bias or uncertainty estimate.
The earlier large open-boundary tails outside these central windows are not
made acceptable by this result. The remaining cases extend the central-window
and observable-sensitivity checks before deciding the production prescription.

Connected values are per-configuration diagnostics. Ensemble analysis must
average the raw L, LS and S before forming the ratio/subtraction and must
account for correlations. Spatial sites, axes and tau values are not independent
ensemble configurations. The saved-mode audit still covers nine selected
slices per configuration; this helper adds no full eigenmode residual audit.

## Local verification scope

The helper is tested locally with simulated Slurm replies and fixture result
trees: sequential admission, resume without duplicate submissions, submit-limit
rejection, ambiguous submission handling, queue failure, mismatched provenance,
mixed-parent collection, missing coverage and archive integrity. It does not
launch a SuperMUC job locally. The unchanged measurement checks run again on
the actual outputs on SuperMUC.

Slurm references:
- https://slurm.schedmd.com/job_array.html
- https://slurm.schedmd.com/sbatch.html
