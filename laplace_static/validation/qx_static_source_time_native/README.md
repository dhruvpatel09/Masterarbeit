# qX source-time native pilot

This is the native implementation following the design and input preflight at
commit `2a84fa84abe68fa7768a54e91d769edd69fa1e50`. It leaves the published
`qx_static_source_time` design, the full71 measurements and their analysis
unchanged. The design's `native_extension_implemented: false` records its
historical stage; this separate package supplies that extension.

The immediate question is whether the measured matrix and the same frozen
source projections change when their physical temporal interval is moved.
The main radii are R/a = 24 and 32; R/a = 4 and 12 are controls. This is a
three-configuration technical pilot, not an ensemble potential determination.

## Fixed scope

The lattice is 192 x 64^3, with open time and periodic space. Source coordinates
are 58, 90 and 122, and every separation tau/a = 0,...,32 is retained. Thus the
largest end coordinate is 154; temporal propagation never wraps around an open
boundary. Spatial half-box distance is 32; temporal half-extent is 96.

| Configuration | Audited slices | New scans, in execution order |
| --- | ---: | --- |
| n170 | 66 | t0=58 p10; t0=122 p10 |
| n640 | 66 | t0=58 p10; t0=122 p10 |
| n650 | 97 | **t0=90 p10 regression first**; t0=58 p10; t0=122 p10; t0=122 p40 control |

There are three array tasks, eight scans and three audits: **eleven MPI
launches**. Each task requests eight nodes, 64 ranks and one thread per rank,
with a one-hour wall limit. Maximum concurrency is three tasks.

Each scan retains all 16 ordered source/sink cutoff pairs (40,60,80,100), four
radii (4,12,24,32), all 33 separations, three spatial axes and their average:
8,448 matrix rows per scan; **67,584 new matrix rows** overall. All 100 saved
modes and the original MPI endpoint payload are retained. Rank blocks now
include R/a=24 alongside 4 and 32, at tau/a=0 and 32: 18,432 rows per scan,
**147,456 rows** overall. Rank blocks are diagnostics, not extra independent
configurations.

The existing t0=90 matrices on all three configurations are reused. The
additional fresh n650 t0=90 scan must reproduce all **2,112 axis-averaged
ordered matrix entries**, including off-diagonal complex values, before its
shifted scans start. The shifted p10/p40 control uses identical t0=122
endpoints and compares all 8,448 entries.

The 36 frozen sources (four radii, three bases and three training times) come
from the published 68-configuration training sample, which excludes n170,
n640 and n650. Exactly the same coefficients are used at every source
position, with sink-row/source-column convention `v^dagger M v`. There is no
pilot retuning. Full collection exports 13,068 fixed-source rows and 9,504
paired differences, distinguishing reused baseline, fresh regression,
source-time shift and projection control.

## Native checks

`scan.c` changes the validated matrix driver only where physical source time
or the additional R=24 rank-block inventory is required. Both endpoints are
written in every matrix, direct-check and block record. `MODE_LOAD` records
the actual HDF5 source and sink coordinates: exactly 33 loads per scan.
The transporter product and SU(3) diagnostics use that same physical window.

The unchanged contraction, I/O, temporal-product and MPI helper functions are
extracted from pinned published sources. The saved-mode audit changes only
the selected-time list relative to the previously validated strict audit.
Each of the 100 modes on every required slice must pass residual and norm
criteria at 1e-9, Gram defects at 1e-9 and imaginary Rayleigh bound at 1e-10.
The checker also recomputes reported residual relations and slice maxima.

Other checks retain the existing direct-sum tolerance
`1e-20 + 5e-12 * scale`, legacy diagonal checks, axis/rank reductions,
equal-time Gram tests and complete periodic payload roundtrips. Temporal HYP
must preserve the other link directions; shifted-window SU(3) maximum
unitarity and determinant defects must be at most 1e-12. The shifted p10/p40
matrix comparison uses `1e-20 + 1e-10 * scale`.

Negative signed correlators and finite-time complex/non-Hermitian matrices
are retained. There is no absolute-value replacement, finite-time
symmetrization, logarithm, seven-point covariance fit, GEVP or final potential
fit in this pilot. Three technical configurations do not establish the
statistical significance of a source-position or boundary effect.

## Local validation and provenance

Run:

```bash
qx_native=laplace_static/validation/qx_static_source_time_native
python3 "$qx_native/verify_files.py"
python3 "$qx_native/test_workflow.py"
python3 "$qx_native/campaign.py" check-inputs
```

The suite covers analytic source-position-dependent output replay, wrong
coordinates and missing slices, all ordered pairs, rank reductions and strict
failure paths. It exercises an entire synthetic prepare/run/retry/reuse/
collect cycle, submission receipts and stopped-job lock recovery. It compiles
the mixed contractions at 4 and 100 modes, periodic payload tests and an
independent noncommuting temporal-path test at every planned separation. The
driver is syntax-checked against the captured interfaces. Python uses 3.6
syntax and the standard library; local C tests require GCC.

These checks do not replace the real SuperMUC MPI/HDF5 build or the native
saved-mode and scan runs. Those have not been performed by the package's
local synthetic validation.

`SOURCE_PINS.json` pins 74 published source/reference files at the design
commit. `FILE_MANIFEST.json` covers this package and those dependencies.
`build.sh` additionally checks the archived clover build pins, records the
actual compiler/modules/linker output and hashes native binaries, generated
sources, headers and libraries. Prepare binds the build and input inventory;
run rechecks them. Source files and binaries must remain unchanged throughout
a prepared campaign.

## Publish from Stromboli

Place `qx_static_source_time_native_workflow.zip` in the repository root.
These commands add only the new implementation directory:

```bash
(
    set -e
    cd /home/m2130292/Masterarbeit
    git pull --ff-only origin main
    unzip -n qx_static_source_time_native_workflow.zip
    qx_native=laplace_static/validation/qx_static_source_time_native
    python3 "$qx_native/verify_files.py"
    python3 "$qx_native/test_workflow.py"
    python3 "$qx_native/campaign.py" check-inputs
    git add -- "$qx_native"
    git commit -m "Add native qX source-time validation pilot"
    git push origin main
    git bundle create qx_source_time_native_sync.bundle main \
        ^2a84fa84abe68fa7768a54e91d769edd69fa1e50
    git bundle verify qx_source_time_native_sync.bundle
    scp qx_source_time_native_sync.bundle \
        di24gus@skx.supermuc.lrz.de:/hppfs/work/pn29se/di24gus/Masterarbeit/
)
```

## Build and prepare on SuperMUC

```bash
(
    set -e
    cd /hppfs/work/pn29se/di24gus/Masterarbeit
    git bundle verify qx_source_time_native_sync.bundle
    git fetch ./qx_source_time_native_sync.bundle main
    git merge --ff-only FETCH_HEAD
    qx_native=laplace_static/validation/qx_static_source_time_native
    python3 "$qx_native/verify_files.py"
    python3 "$qx_native/test_workflow.py"
    git update-ref refs/remotes/origin/main "$(git rev-parse HEAD)" \
        2a84fa84abe68fa7768a54e91d769edd69fa1e50
    rm -- qx_source_time_native_sync.bundle
    bash "$qx_native/build.sh"
    python3 "$qx_native/campaign.py" prepare --name pilot
    git status
)
```

Building and preparing do not submit a job. When ready to allocate the pilot:

```bash
python3 laplace_static/validation/qx_static_source_time_native/campaign.py submit \
    --name pilot --concurrent 3
squeue -u "$USER" -o '%.18i %.20j %.2t %.10M %.6D %R'
```

A planning estimate is roughly **20--45 minutes after allocation** for the
pilot's slowest task, with queue delay additional. Previous p10 scans took
about 3.6--3.9 minutes each; the 66/97-slice audits, new output volume and
filesystem load are not yet timed here. This estimate is provisional; the
native launch receipts will provide measured timings.

## Collect and inspect

After all array tasks have left the queue:

```bash
qx_native=laplace_static/validation/qx_static_source_time_native
python3 "$qx_native/campaign.py" status --name pilot
python3 "$qx_native/campaign.py" collect --name pilot
ls -lh qx_static_source_time_pilot_review.zip
git status
```

Collection replays native stdout, stderr hashes, launch order, physical
endpoint inventories and validated CSV exports. The archive is created at
the repository root and includes source/design/build receipts, case outputs
and failed attempts. Empty stderr files are retained. Binaries are omitted.
Missing cases produce an explicitly incomplete archive and a nonzero exit;
only 3/3 validated cases produce a successful pilot collection.

Inspect the fresh regression, saved-mode checks and projection control first,
then compare signed/complex matrices and frozen-source changes between
positions at each radius. Decide on an ensemble extension only after these
technical checks and timings. No production campaign is submitted
automatically.

Retries use new attempt directories and preserve the old evidence. Duplicate
or ambiguous submissions are blocked. If a wall-time termination leaves a
`RUNNING.json`, do not remove it manually: `campaign.py recover-lock --name
pilot --index INDEX` checks the recorded job with `squeue`, refuses recovery
if it may still be active and archives the stopped-job receipt. Collection
refuses running locks. A failed build remains available for diagnosis; do
not alter a build attached to a prepared campaign.

After the SuperMUC import is verified, remove the staging ZIP and sync bundle
on Stromboli. Retain the review archive until it has been transferred and its
analysis inputs have been verified.
