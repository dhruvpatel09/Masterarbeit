# qX separation and probe-plane scan

Additive extension of spatial averaging validated in job 5903813.

Build and submit from Masterarbeit:

    bash laplace_static/validation/hyp_profile_scan_qx/build.sh
    sbatch laplace_static/validation/hyp_profile_scan_qx/submit.sh

Requires the existing solver library for grid 1x4x4x4 and the baseline.out /
hyp40.out files in hyp_spatial_average_qx/results_5903813. The script performs
only two runs: original-field baseline (mode 3) and one temporal HYP step
with alpha=(1,1,0.5), piter=40 (mode 2). HYP and gauge/eigenvector reads occur
once per run; transport is reused for all separations and probes. Prior
zero-alpha/off/piter80 validation remains recorded in job 5903813.

## Fixed case and output

qX_p2n100, lattice 192x64^3, Nv=100, uniform weights, t0=96, t1=100,
tmid=98. 8 nodes, 64 MPI ranks, one OpenMP thread each, solver_mpiInit,
test partition/QoS, 30-minute limit. Existing eigenvectors are reused.

R=1,2,4, with all 64^3 sources and three positive axes, 786432 combined pairs
for each R/probe. The original raw clover convention is retained:
F=(Q-Qdag)/8, no colour-trace subtraction; S_RAW=sum six Tr(Fdag F), no half.

For each R the probe set contains:
- Regression offset (long,tr1,tr2)=(0,1,2).
- Plane: long=-1..R+1, tr1=-2..2, tr2=0.

Coordinates are in lattice units, relative to the adjointed source endpoint.
Axis a uses spatial directions a,(a+1)%3,(a+2)%3 for long,tr1,tr2, with
periodic spatial wrapping. long_from_midpoint=long-R/2 is exported explicitly.
For R=1 the spatial midpoint is half-integer; the code never rounds it.
The plane is a deliberately small implementation check, not a finalized
physical fit/profile range.

Outputs are raw complex mean L, real mean S_RAW, and complex mean LS,
where LS is multiplied at each source before averaging. No connected
subtraction, ensemble errors, or ensemble averaging is performed.
After comparisons pass, results_JOBID/profile_raw_qX_p2n100_t96_100.csv
contains 664 rows: 332 each for baseline and hyp40. axis=-1 means combined.
Keep separate per-configuration data for the later ratio of ensemble means.

## Endpoint shifts

Each local site is packed into an endpoint payload containing both eigenvector
slices, the temporal product and global coordinates. Rank neighbours are
resolved from MPI_Allgather of actual geometry positions, not rank arithmetic.
A one-cell shift exchanges the negative face with the negative neighbour and
receives the positive neighbour's negative face using MPI_Sendrecv/MPI_BYTE.
Interior cells are shifted locally. Repeating up to four times supplies the
R=1,2,4 endpoints without indexing beyond one-layer qcd halos. Global
coordinates are checked at every site after every shift, including step 3.
Additional payload and face buffers use about 125 MiB per rank at Nv=100.

## Validation gates

- Existing point and R=1 source-average code is retained. Its outputs are
  compared against job 5903813 for baseline and HYP40.
- New packet-shift R=1 averages at offset (0,1,2) must agree with the retained
  one-layer-halo implementation in the same run, for all axes and combined.
- Each R/axis has 262144 sources, R*16384 cross-rank endpoints and R*4096
  periodic wraps; wraps are a subset of crossings for the selected R values.
- Three direct Nv^2 versus colour-factorized comparisons per R/axis exercise
  an interior pair, the first cross-rank pair, and the first global wrap.
- All expected offsets must occur exactly once; values must be finite.
- Mean L must be independent of probe position.
- Every mean S must equal the original spatial clover mean (translation).
- Combined results must equal the mean over the three axis results.
- Baseline and HYP use identical original-field probe densities.

The comparison imports point_compare.py and average_compare.py, which preserve
the previously used checks. Build and submission manifests record sources,
comparators, scripts, libraries and reference outputs. Files from earlier
validation packages are not modified.

## Local verification

GCC C11 -Wall -Wextra -Werror syntax check passed with uploaded project headers
and minimal MPI/HDF5 declaration stubs; this is not production linking.
The actual profile_scan and shift_endpoint functions passed a 64-block
simulated-exchange test with Nv=3 synthetic data. Every outgoing face was
checked against expected coordinate-dependent payload values; neighbours,
tags, received coordinates, and the four sequential shifts were checked.
All 249 per-axis profile rows agreed with an independent NumPy full periodic
lattice calculation, and the comparator checked 83 combined rows. It rejected
wrong crossing counts, a missing PASS and an incorrect offset set.
A separate complex-valued test compared packet_pair and the direct contraction.
Shell syntax and Python compilation checks passed.

These are local mathematical/control tests. SuperMUC must still establish
correct execution with its actual MPI and qcd libraries. Require final PASS
and END in qx_hypscan_JOBID.out and inspect stderr before using the CSV.
