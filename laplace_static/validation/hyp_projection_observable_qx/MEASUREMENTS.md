# Projection and scaling measurements, 6 October 2026

Projection: qX_p2n100, job 5914979, one temporal HYP step,
alpha=(1,1,0.5), 50069504 physical temporal links.

| Diagnostic | piter4 | piter40 |
|---|---:|---:|
| Maximum Frobenius unitarity defect | 1.3828922260540608 | 7.4985014987252858e-16 |
| RMS Frobenius unitarity defect | 4.4105126671388859e-3 | 2.6576002778170504e-16 |
| Maximum determinant defect | 1.7457549354048014 | 7.7718709382519141e-16 |
| RMS determinant defect | 3.6829946759775185e-3 | 1.7697174534518169e-16 |
| Links with unitarity defect >0.01 | 15678 | 0 |
| Links with determinant defect >0.01 | 12682 | 0 |

Round-trip, zero alpha, spatial/absent-boundary changes and independent repeat4
differences were zero. Maximum real/imaginary temporal component difference
between piter4 and 40 was 1.5568237997733603. The structural PASS does not establish
near-SU(3) projection at piter4. Review observable effects and the historical
finite-map convention before committing production settings.

Eigenvector scaling: same n110 input and binary, full 192-timeslice output.
Both completion checks passed at their documented file-inventory/solver-return
scope; independent residual/orthogonality checks are still outstanding.

| Nodes | Ranks | Grid | Observed seconds | Estimated physical core-hours |
|---|---:|---|---:|---:|
| 16 | 768 | 12x4x4x4 | 877 | 187.0933 |
| 32 | 1536 | 24x4x4x4 | 559 | 238.5067 |

The core-hour estimate is nodes*48*seconds/3600 for the Skylake nodes with
48 physical cores per node (LRZ hardware fact sheet). It also equals rank-hours
at the requested 48 single-thread ranks per node. It is an allocated-capacity
estimate, not measured busy CPU time or a verified project-budget debit.
The 16-node 877s includes its completion checker, while 32-node 559s is MPI
application time including output. These ratios are therefore approximate:
speedup 1.56887, parallel efficiency 78.4436%, walltime reduction 36.25998%,
core-hour increase 27.48005%.

For comparable same-configuration/full-write runs, minimizing core-hours per
successful configuration also maximizes throughput at a fixed core budget.
Two 16-node workers are estimated to outperform one 32-node worker in throughput,
but concurrent output effects have not been measured. Four 8-node workers should
be benchmarked before choosing production packing. An isolated 8-node run needs
to finish in roughly <1754s (29m14s) to beat the observed 16-node core-hour cost.

Collect job-level accounting separately; do not sum a job row with its step rows:

    sacct -X -j 5914674,5914980 --format=JobID,State,ElapsedRaw,AllocNodes,AllocCPUS,CPUTimeRAW,TotalCPU,AllocTRES%100

CPUTimeRAW/3600 is allocated Slurm CPU-hours; Slurm CPU units may be hardware
threads, so verify them before equating them with physical core-hours. TotalCPU
reports measured CPU consumption and may depend on how MPI tasks are accounted.
Exact billing follows LRZ's accounting policy. Queue waiting belongs in turnaround
comparisons, and is not execution core-hours.

Hardware reference:
https://jards.gauss-centre.eu/PublicFiles/FactSheets/FactSheet-GCS-SUPERMUC-NG.pdf
Projection reference: arXiv:1012.3037 section 2.1, and hep-lat/0506008 Eq.(2.24-25).
