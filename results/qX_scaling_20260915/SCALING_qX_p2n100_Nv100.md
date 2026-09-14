# qX p2n100 Nv=100 process-grid scaling

Date: 2026-09-15

Configuration:
qX_p2n100

Global lattice:
192 x 64 x 64 x 64

Nv:
100

Lanczos parameters:
res              = 1e-9
redTol           = 1e-12
deg              = 10
upper            = 0.15
tolReortho       = 1e-10
maxRestarts      = 30
nLanczosFactor   = 7
redFreqFactor    = 10

All runs used:
MENTAL_EVAL_ONLY=1
write_ev=0

## Scaling results

| Nodes | MPI ranks | Process grid | Local lattice | Wall time | Speedup | Efficiency | Node-hours |
|------:|----------:|--------------|---------------|-----------|--------:|-----------:|-----------:|
| 4     | 192       | 3 x 4 x 4 x 4  | 64 x 16 x 16 x 16 | 27:38 | 1.000 | 100.0% | 1.842 |
| 8     | 384       | 6 x 4 x 4 x 4  | 32 x 16 x 16 x 16 | 15:14 | 1.814 | 90.7%  | 2.031 |
| 16    | 768       | 12 x 4 x 4 x 4 | 16 x 16 x 16 x 16 | 09:05 | 3.042 | 76.1%  | 2.422 |

All three runs converged on 192 / 192 temporal slices.

Mean per-timeslice eigensolver times:

4 nodes:
23.444 s

8 nodes:
23.029 s

16 nodes:
23.363 s

This confirms that the cost of one 3D Laplacian solve is essentially
unchanged. The wall-time scaling arises primarily from processing more
temporal slices concurrently.

## Numerical consistency against calibration reference

4 nodes:
compared      = 192
missing       = 0
mean abs diff = 1.83230167150050249e-17
max abs diff  = 1.11022302462515654e-16 at t=150

8 nodes:
compared      = 192
missing       = 0
mean abs diff = 2.17201835221262732e-17
max abs diff  = 1.24900090270330111e-16 at t=60

16 nodes:
compared      = 192
missing       = 0
mean abs diff = 0.00000000000000000e+00

The 4- and 8-node differences are at floating-point roundoff level and
are consistent with changed MPI decomposition/reduction order.

The 16-node run uses the same 12 x 4 x 4 x 4 decomposition as the
calibration reference and reproduces all lambda_100 values exactly.

## Conclusion

The Nv=100 Laplacian calculation is decomposition-independent to
machine precision across the tested 4-, 8-, and 16-node process grids.

4 nodes minimizes resource cost.
8 nodes provides a strong efficiency/turnaround compromise.
16 nodes minimizes single-configuration wall time.

## Plot

The measured wall-time scaling compared with ideal scaling from the
4-node baseline is shown in:

`plots/qX_Nv100_scaling.pdf`
