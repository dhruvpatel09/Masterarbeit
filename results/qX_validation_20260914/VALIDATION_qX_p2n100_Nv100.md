# qX p2n100 Nv=100 Laplace eigenvector validation

Date: 2026-09-14

## Ensemble geometry

Global lattice:
192 x 64 x 64 x 64

Process grid:
12 x 4 x 4 x 4

Local lattice:
16 x 16 x 16 x 16

Configuration:
qX_p2n100

## Laplace setup

Nv = 100

Laplace APE:
nAPE = 20
alpha = 0.5

Lanczos:
res              = 1e-9
redTol           = 1e-12
deg              = 10
upper            = 0.15
tolReortho       = 1e-10
maxRestarts      = 30
nLanczosFactor   = 7
redFreqFactor    = 10

## Calibration

All 192 temporal slices converged.

lambda_100 minimum:
4.44931401105784566e-02 at t=43

lambda_100 mean:
4.55467330263308601e-02

lambda_100 maximum:
1.01102449626056037e-01 at t=0

Bulk t=16..175 mean:
4.47961323195665978e-02

Open-boundary values:
t=0   : 1.01102449626056037e-01
t=1   : 5.35372482056311189e-02
t=190 : 5.34733104652426461e-02
t=191 : 1.00164535474082245e-01

## Production eigenvectors

192 / 192 eigenvector HDF5 files written.

Missing timeslices:
0

Duplicate timeslices:
0

Unreadable HDF5 files:
0

Eigenvector dataset:
/disteigvecs = (100, 64, 64, 64, 3)

Datatype:
complex double precision
(r = IEEE F64LE, i = IEEE F64LE)

Each eigenvector HDF5 file:
1258883416 bytes

Total eigenvector directory:
~227 GiB as reported by du

Eigenvalue dataset:
/disteigvals = (192, 100)

Elemental dataset:
/distredelems = (1, 1, 192, 100, 100)

## Independent production/calibration comparison

Compared:
192 / 192 lambda_100 values

Missing:
0

Mean absolute difference:
0.0

Maximum absolute difference:
0.0

The production HDF5 lambda_100 values are bit-for-bit identical,
when decoded as stored little-endian IEEE-754 doubles, to the
canonical calibration result.

## Status

qX_p2n100 Nv=100 eigenvector production VALIDATED.
