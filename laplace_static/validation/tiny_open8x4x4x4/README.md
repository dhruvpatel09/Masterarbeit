# Tiny open-boundary validation

This validation checks the Laplacian-eigenmode static-source correlator on the
provided `open8x4x4x4n1` gauge configuration before further production runs.

## Geometry and observable

- lattice: `Nt x Ns^3 = 8 x 4^3`;
- temporal boundary condition: open;
- spatial boundary conditions: periodic;
- temporal endpoints: `t_src = 2`, `t_sink = 6`;
- temporal extent: `tau/a = t_sink - t_src = 4` links;
- spatial separation: `r/a = 2`;
- constant eigenmode weights: `rho_i = 1`;
- execution: exactly one MPI rank.

For each spatial site `x`, the program constructs

```text
P(x;2,6) = U0(x,2) U0(x,3) U0(x,4) U0(x,5),
tau_ij(x;2,6) = v_i(x,2)^dagger P(x;2,6) v_j(x,6),
L(x,y) = sum_ij tau_ij(y;2,6) tau_ij(x;2,6)^*.
```

There is no temporal modulo operation and no temporal-source averaging.  The
program averages over all 64 spatial source sites for each of the three axial
displacements and also reports the combined 192-term axis average.

## Local data

The ignored `data/` directory contains the selected files extracted from
`tiny_dist.zip`:

- `cnfg/open8x4x4x4n1`;
- HDF5 and DAT eigenvectors for `t=0,...,7`;
- the HDF5 eigenvalues;
- the supplied MATLAB eigenpairs.

The original archive used for this validation has SHA-256

```text
2c563e5ca1341b86800ccfc47b960f4bd85d0d9a0897de7492508c7db10a8aed
```

The HDF5 files contain 16 modes.  The default validation job evaluates both
`Nv=10`, matching production, and `Nv=16`, using every supplied HDF5 mode.

The tiny files store `/time` as an unsigned 32-bit integer.  The production
`DistEigvecsHdf5Reader` uses a signed `int` output and performs strict HDF type
equality before reading, so it rejects these otherwise compatible files.  The
validation executable therefore uses a small serial HDF5 reader local to the
test.  HDF5 performs the safe integer conversion, while the code still checks
the configuration ID, time, run name, spatial dimensions, available modes,
compound-complex layout, and full eigenvector dataset shape.  Neither the
shared production reader nor the provided HDF5 files are modified.

## Build and submit

From `/home/m2130292/Masterarbeit/laplace_static`:

```bash
bash scripts/compile_test_laplace_tiny_open.sh
sbatch scripts/submit_test_laplace_tiny_open.sh
```

The job writes one text result for each requested mode count under `results/`.
Override the defaults at submission when needed, for example:

```bash
sbatch --export=ALL,NVECS_LIST="10",PRINT_SITE_VALUES=1 \
  scripts/submit_test_laplace_tiny_open.sh
```

## Output and checks

The compact output contains:

- `META`, `INPUT`, and `PATH` records documenting the calculation;
- `PAIR` for the source point `(0,0,0)` in each direction;
- `AXIS` for the 64-source average in each direction;
- `COMBINED` for the 192-term average;
- `CHECK` for pair-conjugacy and residual imaginary parts.

Since `r/a=2` is half the spatial extent, reversing the displacement maps each
pair to its conjugate partner.  A complete volume average should therefore be
real up to floating-point roundoff.  This is a structural diagnostic, not a
substitute for comparison with an independent DAT/MATLAB calculation.

The result is a deterministic one-configuration implementation check.  It is
not an ensemble estimate and requires no jackknife, autocorrelation, or fit
analysis.

## Independent DAT reference

`scripts/reference_laplace_tiny_open_dat.py` independently reconstructs the
openQCD gauge field and reads the raw DAT eigenvectors without qcdlib or the
HDF5 reader.  It evaluates the same matrix formula with NumPy:

```bash
python3 -c 'import numpy; print(numpy.__version__)'

python3 scripts/reference_laplace_tiny_open_dat.py --nvecs 10 \
  | tee validation/tiny_open8x4x4x4/results/reference_dat_Nv10.txt

python3 scripts/reference_laplace_tiny_open_dat.py --nvecs 16 \
  | tee validation/tiny_open8x4x4x4/results/reference_dat_Nv16.txt
```

The C/HDF5 `PAIR`, `AXIS`, and `COMBINED` values must agree with the Python/DAT
`REF_PAIR`, `REF_AXIS`, and `REF_COMBINED` values to floating-point precision.
Agreement tests the observable independently of both the HDF5 reader and the
qcdlib openQCD input routine.

## Validated numerical result

The validation completed successfully on 2026-07-13.  The entries below are
the real parts of the 64-source spatial averages in each axial direction and
their combined 192-term average:

| `Nv` | `Re L_x` | `Re L_y` | `Re L_z` | `Re L_axisavg` |
|---:|---:|---:|---:|---:|
| 10 | `3.8752283161748148e-04` | `7.5128727754741735e-04` | `4.8028826296859421e-04` | `5.3969945737783107e-04` |
| 16 | `3.5177887244109839e-04` | `9.1742676444249427e-04` | `3.1640503093513528e-04` | `5.2853688927290924e-04` |

The C/HDF5 and Python/DAT calculations agree to floating-point precision.  The
largest absolute differences over all reported `PAIR`, `AXIS`, and `COMBINED`
values are `4.4e-19` for `Nv=10` and `1.4e-18` for `Nv=16`; the corresponding
largest relative differences are `6.0e-16` and `8.2e-16`.  Pair-conjugacy
residuals vanish, and the largest residual imaginary part of an axial spatial
average is `5.8e-20`.

Changing from `Nv=10` to all 16 supplied modes changes the combined real value
by `-2.0683%`.  This is a useful implementation diagnostic, not a production
mode-truncation systematic.  Likewise, the directional spread on one gauge
configuration must not be interpreted as physical rotational anisotropy.
