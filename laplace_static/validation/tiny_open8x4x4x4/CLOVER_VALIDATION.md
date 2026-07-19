# Tiny open-boundary validation of the clover insertion

## Purpose

This test extends the advisor-verified tiny-lattice calculation of
`L(R, tau)` to the connected clover insertion

```text
Delta S = <L S>/<L> - <S>.
```

The existing `test_laplace_tiny_open.c` validation is not modified. The new
executable first reproduces its full-mode `L(R, tau)` result and then evaluates
the local clover density.

## Default geometry

| Quantity | Value |
|---|---:|
| Configuration | `open8x4x4x4n1` |
| Temporal boundary condition | open |
| Spatial boundary conditions | periodic |
| `t_src`, `t_sink` | `2`, `6` |
| `tau/a` | `4` |
| Temporal links | `t = 2, 3, 4, 5` |
| Insertion time `t_mid` | `4` |
| `R/a` | `2` |
| Mode weights | `rho_i = 1` |
| Number of modes | `Nv = 16` |
| Absolute diagnostic point `(t,x,y,z)` | `(4,0,1,2)` |
| Source-relative probe offset | `(long,tr1,tr2) = (0,1,2)` |

For every source position, the source-relative probe offset is rotated
cyclically with the quark-antiquark axis. Results are averaged over all 64
spatial sources separately for the `x`, `y`, and `z` axes and then over all
`3 x 64 = 192` source-axis combinations.

A fixed absolute value of `S(z)` must not be multiplied by every translated
`L` contribution when forming this single-configuration estimator: it would
factor out and make `Delta S` identically zero. The translated probe position
preserves the displacement of the insertion relative to each source pair,
matching the production-probe construction.

## Clover convention

The code uses the same convention as the production clover probe:

```text
F_munu = (Q_munu - Q_munu^dagger)/8,
plane_density = Tr(F_munu^dagger F_munu),
E2 = plane(01) + plane(02) + plane(03),
B2 = plane(12) + plane(13) + plane(23),
S  = (E2 + B2)/2.
```

The color trace is removed from each `F_munu`. The result file prints all six
plane densities and the derived quantities at the absolute diagnostic point,
so signs and factors of two can be compared directly with an independent
calculation.

For consistency with the production analysis, the reported estimator is

```text
delta_S = Re(<L S>)/Re(<L>) - <S>.
```

The file also reports `delta_S_complex = Re(<L S>/<L>) - <S>` as a diagnostic.
These agree when the residual imaginary part of the averaged `L` is negligible.

## Build and run

From `/home/m2130292/Masterarbeit/laplace_static`:

```bash
bash scripts/compile_test_laplace_clover_tiny_open.sh
jobid=$(sbatch --parsable scripts/submit_test_laplace_clover_tiny_open.sh)
echo "Submitted job ${jobid}"
```

The default output is

```text
validation/tiny_open8x4x4x4/results/
laplace_clover_tiny_open_cfg1_t2-6_r2_Nv16_probe0_1_2.txt
```

Set `PRINT_SITE_VALUES=1` when submitting to retain all 192 individual
source-axis contributions. The default concise output is sufficient for the
primary validation.

## Acceptance checks

The run is accepted only if the combined real part of `L` reproduces the
advisor-verified value

```text
5.2853688927290924e-04
```

within the embedded absolute tolerance `1e-15`. The `CHECK` record also reports
the largest pair-conjugacy residual, residual imaginary part of the directional
`L` averages, spread of the vacuum `S` averages among axes, and difference
between the two ratio conventions.

This remains a one-configuration implementation validation. It does not carry
a statistical uncertainty and must not be interpreted as an ensemble flux-tube
measurement.
