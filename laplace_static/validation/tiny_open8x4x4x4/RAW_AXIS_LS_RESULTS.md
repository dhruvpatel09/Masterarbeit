# Raw full-clover axis-resolved `<LS>` cross-check

## Purpose

This note records the one-configuration raw full-clover cross-check for the
tiny open-boundary configuration `open8x4x4x4n1`.  The primary outputs are the
three separate spatial averages of `<LS>` for the spacetime planes `(1,0)`,
`(2,0)`, and `(3,0)`.  Each value is an average over 64 spatial origins; the
three directions are recorded separately and are not combined in the
axis-resolved result.

## Parameters and conventions

- Lattice: `Nt x Ns^3 = 8 x 4^3`.
- Configuration: `open8x4x4x4n1` (`cfg = 1`).
- Temporal boundary condition: open; spatial boundary conditions: periodic.
- Laplacian modes: `Nv = 16`, with constant weights `rho_i = 1`.
- Temporal endpoints: `t_src = 2`, `t_sink = 6`, hence `tau/a = 4`.
- Temporal Wilson-line links: `t = 2, 3, 4, 5` (no temporal wrapping).
- Insertion time: `t_mid = 4`.
- Spatial separation: `R/a = 2`.
- Local origin line `(0,0,0)`: adjointed,
  `[prod_{t=2}^{5} U_0(t,x)]^dagger`.
- Local displaced line `(2,0,0)`: unadjointed.
- Local probe offset from the adjointed line: `(0,1,2)`.
- No HYP smearing is used in this cross-check.

The clover field and full local action density are

```text
F_munu = (Q_munu - Q_munu^dagger)/8,
S_raw  = -sum_{mu<nu} Tr(F_munu F_munu)
       =  sum_{mu<nu} Tr(F_munu^dagger F_munu).
```

No colour-trace subtraction is applied to `F_munu`, there is no additional
factor `1/2`, and no electric/magnetic separation is made.

The cyclic rotation of the local geometry gives:

| Spacetime plane | Axis | Global displaced-line offset | Global probe offset |
|---|---:|---:|---:|
| `(1,0)` | x | `(2,0,0)` | `(0,1,2)` |
| `(2,0)` | y | `(0,2,0)` | `(2,0,1)` |
| `(3,0)` | z | `(0,0,2)` | `(1,2,0)` |

## Absolute-position clover check

At the absolute position `(t,x,y,z) = (4,0,1,2)`, the six contributions to
`-Tr(F_munu F_munu)` are:

| Plane | Contribution |
|---:|---:|
| `01` | `1.4213589746933711e-01` |
| `02` | `2.2297533246566670e-01` |
| `03` | `3.3604296576313342e-01` |
| `12` | `3.4731939840430276e-01` |
| `13` | `2.1505154551205338e-01` |
| `23` | `2.0472953564018997e-01` |

Their sum is

```text
S_raw(4,0,1,2) = 1.4682546752546832e+00.
```

The difference from the reference value `1.4682546752546830e+00` is
`2.2204460492503131e-16`; the largest individual-plane difference is
`1.1102230246251565e-16`.  Both are below the regression tolerance `1e-14`.

## Axis-resolved results

For each row,

```text
<LS>_axis = (1/64) sum_over_spatial_origins L_axis S_raw,
```

with the origin line adjointed and the displaced line unadjointed.

| Spacetime plane | Axis | `Nsrc` | `Re <L>` | `Im <L>` | `Re <LS>` | `Im <LS>` |
|---|---:|---:|---:|---:|---:|---:|
| `(1,0)` | x | 64 | `3.5177887244109839e-04` | `5.7598240413292423e-20` | `5.2123931442376432e-04` | `-1.3283571590837837e-05` |
| `(2,0)` | y | 64 | `9.1742676444249427e-04` | `1.0164395367051604e-20` | `1.3031256847886267e-03` | `4.9513870568180714e-05` |
| `(3,0)` | z | 64 | `3.1640503093513528e-04` | `0.0000000000000000e+00` | `4.9346521324629140e-04` | `2.1346564030651379e-05` |

The vacuum spatial average is the same for all three rotated probe offsets,
up to roundoff:

```text
<S_raw> = 1.4488716275970441e+00,
maximum directional spread = 2.2204460492503131e-16.
```

## Additional diagnostics

The combined 192-term average is retained only as an internal diagnostic; it
is not a replacement for the three separate 64-term averages:

```text
<L>_combined  = 5.2853688927290924e-04
                + 2.2587545260114674e-20 i,
<LS>_combined = 7.7261007081956079e-04
                + 1.9192287669331419e-05 i.
```

The combined real part of `<L>` reproduces the previously validated reference
exactly at the printed precision.  The maximum conjugacy residual of a static
pair is zero at the printed precision, and the maximum absolute imaginary part
of an axis-averaged `<L>` is `5.7598240413292423e-20`.

The nonzero imaginary parts of the axis-resolved `<LS>` values are not a failed
reality check.  At `R = Ns/2`, translating the spatial origin by `R` conjugates
the static pair, which makes the complete `<L>` average real.  With `S_raw`
inserted at a fixed asymmetric offset relative to the designated adjointed
line, the translated partner samples `S_raw` at a different relative point;
therefore `<LS>` is not required to be real on a single configuration.

Because this calculation uses only one gauge configuration, no statistical
uncertainties are assigned and the directional differences are not interpreted
physically.

## Reproducibility record

- Slurm job: `156963` (`COMPLETED`, exit code `0:0`, elapsed `00:00:02`).
- Raw output:
  `validation/tiny_open8x4x4x4/results/laplace_clover_tiny_open_raw_axisLS_cfg1_t2-6_r2_Nv16_probe0_1_2.txt`.
- Source: `src/test_laplace_clover_tiny_open.c`.
- Compile script: `scripts/compile_test_laplace_clover_tiny_open.sh`.
- Submit script: `scripts/submit_test_laplace_clover_tiny_open.sh`.

This result supersedes the earlier traceless, half-normalized clover-insertion
result for this normalization-and-orientation cross-check.
