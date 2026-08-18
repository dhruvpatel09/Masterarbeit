# Unsummed raw full-clover Laplace-static term

## Scope

This note records one unsummed spatial term on the tiny open-boundary
configuration `open8x4x4x4n1`.  The primary result uses the first spacetime
plane `(1,0)`, with the adjointed static line starting at the absolute
coordinate `(t,x,y,z) = (2,0,0,0)`.  No spatial-origin or axis average is
applied to this term.

The label `spacetime_plane=(1,0)` denotes the orientation of the static-line
pair in the `x` direction.  It is distinct from the clover-plane label `01`
appearing in the decomposition of the action density.

## Parameters and conventions

- Lattice: `Nt x Ns^3 = 8 x 4^3`.
- Configuration: `open8x4x4x4n1` (`cfg = 1`).
- Temporal boundary condition: open; spatial boundary conditions: periodic.
- Laplacian modes: `Nv = 16`, with constant weights `rho_i = 1`.
- Temporal endpoints: `t_src = 2`, `t_sink = 6`, hence `tau/a = 4`.
- Temporal Wilson-line links: `t = 2, 3, 4, 5`.
- Insertion time: `t_mid = 4`.
- Spatial separation: `R/a = 2`.
- Local probe offset from the adjointed line: `(0,1,2)`.
- No HYP smearing is used.

The static-line factor for the first spacetime plane is

```text
L = sum_ij tau_ij(x + R e_x) tau_ij(x)^*,
```

so the line at the spatial origin is adjointed and the line displaced by
`R e_x` is unadjointed.

The raw full-clover action density is

```text
F_munu = (Q_munu - Q_munu^dagger)/8,
S_raw  = -sum_{mu<nu} Tr(F_munu F_munu)
       =  sum_{mu<nu} Tr(F_munu^dagger F_munu).
```

No colour-trace subtraction, additional factor `1/2`, or electric/magnetic
separation is applied.

## Absolute geometry and internal indices

| Object | Absolute coordinate or path | Flat index |
|---|---:|---:|
| Adjointed-line start | `(2,0,0,0)` | spatial site `0` |
| Adjointed-line end | `(6,0,0,0)` | spatial site `0` |
| Unadjointed-line start | `(2,2,0,0)` | spatial site `32` |
| Unadjointed-line end | `(6,2,0,0)` | spatial site `32` |
| Action-density probe | `(4,0,1,2)` | four-dimensional site `262` |

## Static-line value

The unsummed complex static-line factor is

```text
L = +1.3622847961988062e-03
    +1.1536528246516013e-04 i.
```

Interchanging the adjointed and unadjointed spatial sites gives

```text
L_reversed = +1.3622847961988062e-03
             -1.1536528246516013e-04 i.
```

The reversed value is the exact complex conjugate at the printed precision:

```text
L_conjugacy_abs = 0.0000000000000000e+00.
```

## Action-density value

At the absolute probe position `(4,0,1,2)`, the six contributions to
`-Tr(F_munu F_munu)` are:

| Clover plane | Contribution |
|---:|---:|
| `01` | `1.4213589746933711e-01` |
| `02` | `2.2297533246566670e-01` |
| `03` | `3.3604296576313342e-01` |
| `12` | `3.4731939840430276e-01` |
| `13` | `2.1505154551205338e-01` |
| `23` | `2.0472953564018997e-01` |

Their sum is

```text
S_raw = +1.4682546752546832e+00.
```

The difference from the stored point reference is
`2.2204460492503131e-16`, and the largest individual-plane difference is
`1.1102230246251565e-16`; both are below the tolerance `1e-14`.

## Unsummed product

Multiplying the complex static-line factor by the real local action density
gives

```text
L S_raw = +2.0001810210472704e-03
          +1.6938561534154847e-04 i.
```

Both components agree exactly at the printed precision with the product of the
separately printed `L` and `S_raw` values.

## Complete-dump checks

The raw output also contains all `3 x 4^3 = 192` unsummed terms.  A
term-by-term audit gives:

- 64 unique spatial origins for each of the `x`, `y`, and `z` orientations;
- valid absolute endpoints, periodic displacements, rotated probes, and flat
  site indices for every record;
- maximum six-plane summation residual `2.2204460492503131e-16`;
- maximum `L S_raw` multiplication residual zero at the printed precision;
- maximum reversed-line conjugacy residual zero at the printed precision;
- exact reproduction of every printed axis sum from its 64 unsummed records;
- successful clover-point and previously validated `<L>` regressions.

The complete dump is retained for numerical diagnosis.  The primary result of
this note is the single term at the absolute geometry stated above.

Because this calculation uses one gauge configuration, no statistical
uncertainties are assigned and the values are not interpreted physically.

## Reproducibility record

- Slurm job: `156964` (`COMPLETED`).
- Raw output:
  `results/laplace_clover_tiny_open_raw_unsummed_all_sites_cfg1_t2-6_r2_Nv16_probe0_1_2.txt`.
- Raw-output size: `133471` bytes (`218` lines).
- Raw-output SHA-256:
  `d4102dc4598cef5e0c3a636a4ac9bad733e1f6b8ea264abaa5de7e6c445ecc81`.
- Source: `src/test_laplace_clover_tiny_open.c`.
- Compile script: `scripts/compile_test_laplace_clover_tiny_open.sh`.
- Submit script: `scripts/submit_test_laplace_clover_tiny_open.sh`.
