# Tiny open-boundary clover-insertion result

## Setup and observable

The clover insertion was evaluated on the supplied configuration
`open8x4x4x4n1` using a single MPI rank and all 16 supplied Laplacian modes.
The temporal and spatial parameters were

| Quantity | Value |
|---|---:|
| Lattice | `8 x 4 x 4 x 4` |
| Temporal boundary condition | open |
| Spatial boundary conditions | periodic |
| `t_src`, `t_sink` | `2`, `6` |
| `tau/a` | `4` |
| Temporal links | `t = 2, 3, 4, 5` |
| Insertion time `t_mid` | `4` |
| `R/a` | `2` |
| Mode weights | `rho_i = 1` |
| Number of modes | `Nv = 16` |
| Source-relative probe offset | `(long,tr1,tr2) = (0,1,2)` |

For every source position, the probe offset was rotated cyclically with the
quark-antiquark axis. The calculation was averaged over all 64 spatial source
positions separately for the `x`, `y`, and `z` axes and then over all
`3 x 64 = 192` source-axis combinations.

The production real-part convention was used for the connected estimator,

```text
Delta S = Re(<L S>)/Re(<L>) - <S>.
```

## Clover convention

The field strength and positive plane density were defined by

```text
F_munu = (Q_munu - Q_munu^dagger)/8,
plane_density = Tr(F_munu^dagger F_munu),
E2 = plane(01) + plane(02) + plane(03),
B2 = plane(12) + plane(13) + plane(23),
S  = (E2 + B2)/2.
```

The color trace was removed from each `F_munu`.

## Local clover check

At the absolute diagnostic point `(t,x,y,z) = (4,0,1,2)`, the six positive
plane densities were

| Plane | Density |
|---:|---:|
| `01` | `1.4209974408431977e-01` |
| `02` | `2.2190367299478794e-01` |
| `03` | `3.3388571617157603e-01` |
| `12` | `3.4677403116950495e-01` |
| `13` | `2.1505118531996154e-01` |
| `23` | `2.0268486069921707e-01` |

Their derived values are

| Quantity | Value |
|---:|---:|
| `E2` | `6.9788913325068380e-01` |
| `B2` | `7.6451007718868358e-01` |
| `E2 + B2` | `1.4623992104393673e+00` |
| `S = (E2 + B2)/2` | `7.3119960521968363e-01` |
| Positive ordered-`munu` sum `2(E2+B2)` | `2.9247984208787345e+00` |

The average plaquette of the configuration is
`6.5644455389115142e-01`.

## Connected result

The directional and combined real-part estimators are

| Axis | `Re <L>` | `Delta E2` | `Delta B2` | `Delta S` |
|---:|---:|---:|---:|---:|
| `x` | `3.5177887244109839e-04` | `4.0860396063212567e-03` | `2.9592856597030037e-02` | `1.6839448101676258e-02` |
| `y` | `9.1742676444249427e-04` | `-6.2706931830933677e-03` | `-2.2089677245517247e-02` | `-1.4180185214304974e-02` |
| `z` | `3.1640503093513528e-04` | `-3.4159014584771552e-02` | `1.4508078208785324e-01` | `5.5460883751541235e-02` |
| combined | `5.2853688927290924e-04` | `-9.5380308203456554e-03` | `2.2734944131823198e-02` | `6.5984566557391044e-03` |

For the combined 192-term average, complex quantities are written as
`(Re,Im)`:

```text
<L>      = (5.2853688927290924e-04, 2.2587545260114674e-20),
<S>      =  7.2229780049854819e-01,
<L S>    = (3.8524856035899351e-04, 9.5205216855528757e-06),
Re(<LS>)/Re(<L>) = 7.288962571542873e-01,
Delta S  =  6.5984566557391044e-03.
```

If no real-part projection is applied, the literal single-configuration ratio
`<LS>/<L> - <S>` also has the imaginary part
`1.801297483445276e-02`. The quoted physical estimator is its real part, in
accordance with the production convention.

## Checks and interpretation

The previously validated value of `Re <L>` was reproduced exactly:

```text
expected = 5.2853688927290924e-04,
absolute difference = 0,
L regression check = passed.
```

The maximum pair-conjugacy residual is zero. The largest residual imaginary
part of a directional `<L>` average is `5.8e-20`; the vacuum-`S` averages for
the three rotated axes agree within `3.3e-16`; and the two real-ratio
implementations agree within `1.1e-16`. Also,
`Delta S = (Delta E2 + Delta B2)/2` to floating-point precision.

This is a deterministic one-configuration implementation check. No statistical
uncertainty is assigned, and the directional spread must not be interpreted as
a physical rotational anisotropy. The `L(R,tau)` component has already been
confirmed independently; the clover normalization and connected result remain
to be compared with an independent calculation.

## Detailed output file

- `results/laplace_clover_tiny_open_cfg1_t2-6_r2_Nv16_probe0_1_2.txt`
