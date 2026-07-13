# Tiny open-boundary validation of `L(R, tau)`

## Setup

The Laplace-static observable was evaluated on the provided gauge
configuration `open8x4x4x4n1` using a single MPI rank and the corresponding
HDF5 eigenvectors.

| Quantity | Value |
|---|---:|
| Lattice | `8 x 4 x 4 x 4` |
| Temporal boundary condition | open |
| Spatial boundary conditions | periodic |
| Source time `t_src` | `2` |
| Sink time `t_sink` | `6` |
| Temporal separation `tau/a` | `4` |
| Temporal links | `t = 2, 3, 4, 5` |
| Spatial separation `R/a` | `2` |
| Mode weights | `rho_i = 1` |
| Numbers of modes | `Nv = 10, 16` |

For each spatial direction, the observable was averaged over all `4^3 = 64`
source positions. The combined result is the average over the three axial
directions and therefore contains `3 x 64 = 192` terms.

## Results

The real parts of the directional and combined averages are

| `Nv` | `Re L_x` | `Re L_y` | `Re L_z` | `Re L_axisavg` |
|---:|---:|---:|---:|---:|
| 10 | `3.8752283161748148e-04` | `7.5128727754741735e-04` | `4.8028826296859421e-04` | `5.3969945737783107e-04` |
| 16 | `3.5177887244109839e-04` | `9.1742676444249427e-04` | `3.1640503093513528e-04` | `5.2853688927290924e-04` |

The largest residual imaginary part of a directional average is
`6.8e-21` for `Nv=10` and `5.8e-20` for `Nv=16`. These values are consistent
with zero at floating-point precision. The maximum pair-conjugacy residual is
zero for both mode counts.

Because this calculation uses only one gauge configuration, no statistical
uncertainty is assigned.

## Detailed output files

- `results/laplace_tiny_open_cfg1_t2-6_r2_Nv10.txt`
- `results/laplace_tiny_open_cfg1_t2-6_r2_Nv16.txt`

The `PAIR` records give the correlator for a specified source position and axial
displacement. The `AXIS` records give the averages over all 64 source positions
for one direction, and `COMBINED` gives the 192-term average over all three directions.
