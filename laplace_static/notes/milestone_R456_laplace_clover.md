# Milestone: axis-averaged Laplace-clover action-density profile, R/a=4,5,6

## Setup

- Ensemble: Em1
- Lattice: 48 x 24^3
- Ncfg = 100
- Nv = 10 Laplace eigenvectors
- Temporal source averaging: t0 = 0..47
- Spatial source orientations: x, y, z, then axis-averaged
- Observable:
  rho_S = <L S>/<L> - <S>, with S = (E^2 + B^2)/2

## Main result

At tau/a = 4, the axis-averaged rho_S profile shows a robust negative double-dip structure along the static-source axis.

Minima:

- R/a = 4: dx = 1 and 3, i.e. x_mid = -1 and +1
- R/a = 5: dx = 1 and 4, i.e. x_mid = -1.5 and +1.5
- R/a = 6: dx = 1 and 5, i.e. x_mid = -2 and +2

Thus, the strongest negative response remains about one lattice spacing inside each static source, and the structure broadens geometrically with R.

## Tau-stability check

tau/a = 6 is noisier but statistically compatible with tau/a = 4.

Using pyerrors/Gamma-method errors with correlated ensemble labels:

- R/a = 4: max pull inside source interval = 0.688
- R/a = 5: max pull inside source interval = 0.940
- R/a = 6: max pull inside source interval = 1.646

## Interpretation

The clean tau/a=4 profiles should be used as the main spatial results. The tau/a=6 profiles serve as stability diagnostics.

## Key files

- results/laplace_clover_R456_axisavg_advisor_summary.txt
- results/figures/laplace_clover_R4_R5_R6_tau4_rho_S_axis_profile_midpoint_axisavg_jk_ratio_compare.png
- results/figures/laplace_clover_R4tau4_rho_S_heatmap_t0avg48_n100_axisavg_jk_ratio.png
- results/figures/laplace_clover_R5tau4_rho_S_heatmap_t0avg48_n100_axisavg_jk_ratio.png
- results/figures/laplace_clover_R6tau4_rho_S_heatmap_t0avg48_n100_axisavg_jk_ratio.png
