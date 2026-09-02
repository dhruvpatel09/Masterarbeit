# Weighted Laplace-eigenmode profile analysis

## Scope

This directory summarizes the fixed-source-time profile-selection study for
the Laplace-static correlator.  The input is the 100-configuration,
mode-resolved `PROFILE_BASIS=delta` scan from job `156696`, with

$$
N_v=10,\qquad t_0/a=0,\qquad R/a=1,\ldots,6,\qquad
\tau/a=1,\ldots,8.
$$

The calculation uses one on-axis, $x$-oriented separation and averages over
all $24^3=13\,824$ spatial source positions.  It is an offline
profile-selection analysis, not the final temporal-source- and axis-averaged
production calculation.

For a Gaussian Laplacian-eigenmode profile,

$$
\rho_i(\sigma)=\exp\!\left[-\frac{\lambda_i^2}{4\sigma^2}\right],
\qquad
w_i(\sigma)=\rho_i^2
=\exp\!\left[-\frac{\lambda_i^2}{2\sigma^2}\right].
$$

The constant profile is the $\sigma\to\infty$ limit, $w_i=1$.

## Validation

- All 100 raw files and their mode-resolved matrix dimensions were validated.
- The constant profile reconstructed the original `DATA` channel with a
  maximum absolute residual of $2.038\times10^{-21}$.
- The $\sigma=0.05$ diagonal agreed with the earlier seven-profile
  projection implementation.
- Delete-one jackknife samples are propagated through scalar fits, GEVP
  solutions, and fixed-vector projections.
- A GEVP estimate is treated as reportable only when the full-sample metric
  and every delete-one metric are accepted, the relevant principal
  correlators are positive, and the jackknife error is finite.

## Direct scalar profiles

Correlated constant fits use the established window
$\tau/a=3,\ldots,5$.

| $R/a$ | Constant profile | $\sigma=0.05$ | Paired shift |
|---:|---:|---:|---:|
| 1 | $0.481552(1266)$ | $0.478814(2501)$ | $-1.40\sigma$ |
| 2 | $0.696894(3998)$ | $0.694269(4773)$ | $-1.05\sigma$ |
| 3 | $0.834223(7339)$ | $0.832257(7714)$ | $-0.75\sigma$ |
| 4 | $0.924332(11326)$ | $0.922736(11800)$ | $-0.46\sigma$ |
| 5 | $1.026585(20660)$ | $1.031387(20422)$ | $+1.08\sigma$ |
| 6 | $1.116271(27947)$ | $1.120689(26809)$ | $+0.58\sigma$ |

The $\sigma=0.05$ profile is statistically compatible with the constant
profile at every radius.  Its fit uncertainty ranges from $0.96$ to
$1.98$ times the constant-profile uncertainty and does not provide a
consistent precision or plateau improvement.  The broader Gaussian profiles
approach the constant profile and likewise provide no material improvement.

For $\sigma=0.01$, the endpoint weights span approximately
$1.15\times10^{-15}$ to $8.01\times10^{-5}$ over the retained modes.  The
resulting fit uncertainties are $7.7$ to $131$ times larger than the
constant-profile uncertainties.  The $R/a=6$ fit is invalid, and the
$R/a=4$ fit has $p=0.014$.  This profile is therefore unsuitable as a
direct production channel for the present $N_v=10$ data.

## Full three-profile GEVP

The basis is

$$
\{\sigma=0.01,\ \sigma=0.05,\ \sigma=\infty\}.
$$

| $\tau_{\rm ref}/a$ | Metric result | Delete-one support | Consequence |
|---:|---|---|---|
| 1 | Condition number $3.66\times10^3$ to $1.20\times10^4$ | $99/100$ at every radius | No complete jackknife errors |
| 2 | Condition number $6.17\times10^4$ to $2.99\times10^5$ | $100,100,100,87,97,73$ out of 100 | Only $R/a=1,2,3$ have complete metrics |
| 3 | Negative minimum metric eigenvalue at every radius | $0,0,0,0,0,1$ out of 100 | Full 3x3 GEVP invalid |

The third principal channel is non-positive or undefined throughout the
relevant range.  The second channel has no reference-time-stable plateau with
complete uncertainties.  Neither channel supports an excited-static-potential
interpretation.

The diagnostic-exclusion count in the text report counts repeated failed
evaluations over radii, times, and delete-one samples.  It is not a count of
independent physical failures.

## Reduced two-profile gate

At $\tau_{\rm ref}/a=2$, the normalized correlation between the
$\sigma=0.05$ and constant operators is $0.999616$ to $0.999854$.  Their
near duplication produces the smallest direction of the full metric.

Two non-duplicate bases are therefore evaluated separately:

$$
\{\sigma=0.01,\sigma=0.05\},\qquad
\{\sigma=0.01,\sigma=\infty\}.
$$

Their $\tau_{\rm ref}/a=2$ central-metric condition numbers are approximately
$22.5$ to $26.3$ and $20.5$ to $23.0$, respectively.  The updated
analysis script evaluates both bases at $\tau_{\rm ref}/a=1,2,3$, including
all delete-one samples and fixed vectors defined by
$(\tau_{\rm ref}/a,\tau_d/a)=(2,3)$.

The reduced-basis branch is accepted only if the candidate $n=2$ energy has

1. a positive and accepted metric in all delete-one samples;
2. positive adjacent principal correlators and a finite jackknife error;
3. stability across both reduced bases and admissible reference times; and
4. an effective-energy plateau.

The reference-time regime follows $t_0\ge t/2$.  For an adjacent-time ratio
$\lambda(\tau)/\lambda(\tau+a)$, the output column
`reference_time_condition_met` records
$\tau_{\rm ref}\ge(\tau+a)/2$.

### Reduced-basis outcome

Both reduced metrics are positive in the full sample and all 100 delete-one
samples at $\tau_{\rm ref}/a=1,2$.  At $\tau_{\rm ref}/a=2$, their condition
numbers are $22.5$--$26.3$ for $\{0.01,0.05\}$ and $20.5$--$23.0$ for
$\{0.01,\infty\}$.  The near-singularity of the full three-profile metric is
therefore removed at this reference time.

The later reference time does not remain admissible.  At
$\tau_{\rm ref}/a=3$, neither reduced basis has complete delete-one metric
support at any radius.  The support ranges from $1/100$ to $96/100$ for
$\{0.01,0.05\}$ and from $2/100$ to $98/100$ for
$\{0.01,\infty\}$; several central metrics are already non-positive.

At $(\tau_{\rm ref}/a,\tau/a)=(2,3)$, the ground level agrees centrally with
the constant-profile plateau at every radius.  Its uncertainty is nevertheless
$16.8$ to $224$ times the corresponding constant-profile fit uncertainty.
For the candidate second level, neither basis has a single finite result with
$100/100$ jackknife support anywhere in the reference-time regime.  The
fixed-vector candidate also has no plateau: over $\tau/a=2,\ldots,5$, most
complete estimates are negative or have order-one-to-several-unit
uncertainties.  The only positive complete estimates are isolated
$R/a=3$, $\tau/a=5$ points,
$0.591\pm0.521$ and $0.643\pm0.446$ for the two bases.

The reduced-basis gate is therefore not passed.  The present data do not
support a first-excited static potential from these Gaussian operators.

## Current production choice

The supported source profile remains

$$
N_v=10,\qquad \rho_i=1.
$$

The scalar comparison and both full and reduced variational diagnostics lead
to the same production choice.  No additional weighted production run is
indicated by this study.

## Output files

The updated `analyze_laplace_gaussian_gevp.py` writes the original scalar and
full-basis outputs together with

- `gaussian_gevp_156696_reduced_conditioning.csv`;
- `gaussian_gevp_156696_reduced_principal.csv`;
- `gaussian_gevp_156696_reduced_projected.csv`;
- `gaussian_gevp_156696_reduced_vectors.csv`;
- `gaussian_gevp_156696_reduced_levels_tref2.png`; and
- `gaussian_gevp_156696_reduced_projected_tref2_td3.png`.

Hollow plot markers denote finite central estimates with incomplete
delete-one jackknife support.  Rows without finite central estimates are not
plotted.

## Reproduction

```bash
cd /home/m2130292/Masterarbeit

python3 -B laplace_static/scripts/analyze_laplace_gaussian_gevp.py --self-test

delta100_dir="laplace_static/results/weighted_pilot/delta_n1-100_R0-6_T1-8_t0fixed0"
gevp_dir="laplace_static/results/weighted_profile_gevp"

python3 -B laplace_static/scripts/analyze_laplace_gaussian_gevp.py \
  "$delta100_dir" 156696 \
  --eigenvalue-root mental/runs_Em1p4_Nv10_qcdnew_full \
  --output-dir "$gevp_dir" \
  --gevp-reference-times 1 2 3 \
  --plot-reference-time 2 \
  --optimization-time 3
```

No new gauge-field or eigenvector scan is required.

## Excluded scope

This analysis does not include a Gaussian-weighted temporal-source-averaged
production run, a weighted clover/action-density calculation, or a
$d_z$-resolved three-dimensional clover profile.

## Provenance

The source snapshot used for the original full-basis output has SHA-256
`b817b5d2920779ad6a9334dd6240ed89938012946f442377ee28823f82e812b1` for
`analyze_laplace_gaussian_gevp.py`.  The full- plus reduced-basis version used
for the reported numerical output has SHA-256
`255acbf27b28c2d5a9300ed26e6ff13b0a9d47ae7ac30d354353f1989c27151c`.
The final LF-only CSV-output revision has SHA-256
`d694b3dce08774e91f6349c8e9135e284afe30e2f3accbcfb1c30698e66f6e3c`;
this revision changes line endings only.

The plateau-comparison source snapshot has SHA-256
`8b48e0e4f5896cb82c360b04927b9f60b69aa69fdd54908c7f346b0612078bd5`;
the neutral-text revision has SHA-256
`f990e6cd63acbe282935642f1b8b4807cc9aed97096960761f7ec2bd8f288ce4`,
and the LF-only CSV-output revision has SHA-256
`0c8e047c6714b54192af78fe364596222f3e164e663e32a5d63864e63a76f7c6`.

## References

- R. Höllwieser et al., *Constructing static quark-antiquark creation
  operators from Laplacian eigenmodes*, arXiv:2212.08485.
- B. Blossier et al., *On the generalized eigenvalue method for energies and
  matrix elements in lattice field theory*, arXiv:0902.1265.
