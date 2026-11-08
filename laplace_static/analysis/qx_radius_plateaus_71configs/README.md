# qX radius and plateau diagnostics

Offline analysis of the existing 71-configuration static-correlator measurements: complete measured-time plots at each distance, correlated start/stop window scans, focused checks near \(R/a=24\), and leading exponential approach diagnostics. No new MPI measurements are introduced.

Dependencies are pinned at [2c879d6](https://github.com/dhruvpatel09/Masterarbeit/tree/2c879d63fcd0b2dd4aa413ecd5324aa75add126d). Labels are 100 through 800 in steps of 10, including n640. The statistical sample is 71 configurations after source and axis averaging; sites and axes are not extra independent configurations.

## Coverage

The lattice is \(192\times64^3\), with open temporal and periodic spatial boundaries.

| Quantity | Coverage |
|---|---|
| Spatial half-extent | \(R/a=32\) |
| Temporal half-extent | 96 |
| Dense scan, \(N_v=100\) | Every integer \(R/a=1,\ldots,32\); correlator times 0–16 |
| Extended scan, \(N_v=40,60,80,100\) | \(R/a=1,2,4,8,12,16,20,23,24,25,28,32\); correlator times 0–32 |
| Source time and projection | \(t_0/a=90,\ p_{\rm iter}=10\) |

Time 32 is the measured cutoff, not half the temporal lattice. Correlators through 32 give effective energies with lower endpoint 0–31, plotted at midpoints 0.5–31.5. Longer-time measurements at the other integer radii are unavailable.

There are 68 measured \((N_v,R)\) series, 1,924 correlator means and 1,856 possible effective energies. The 42 ratios with a nonpositive ensemble-mean endpoint remain explicitly undefined.

## Statistical method

The independent-configuration baseline follows the ensemble assessment supplied by the advisors; this is not a new thermalization/autocorrelation determination.

\[
aE_{\rm eff}(n)=\log[\overline C(n)/\overline C(n+1)],\qquad
\delta E_i(n)=\frac{C_i(n)-\overline C(n)}{\overline C(n)}
-\frac{C_i(n+1)-\overline C(n+1)}{\overline C(n+1)}.
\]

Covariances use the sample covariance of these paired fluctuations divided by 71 (sample denominator 70). Logs follow ensemble averaging, with all signed configuration values retained. Normalized-correlator errors retain covariance with the equal-time denominator. No positive-configuration filtering, absolute-value logarithms or replacement of undefined ratios by zero is used.

Constant fits use full time covariance and a strictly positive definite Cholesky solve. There is no pseudoinverse, ridge, diagonal fallback or chi-square error inflation. Delete-one checks hold full-sample GLS weights fixed, matching previous work. A delete-one error is blank if a necessary leave-one-out mean is nonpositive. Window, cutoff and axis differences preserve configuration pairing.

The catalogue retains all 2,068 prescribed attempts: 1,957 defined, 111 with nonpositive endpoint means. The exploratory display screen requires minimum endpoint mean/SE at least 3, approximate Gaussian chi-square tail at least 0.05, and covariance condition number at most \(10^{12}\). It is not a plateau acceptance rule, a multiple-comparison-adjusted test, or proof of negligible excited-state contamination. Estimated covariance and nonlinear weak-signal ratios limit the interpretation of the p-values.

## Findings near R/a = 24

At correlator time 24, signed mean/SE is 14.78, 9.23, 3.23 and 0.50 for \(N_v=40,60,80,100\), respectively. More modes do not uniformly improve late-time signal. Mode cutoff changes the source operator; piter controls numerical SU(3) projection. \(N_v=100\) is not an exact physical reference.

Representative constant-fit diagnostics:

| \(N_v\) | Ratio window | \(aE\), independent SE | Minimum endpoint mean/SE | Approximate p |
|---|---|---|---:|---:|
| 40 | 6–11 | \(0.26982\pm0.00219\) | 67.93 | 0.492 |
| 40 | 14–19 | \(0.26564\pm0.00512\) | 25.28 | 0.080 |
| 40 | 22–27 | \(0.27953\pm0.02327\) | 4.24 | 0.025 |
| 40 | 8–31 | \(0.26271\pm0.00174\) | 1.67 | \(1.6\times10^{-6}\) |
| 100 | 6–11 | \(0.47489\pm0.01434\) | 8.44 | 0.272 |
| 100 | 14–19 | \(0.21587\pm0.09777\) | 0.23 | 0.878 |
| 100 | 22–27 or 8–31 | Undefined | Negative mean | — |

The \(N_v=100\) 14–19 fit demonstrates why a high p-value alone is insufficient: the signal is weak and its paired delete-one error is undefined. Extending the \(N_v=40\) fit through 31 gives a small formal error but poor constant-model agreement and weak endpoints.

For \(N_v=40\), the paired change from 6–11 to 14–19 is \(-0.00418\pm0.00530\), compatible with zero. This does not validate all later windows: the 22–27 fit has poor approximate constant-model agreement. Compare neighboring distances 20, 23, 25 and 28 alongside the signal plots.

At \(R/a=12\), window 14–19, the paired \(N_v=40\) minus \(100\) difference reproduces \(-0.0016367\pm0.0004814\), or \(-3.40\) pointwise SE. Source/cutoff dependence calls for excited-state and overlap studies, not a claim of different physical potentials.

## Leading approach diagnostic

\[
aE_{\rm eff}(n)=A+B\exp[-c(n-n_{\min})].
\]

This is an asymptotic approach ansatz, not the exact effective energy of a two-exponential correlator. Effective energies contain no independent information beyond the correlators. Signed \(A,B\) avoid imposing a positive spectral representation on these source diagnostics.

At 65 logarithmically spaced rates \(0.005\le c\le2\), conditional linear solves determine \(A,B\); every grid local minimum is refined. Covariance is held fixed. The 164 attempts contain 134 defined and 30 undefined fits; 17 pass the additional exploratory signal, model, rate/amplitude, curvature and stationarity checks. These are not 17 accepted ground-state determinations.

Parameter errors use the optimum's response to observations, including residual terms in the observed Hessian. They are conditional on covariance/model/window, omit covariance-estimation, truncation-model and selection uncertainty, and remain blank for nonpositive local curvature. Boundary rates and unresolved amplitudes/rates remain flagged.

## Source shape and spatial interpretation

The equal-time plot includes the source-shape minimum near sampled \(R/a=25\) for \(N_v=100\). Periodic spatial half-box symmetry, temporal flattening, a source-shape minimum and physical string breaking are distinct. The spatial-window comparison assembles candidate energies without selecting \(E_0(R)\), fitting a new Cornell model or inferring a string tension.

## Figures

Each has PNG and SVG versions.

| Plot | Purpose |
|---|---|
| [Coverage](results/measurement_coverage.png) | Unmeasured geometry versus undefined observations |
| [R1–R8](results/effective_energy_dense_R01_R08.png), [R9–R16](results/effective_energy_dense_R09_R16.png), [R17–R24](results/effective_energy_dense_R17_R24.png), [R25–R32](results/effective_energy_dense_R25_R32.png) | Every integer radius in the common short scan |
| [Extended 1](results/effective_energy_extended_page1.png), [2](results/effective_energy_extended_page2.png), [3](results/effective_energy_extended_page3.png) | All four cutoffs, all 12 radii, complete time range and detail panels |
| [Signed neighboring signals](results/signed_signal_R20_R28.png), [R24 audit](results/R24_signal_audit.png) | Signed means, normalization and signal/error |
| [Window scan 1](results/constant_window_scan_page1.png), [2](results/constant_window_scan_page2.png) | Six-ratio windows at every extended radius |
| [Endpoint sensitivity](results/window_endpoint_sensitivity.png) | Vary stopping time at a fixed start |
| [Paired changes](results/paired_window_changes.png) | Compare correlated windows |
| [Spatial comparison](results/spatial_window_comparison.png) | Candidate energies versus radius |
| [Half-box source shape](results/half_box_source_shape.png) | Equal-time behavior near spatial half-extent |
| [Approach checks](results/approach_model_checks.png), [asymptote sensitivity](results/approach_asymptote_sensitivity.png) | Approximate model and conditional errors |
| [Paired axes](results/R24_paired_axis_checks.png) | Direction differences at R24 |

Hollow energy markers denote endpoint signal/SE below three. Crosses on a separate validity strip denote undefined ratios, not zero energy. Filled/hollow fit markers denote the exploratory screen. Detail panels clip error bars; edge triangles identify off-scale means. Full-range panels retain those values. Cutoff curves can overlap when their values agree.

## Tables

- correlator_means.csv, effective_energies.csv, data_coverage.csv preserve primary decimal strings and explicit coverage/validity.
- constant_fit_catalogue.csv contains every start/stop attempt, errors, fit statistics and flags.
- paired_window_contrasts.csv and paired_cutoff_contrasts.csv preserve correlations.
- axes_near_R24.csv includes all pairwise directions at the five focused distances.
- spatial_window_estimates.csv explicitly labels outside_measured_time.
- approach_fits.csv, approach_solutions.json retain model attempts, candidates, errors and flags.
- SUMMARY.json records cohort, coverage and regression totals.

## Reproduction

From the repository root on Stromboli:

    qx_analysis=laplace_static/analysis/qx_radius_plateaus_71configs
    python3 "$qx_analysis/verify_files.py"
    python3 "$qx_analysis/test_analysis.py"
    python3 "$qx_analysis/reproduce.py" --verify

Analysis requires NumPy, SciPy and Matplotlib (tested versions in inputs/ENVIRONMENT.json). The integrity check uses only the standard library and runs on SuperMUC. Source syntax is compatible with Python 3.6.

Replay verifies exact bytes of two raw tables and eight numerical dependencies; requires all 71 labels at every measured geometry; checks 32,589 retained raw overlap records; reconstructs means/energies/errors; audits 204,800 published combined-axis covariance entries; and reproduces 544 previous short-scan fits and 576 extended fit attempts before checking the new exports.

The verifier re-evaluates the archived complete rate grid, conditional linear solves, objectives, observed curvature, errors and flags. It does not repeat optimizer searches or regenerate figures. Integrity receipts cover plots and outputs. Synthetic tests also exercise fresh searches, avoiding dependence on optimizer last digits across platforms.

The approach-fit curvature condition number is a numerical diagnostic of a
3-by-3 Hessian. For ill-conditioned cases, its relative last digits can vary
between LAPACK implementations. The verifier compares its reciprocal
\(\sigma_{\min}/\sigma_{\max}\), using relative tolerance \(3\times10^{-7}\)
and an absolute floor of eight float64 machine epsilons. It separately requires
the same classification relative to the existing \(10^{12}\) condition limit.
Fit parameters, propagated errors, objectives, data, validity and resolution
flags keep their original checks. No output estimates or plots are changed
by this portability correction.

To regenerate in a separate empty directory:

    OPENBLAS_NUM_THREADS=1 python3 "$qx_analysis/run_analysis.py" \
      --output temp/qx_radius_plateaus_check

An existing nonempty destination requires explicit --replace. Run numerical replay on Stromboli; use the lightweight integrity check for SuperMUC synchronization.

## Remaining work

Review complete time plots with signed signal and paired window/cutoff changes. Select distances and source choices where more measurements or an improved operator basis answer an unresolved question. Filling unmeasured late-time integer radii needs new MPI work.

The mixed-source matrix pilot remains the existing three-configuration validation, not a 71-configuration GEVP result. A full matrix campaign and smooth/Gaussian Laplacian-eigenvalue profiles require separate planning and reuse decisions. Optimized clover profiles follow a validated source/energy determination.
