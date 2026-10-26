# Complete qX ensemble: 71 configurations

The matched cohort comprises n100–800 in steps of 10, including the recovered n640 configuration. Main raw profiles use piter=4,10,20,40. The matched raw/traceless comparison uses piter=4,40. Every prescription uses the same configurations and the three central separations tau/a=2,4,6, midpoint t/a=98, Nv=100 and R/a=1,2,4.

The inputs combine the preserved 70-configuration tables with the independently replayed n640 measurement archive from job 5946008. `inputs/SOURCES.json` pins the earlier tables, estimators, parsers and all archived n640 files. Original 70-configuration numerical results remain unchanged and their profile and temporal marginals are reproduced exactly during verification. `results/n640_inclusion_comparison.csv` isolates the correlated change from adding n640.

The estimator is the ratio of configuration means, Re(mean(LO))/Re(mean(L)) - mean(O), with equal configuration weights. Raw clover is the six-plane sum Tr(F†F), with F=(Q-Q†)/8 on the original gauge field. The traceless side observable subtracts sum |Tr(F)|²/3; neither convention has an additional factor 1/2. Ratios, temporal contrasts and trace fractions are recomputed within synchronized jackknife omissions. Errors remain preliminary without an autocorrelation correction.

Main estimates use all 71 configurations. Conditional binning uses complete equal-size bins anchored at n100. Widths 1,2,4,8 retain 71,70,68,64 configurations, respectively; all omitted labels are listed in `results/binning_plans_71configs.json`. Error ratios compare blocked and unblocked estimates on exactly the same retained samples. With n640 restored, bins after the former gap are regrouped. Increasing labels remain a provisional ordering; no MDU conversion or replica assignment is assumed.

Traceless piter10/20 primaries measured on n640 alone do not complete those prescriptions on the full cohort. The n100/n170/n670/n740 piter6/8/80 study remains a separate four-configuration diagnostic. No temporal plateau, ground-state saturation, physical normalization or statistical independence is inferred here.

From the repository root:

```bash
python3 laplace_static/analysis/qx_complete_ensemble_71configs/test_analysis.py
python3 laplace_static/analysis/qx_complete_ensemble_71configs/reproduce.py --verify
```

Verification uses Python's standard library and requires no HPC gauge field or eigenvector files. Numerical regeneration uses `reproduce.py`; add `--plots` to regenerate PNG/SVG figures with Matplotlib and NumPy. Output figures separate the piter4 differences from the piter10/20 rounding scale and include a zoom of the near-roundoff link diagnostics.
