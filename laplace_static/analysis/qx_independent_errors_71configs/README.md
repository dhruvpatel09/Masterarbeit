# Independent-configuration errors for qX

The sample contains all 71 supplied configurations, n100 through n800 in steps of 10, including recovered n640. The ensemble assessment establishes thermalization and approximate independence of these saved configurations. One spatially averaged configuration is one statistical sample. There is no new thermalization cut or compulsory binning correction.

For the five matched primary quantities `(L, LS, S, T, LT)`, the connected observables are

```
C_RAW = mean(LS)/mean(L) - mean(S)
C_TRACE = mean(LT)/mean(L) - mean(T)
C_TL = (mean(LS)-mean(LT))/mean(L) - (mean(S)-mean(T))
```

The naive error propagates the complete covariance of these sample means through the nonlinear estimator. It includes correlations among the primary quantities, projection prescriptions and temporal separations. It does not average per-configuration ratios or count translated sources as independent configurations.

`reproduce.py` uses Python's standard library and Decimal arithmetic. All existing full71 point estimates are preserved exactly. The configuration jackknife remains a nonlinear-estimator cross-check. Across the 249 combined-axis positions, the raw and traceless p40 error differences are below 0.012%.

```
python3 laplace_static/analysis/qx_independent_errors_71configs/test_analysis.py
python3 laplace_static/analysis/qx_independent_errors_71configs/reproduce.py --verify
python3 laplace_static/analysis/qx_independent_errors_71configs/render_plots.py
```

`pyerrors_check.py` is an optional, separately recorded automatic-differentiation cross-check in a modern Python environment with `pyerrors==2.17.0`. It uses the same ensemble name and all 71 configuration labels for every primary observable. `gamma_method(S=0)` selects the independent-sample error and reproduces both profile errors and paired temporal-contrast errors. It is the PyErrors counterpart of the `Stau=0` setting in UWerr. The earlier Em1 scripts used nonzero window parameters; their paths and hashes are preserved in `inputs/SOURCES.json`.

```
python3 laplace_static/analysis/qx_independent_errors_71configs/pyerrors_check.py
```

PyErrors is not required for the standard-library verification. Its version record and computed comparison are preserved under `results/PYERRORS_CHECK.json` and `results/pyerrors_S0_checks_71configs.csv`. The latter can be checked against the analytic exports without installing PyErrors.

The present clover measurements retain their original central geometry, R/a=1,2,4 and tau/a=2,4,6. The raw convention remains the primary observable; the traceless convention is a matched side study. Physical normalization, continuum interpretation and a ground-state plateau require separate checks. The new fixed-source static-correlator workflow is under `laplace_static/validation/qx_static_potential`.
