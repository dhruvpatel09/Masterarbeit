# Complete qX static-correlator study: 71 configurations

All supplied configurations n100–800 in steps of 10 are included, including recovered n640. The production baseline is piter=10, Nv=100, rho=1, one temporal HYP step with alpha=(1,1,0.5), a fixed source t0/a=90, every tau/a=0–16 and every R/a=1–32. The lattice is 192 x 64^3, with open temporal and periodic spatial boundaries. The three validated pilot p10 measurements are reused; 68 configurations were measured by array 5955993.

The correlator is the existing gauge-invariant Laplacian static correlator. It is not a thin-spatial-string Wilson loop. Spatial source translations and the three axes are averaged within each configuration. There are 71 statistical configuration samples.

## Data and validation

The archive contains 157,052 production correlator rows and 272,640 regional records. Reproduction reparses all 71 cases and 145 MPI launches, including the pilot p40 controls, and checks the assembled tables against the logs. The selected saved-mode audits cover 1,207 time slices and 120,700 modes. The largest relative residual is 1.530e-11, below the unchanged 1e-9 limit; the largest production p10 central unitarity defect is 7.956e-16.

Every original review-archive file is preserved under `inputs/`, with a path mapping and SHA256 hashes. The original archive hash is recorded in `SOURCES.json`. The reused pilot attempts are identified by the immutable full-stage preparation record: the full review archive includes their attempt directories but not the pilot case-level LATEST markers. No missing marker is fabricated. Empty stderr files are explicitly retained by Git rules.

## Estimates and errors

The ensemble assessment establishes thermalization and approximate independence of the supplied saved configurations. The baseline error is therefore the independent-configuration sample error. The effective energy is

```
aE_eff(R,tau) = log(mean(L(R,tau)) / mean(L(R,tau+1)))
```

It is not the average of per-configuration logarithms. The linear propagation retains the full covariance between the two primary means. Paired delete-one-configuration errors provide a nonlinear cross-check. Nonpositive means and unsupported leave-one-out logarithms are flagged; no configurations are discarded, and no absolute values or clipping are applied to the correlators.

All 544 combined-axis production means are positive, although 197 individual combined-axis measurements are nonpositive. The weakest combined-axis mean/SE is about 1.83, near R/a=24 at late time. Positivity of a mean alone does not establish a reliable fit region.

The optional PyErrors `gamma_method(S=0)` check covers all 512 combined-axis effective-energy errors with the same ensemble label and all 71 configuration IDs. Its measured version record and comparisons are included under `results/`. No additional autocorrelation inflation is imposed.

## Temporal and spatial fit scans

The complete settings are in `ANALYSIS_PLAN.json`. The reference constant-energy window is tau/a=6–11. A fixed grid of 17 temporal windows is exported for each R and each axis channel. Correlated generalized least squares uses the full time covariance. The covariance weights are held fixed for the paired delete-one error check. Fits report covariance conditioning and approximate Gaussian chi-square tail probabilities. These probabilities are diagnostics for nonlinear means and an estimated covariance; they are not exact finite-sample tests or multiple-comparison-corrected plateau decisions.

`diagnostic_compatible` requires approximate p >= 0.05, minimum primary mean/SE >= 5 throughout the fitted interval, and covariance infinity-norm condition <= 1e12. These explicit screening settings do not define a physics acceptance criterion. Window differences retain their paired configuration covariance.

The full sample resolves later-time dependence. For example, at R/a=12 the reference estimate is approximately 0.198207(678), while the 8–13 window gives approximately 0.196452. The paired shift is about 3.19 standard errors. A constant fit can have an acceptable p value within one interval while still shifting in a later interval.

Cornell fits use `aV(R)=V0 + sigma_a2*R - alpha/R`, with covariance across R. Eleven spatial ranges are scanned in four temporal windows, giving 44 candidate fits. The fitted coefficient and statistical error are exported for transparency, together with goodness of fit and temporal-window flags. Strong short-distance model failures and temporal/range dependence prevent quoting a finalized string tension from this scan. An isolated compatible candidate is not a final selection after scanning many ranges.

The next measurement priority is a controlled longer-time check at moderate separations whose signals remain strong but whose window estimates drift. The weak region near R/a=24 also merits an operator-overlap study. Physical comparison with Em1 requires the qX scale and a common additive-potential convention.

## Figures

- [Correlators and signal map](results/correlator_signals_71configs.png)
- [Effective energies](results/effective_energies_71configs.png)
- [Temporal-window stability](results/temporal_window_stability_71configs.png)
- [Potential and Cornell diagnostics](results/potential_and_cornell_diagnostics_71configs.png)
- [Fit-quality scans](results/fit_quality_scans_71configs.png)
- [Regional covariance diagnostics](results/regional_covariance_diagnostics_71configs.png)

The regional comparison pools 64 source cubes per axis, or 192 region/axis channels in the combined average. It compares the configuration error with a calculation that artificially omits covariance among these regional channels. Ratios can be above or below one: the complete configuration covariance remains the baseline. Regional channels do not increase the number of independent configurations beyond 71.

## Reproduction

The core requires Python 3.6 or newer and NumPy. Figures additionally require Matplotlib. SciPy is an optional reference for one test. PyErrors 2.17.0 is an optional automatic-differentiation cross-check in a modern Python environment. The saved PyErrors comparison can be verified without installing PyErrors.

```bash
python3 laplace_static/analysis/qx_static_potential_71configs/test_analysis.py
python3 laplace_static/analysis/qx_static_potential_71configs/reproduce.py --verify
python3 laplace_static/analysis/qx_static_potential_71configs/render_plots.py
```

Primary means, nonlinear estimates and temporal covariance exports use Decimal arithmetic at precision 50 and are reproduced byte-for-byte. Fit and regional diagnostics use floating-point linear algebra and are compared with relative tolerance 3e-9 and absolute tolerance 1e-14. Input receipts pin the executed workflow and all 45 source/provenance files. Large gauge/eigenvector files and compiled executables are represented by the original receipts rather than copied into Git.

For error propagation and covariance-fit context, see the PyErrors documentation and Bruno–Sommer, arXiv:2209.14188. The present fit-quality scan uses the explicitly labelled Gaussian approximation; it does not implement the paper's more general p-value construction.

## Publishing from Stromboli

Place `qx_static_potential_71configs_results.zip` in the Stromboli repository root:

```bash
(
    set -e
    cd /home/m2130292/Masterarbeit
    git pull --ff-only origin main
    unzip -n qx_static_potential_71configs_results.zip
    python3 laplace_static/analysis/qx_static_potential_71configs/test_analysis.py
    python3 laplace_static/analysis/qx_static_potential_71configs/reproduce.py --verify
    git add laplace_static/analysis/qx_static_potential_71configs
    git commit -m "Add complete qX static-correlator and fit diagnostics"
    git push origin main
    git bundle create qx_static_full71_results_sync.bundle main \
      ^d1d4546328de88645592edde04ad8bf9bfbd60f6
    scp qx_static_full71_results_sync.bundle \
      di24gus@skx.supermuc.lrz.de:/hppfs/work/pn29se/di24gus/Masterarbeit/qx_static_full71_results_sync.bundle
)
```

After the push and transfer succeed, import on SuperMUC. The copy-integrity check uses only the Python standard library:

```bash
(
    set -e
    cd /hppfs/work/pn29se/di24gus/Masterarbeit
    git bundle verify qx_static_full71_results_sync.bundle
    git fetch ./qx_static_full71_results_sync.bundle main
    git merge --ff-only FETCH_HEAD
    git update-ref refs/remotes/origin/main "$(git rev-parse HEAD)" \
      d1d4546328de88645592edde04ad8bf9bfbd60f6
    python3 laplace_static/analysis/qx_static_potential_71configs/verify_files.py
    rm -- qx_static_full71_results_sync.bundle
    git status
)
```

After successful destination verification, remove the Stromboli transfer copies with `rm -- qx_static_full71_results_sync.bundle qx_static_potential_71configs_results.zip`. Retain campaign case directories and saved-mode data for future measurements.
