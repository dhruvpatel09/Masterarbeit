# qX static correlators: all 71 configurations, four mode cutoffs

This analysis assembles the temporal-extension campaign on all qX configurations
n100 through n800 in steps of ten, including recovered n640. Production array
5958012 contributes 68 configurations; the three validated p10 pilot cases from
array 5957304 are reused. Matched p40 controls on n170, n640 and n650 are replayed
as technical controls and do not add configuration samples.

## Geometry and observable

The lattice is 192 x 64^3 with open temporal and periodic spatial boundaries.
The static correlator uses the saved spatial Laplacian modes and temporal HYP
transporters. The source time is fixed at t0/a=90; temporal separations are
tau/a=0,...,32. The measured spatial distances are
R/a=1,2,4,8,12,16,20,23,24,25,28,32, along all three axes and with every spatial
source included. These are twelve selected distances, not a dense R=1,...,32
scan at every new temporal separation.

The same 100 saved modes are read at each endpoint. The p10 launch evaluates
projectors with Nv=40,60,80,100 using the lowest stored modes; their eigenvalue
ordering is audited. One temporal HYP step uses alpha=(1,1,0.5). The value
`piter=10` counts projection iterations and is distinct from the mode count.
Changing Nv changes the source operator and its state overlaps. It is not a
sequence of estimates with independent errors, nor an extrapolation to an exact
answer at Nv=100. These static measurements contain no clover insertions.

## Completeness and technical validation

Every admitted case is re-parsed from stdout. Launch hashes, exit codes,
prepared-case receipts and collected tables are checked. The complete archive
contains 831 files; all are preserved byte for byte, with original-to-stored
paths and SHA256 hashes in `inputs/ARCHIVE_FILES.json`. Archived build files are
stored in `build_provenance`; original Git attribute/ignore files are stored as
text data so they do not affect tracking of archived evidence. Executables and
linked libraries are represented by their original build hashes.

Validation covers:

- all 71 configurations, 145 MPI launches, 452,412 production correlator rows
  and 340,800 correlated regional rows;
- 234,300 saved modes on 33 selected time slices per configuration;
- all 63,048 overlapping reference records, reproduced exactly;
- direct versus factorized contractions for every cutoff, spatial coordinate
  shifts, source counts, axis reductions and complete periodic payload shifts;
- temporal transporter links at t/a=90,...,121, plus preservation of the other
  field components.

The maximum audited relative mode residual is 1.5842035e-11, below the unchanged
1e-9 limit. The maximum central unitarity defect is 7.9557074e-16. The largest
p10/p40 correlator component difference on the matched controls is 1.3234890e-23.
These are implementation checks; physical plateau stability is assessed
separately.

## Statistical treatment

The supplied configuration sequence is treated as thermalized and approximately
independent, following the ensemble information. One statistical sample is one
configuration after spatial-source and axis averaging. There are 71 samples;
spatial sources, axes, times and cutoffs do not multiply that count.

Primary mean and log arithmetic uses 50 decimal digits. The naive error of a
mean is the sample standard deviation divided by sqrt(71). Normalized
correlators and effective energies use ratios of ensemble means. Full paired
covariance is propagated through these nonlinear functions. Paired delete-one
jackknife errors are exported as a cross-check; a nonlinear jackknife result is
reported only if every omitted-configuration replicate is defined.

All signed configuration measurements are retained. A nonpositive ensemble
mean at either endpoint makes an effective energy undefined; its CSV fields
are empty. No absolute values or positive-configuration selection are used.
There are 1,536 possible effective-energy points, of which 1,494 are defined.
The 42 undefined points remain in the exports with explicit flags.

PyErrors 2.17.0 at S=0 reproduces the errors of all 1,494 defined energies and
1,102 defined paired cutoff contrasts. Maximum relative discrepancies are
1.45e-15 and 4.22e-15 respectively. The optional check requires modern Python;
the core numerical analysis supports Python 3.6 with NumPy. S=0 checks the
independent-sample covariance algebra, not autocorrelation independence or
Gaussian behaviour of weak logarithms.

For points whose two endpoint means have mean/SE at least three, the observed
jackknife/naive error ratios lie between approximately 0.990 and 1.012. The
mean/SE=3 threshold is an illustrative plotting flag, not an agreed physics
acceptance criterion.

## Correlated fits and interpretation

The report scans twelve six-point effective-energy windows, starting at
tau/a=4,6,...,26. A window 6...11 uses correlator means at times 6...12.
Constant fits retain the full covariance between their six effective energies.
The full-sample GLS weights are held fixed in the paired jackknife. Differences
between windows and cutoffs retain their configuration covariance.

Cornell diagnostics fit aE(R)=V0-alpha/R+(sigma a^2)R to five explicitly listed
sparse spatial sets. They retain covariance between spatial distances and use
the same fixed-weight jackknife convention. The fitted parameter is a lattice
quantity; a physical scale and controlled fit systematics would be required
before a physical string-tension comparison with another ensemble.

The chi-square tail probabilities are approximate Gaussian reference values
with a sample-estimated covariance. They are exploratory and unadjusted across
the many windows, cutoffs and ranges. An acceptable value in one window does
not establish a ground-state plateau. Small-R lattice effects, periodic images,
source overlaps and open-boundary/source-placement effects remain relevant to
physical interpretation.

The full sample confirms a stronger R=24 signal for lower Nv. At tau/a=24,
mean/SE is about 14.8 for Nv=40 and 0.50 for Nv=100. At R=32 the corresponding
values are about -1.07 and 3.52: the improvement is not uniform in distance.
The negative value is a weak signed mean, not a negative physical potential.

Later windows bring several cutoff estimates closer. Residual paired
differences remain, for example at R=12 for Nv=40 versus Nv=100 in window
14...19. This is consistent with remaining operator-overlap or excited-state
effects and does not identify a unique preferred cutoff. The fitted energies
and Cornell slopes in this package are diagnostics, not a final potential or
string tension.

## Plots and tables

Seven PNG/SVG figure pairs are provided in `results`:

- `normalized_correlators_71configs`: common-source normalized signals;
- `effective_energies_71configs`: logarithmic energy estimates and signal flags;
- `signed_signal_support_71configs`: signed mean/SE over all selected geometries;
- `fit_window_stability_71configs`: correlated constant-fit window scans;
- `fitted_energies_vs_distance_71configs`: comparison of three displayed windows;
- `paired_cutoff_differences_71configs`: covariance-resolved cutoff differences;
- `cornell_window_diagnostics_71configs`: time/range dependence of Cornell slopes.

Hollow markers flag weak endpoint means, gaps identify undefined inputs, and
triangles mark values beyond a chosen display range. Complete values and flags
are retained in the CSV files. Close or overlapping curves contain all four
cutoffs; the paired-difference plot resolves their smaller differences.

The energy covariance table stores all cutoff/time pairs at each distance.
Its 128-component covariance has rank at most 70; fits use smaller blocks.
The fitted-energy covariance table additionally preserves correlations between
distances and cutoffs at each fit window. These tables and the explicit paired
contrast tables allow comparisons without adding correlated errors as if
they were independent.

## Reproduction and publication

On Stromboli or another Python/NumPy environment:

```bash
python3 laplace_static/analysis/qx_static_extension_71configs/test_analysis.py
python3 laplace_static/analysis/qx_static_extension_71configs/reproduce.py --verify
```

This replays every case and recomputes means, covariance, paired errors,
constant fits, window/cutoff contrasts and Cornell diagnostics. The estimator
module is reused unchanged from the published full71 static-correlator
analysis; its upstream receipt is recorded in `inputs/ESTIMATOR_SOURCE.json`.
The statistical tests include closed-form correlated-fit weights, cancellation
of common fluctuations, signed inputs and singular-covariance rejection.
Captured C contraction and complete periodic-shift kernels are tested when GCC
is available. The actual MPI/HDF5 execution is the validated SuperMUC campaign.

For package integrity on a SuperMUC login without NumPy:

```bash
python3 laplace_static/analysis/qx_static_extension_71configs/verify_files.py
```

Add `--replay` to repeat standard-library input replay. Plot regeneration needs
Matplotlib; the optional PyErrors cross-check needs a modern Python environment:

```bash
python3 laplace_static/analysis/qx_static_extension_71configs/render_plots.py
python3 laplace_static/analysis/qx_static_extension_71configs/pyerrors_check.py
```

Publish from Stromboli and synchronize the commit to SuperMUC using a bundle.
Keep transfer files inside the repository root and remove them after successful
publication and destination verification. Keep the underlying saved-mode,
n640 recovery and campaign directories intact.
