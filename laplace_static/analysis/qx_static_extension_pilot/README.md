# qX temporal-extension and mode-cutoff pilot

Job 5957304 measured configurations n170, n640 and n650. The scan uses fixed
source time 90, temporal separations 0 through 32, and spatial distances
1,2,4,8,12,16,20,23,24,25,28,32 on the 192 x 64^3 qX lattice. Every spatial
source and all three spatial axes contribute. One temporal HYP step uses
alpha=(1,1,0.5); piter counts projection iterations. The p10 scan includes
Nv=40,60,80,100, and the matched p40 control includes Nv=100.

All 77 review-archive files are retained byte for byte. Their original paths,
stored paths and SHA256 hashes are in `inputs/ARCHIVE_FILES.json`. Archived
build output is stored under `build_provenance`; the original Git attribute
and ignore files are stored as text data so that they do not change tracking
of the archived files. Linked executables and libraries are represented by the
original build hashes; they are not bundled in this analysis package.

## Validated measurements

- 9 MPI launch logs, 9,900 mode audit records on 33 time slices per configuration,
  23,976 correlator records, and 28,800 correlated regional records replayed.
- All 5,112 overlapping archived reference records are reproduced exactly.
- Maximum relative mode residual: 8.2343651e-12, below the unchanged 1e-9 limit.
- Maximum p10/p40 correlator component difference: 1.3234890e-23.
- Central transporter defects remain at roundoff. Direct contractions, source
  counts, axis reductions and complete periodic payload roundtrips pass.

The three configurations are selected technical controls. They establish that
the extended implementation reproduces the previous scan and can be used for
the full-71 campaign. They do not establish population-level errors, a final
cutoff choice, a ground-state plateau, or a string tension.

## Diagnostic exports and plots

`correlator_diagnostics_pilot.csv` forms ensemble means over all three
configuration samples. The normalized correlator is the ratio of ensemble
means, not the mean of per-configuration ratios. `effective_energies_pilot.csv`
contains log ratios of adjacent ensemble means. Errors propagate the full
paired covariance of the configuration samples. Delete-one jackknife errors
are also exported when all three replicates have positive mean endpoints.
`cutoff_differences_pilot.csv` preserves the covariance between cutoffs.

Negative measurements are retained. A nonpositive endpoint mean makes its
effective energy undefined; that value is an empty CSV field. No absolute
values, positive-configuration selection, or silent deletion of failed
jackknife replicates is used. A mean/SE value below three is an illustrative
plotting flag, not a physics acceptance criterion. With only three selected
configurations, these error and signal ratios can fluctuate substantially.

The three figure pairs show:

- `normalized_correlators_mode_cutoffs_pilot`: signed normalized correlators,
  including negative means, with illustrative paired naive errors.
- `effective_energies_mode_cutoffs_pilot`: log ratios, with gaps and flags for
  weak or nonpositive endpoint means. Display ranges emphasize the central
  trends; the CSV retains every numerical value.
- `signed_signal_support_pilot`: signed mean/SE maps at all selected distances.
  Crosses mark nonpositive ensemble means; color values saturate at +/-5.

Lower Nv gives a larger R=24 signal in this pilot, while the improvement does
not extend uniformly to other distances such as R=32. Changing Nv changes the
source operator and its state overlaps. The relevant comparison is the full
ensemble's paired effective energies and their temporal-window stability, not
which cutoff gives the most attractive pilot curve. All four cutoffs therefore
remain in production. The statistical unit is one configuration, not a spatial
source, axis, time separation or cutoff.

## Production timing and continuation

Audit plus p10 scan took 580.21, 581.17 and 580.87 seconds for n170, n640 and
n650. The existing planner therefore admits two configurations per task.
After reusing the three pilot p10 cases, the remaining 68 configurations require
34 tasks. The observed estimate is 19.4 minutes per two-configuration task.
Nine waves with four tasks continuously allocated take about 2.9 hours before
queue gaps, launch overhead and runtime variation. This is not a start-time or
completion-time forecast. Every allocation retains the one-hour limit.

On SuperMUC, using the existing validated workflow and campaign directories:

```bash
python3 laplace_static/validation/qx_static_extension/campaign.py prepare \
  --name late_modes_full71 --stage full --pilot late_modes_pilot
python3 laplace_static/validation/qx_static_extension/campaign.py submit \
  --name late_modes_full71 --start 0 --count 34 --concurrent 4
```

After the array finishes, collect with `campaign.py collect --name
late_modes_full71`. Collection must validate all 71 configurations, including
the three reused pilot cases. Keep the underlying saved-mode, recovery and
campaign directories intact. The production output will support matched
cutoff/time comparisons before new plateau or Cornell fits.

## Reproduction

Core verification needs Python 3.6 or newer and the standard library:

```bash
python3 laplace_static/analysis/qx_static_extension_pilot/test_analysis.py
python3 laplace_static/analysis/qx_static_extension_pilot/reproduce.py --verify
```

The tests additionally replay the captured C contraction and periodic-shift
kernels when GCC is available. The actual MPI/HDF5 execution is the validated
SuperMUC run. Re-rendering figures needs Matplotlib and NumPy:

```bash
python3 laplace_static/analysis/qx_static_extension_pilot/render_plots.py
```

The package integrity inventory checks all preserved inputs, code and published
outputs. Raw stdout is re-parsed, admission markers and launch hashes are
checked, and every CSV diagnostic is recomputed during numerical verification.
