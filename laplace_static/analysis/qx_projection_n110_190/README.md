# Reproduce the nine-configuration analysis

Read REPORT.md for the result and its scope. SUMMARY.json locates each extremum.
The two completed raw campaigns are included; full HDF5 data are not.

From this directory, run:

```bash
python3 test_nine_config_analysis.py
python3 analyze_nine_configurations.py inputs/campaign_5918419 inputs/campaign_5920223 recomputed
```

The analysis uses the Python standard library. Decimal precision=50 protects
the small difference after vacuum subtraction. The estimator averages raw
configuration measurements before ratio/subtraction; both prescriptions share
each jackknife deletion. There is no bias-corrected central estimate or inferred
autocorrelation time. The supplied uncertainty is explicitly preliminary.

Optional standalone figures require matplotlib:

```bash
python3 plot_nine_configurations.py recomputed
```

The pre-generated PNGs show longitudinal profiles, paired projection shifts,
and per-configuration central-path SU(3) maxima. All scanned points (including
transverse offsets and axis diagnostics) remain in the CSVs. The first two
figures only plot tr1=tr2=0; table maxima cover all combined-axis scan points.

campaign_audit_*.json records actual stdout/export/provenance revalidation.
audit_profile_campaign.py can repeat the full audit if the original review ZIPs
are extracted; the inputs here contain only the collected scientific tables.

From Masterarbeit, verify the delivery's analysis package with:

```bash
sha256sum -c laplace_static/analysis/qx_projection_n110_190/package.sha256
```
