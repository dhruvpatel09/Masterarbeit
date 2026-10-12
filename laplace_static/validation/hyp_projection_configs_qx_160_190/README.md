# Paired temporal-HYP profiles: qX n160–n190

Saved-mode audit 5920205 passed for 5200 modes on 52 selected slices, with the
full eigenvalue table checked at all 192 times. The largest relative residual
was 7.705107848820442e-11; none exceeded the requested 1e-9. Maximum norm defect
was 1.887379141862766e-15; maximum Gram diagonal defect 3.774758283725532e-15;
maximum off-diagonal overlap 4.625071495221476e-11. All four MPI exits were zero
and application stderr files were empty. Actual output records were re-parsed
and matched against the exported CSVs. This is sampled saved-data validation
using the existing qcd operator, not full eigenvector coverage or an independent
operator implementation.

## Measurements and scope

Twelve cases cover configurations 160,170,180,190 at tau=2,4,6, midpoint98.
Each case starts independently from the original gauge field for piter4 and
piter40: one temporal HYP step, alpha=(1,1,0.5), Nv100, rho=1, R=1,2,4 and the
same 332 probe/axis records as the completed n110–n150 campaign. All six endpoint
slices (95,96,97,99,100,101) and midpoint98 were included in audit5920205.

The clover probe retains the current RAW convention on the original field:
no traceless subtraction and no extra factor1/2. Spatial APE20, alpha0.5 defines
the eigenvector operator and is separate from temporal HYP. The numerical C
source, SU(3) metrics header, measurement launcher and scientific checkers are
byte-identical to those in the completed 15-case campaign. Changes select new
configurations, pin the current quality evidence and adapt case bookkeeping.

Each measurement uses the previously validated observable grid1x4x4x4, 64 MPI
ranks on eight nodes, one thread per rank and solver_mpiInit/FUNNELED. This
differs from the eigenvector-generation/audit grid6x4x4x4. The eight-node
observable allocation and launcher are retained from the successful campaign.

Case0 additionally repeats n100 tau4 at piter4 and piter40 and checks all 664
raw reference records from job5915129. Its success is required before later
cases are admitted. Outputs enforce profile coverage, counts, translations,
halo checks, direct samples, MPI reductions, transporter factorization,
unchanged spatial/boundary links and unchanged original-field clover probes.
SU(3) diagnostics cover the physical temporal links in each central transporter
window, not boundary tails elsewhere. piter40 must meet the existing 1e-12
control tolerance. piter4 defects and observable sensitivity are reported
without an automatic physics-acceptance threshold.

Connected CSV columns are per-configuration diagnostics. Form ensemble raw
averages of L, LS and S before computing mean(LS_RE)/mean(L_RE)-mean(S_RAW).
Configurations supply the ensemble samples; sites, axes and tau values do not
provide additional independent configurations. Assess correlated uncertainty
and absolute differences along with relative changes near small signals.

## Start on SuperMUC

Place hyp_projection_configs_qx_160_190_integration.zip in Masterarbeit:

```bash
(
    set -euo pipefail
    cd /hppfs/work/pn29se/di24gus/Masterarbeit

    profiles=laplace_static/validation/hyp_projection_configs_qx_160_190
    controller=laplace_static/validation/hyp_projection_configs_qx_160_190_queue
    test ! -e "$profiles"
    test ! -e "$controller"
    unzip hyp_projection_configs_qx_160_190_integration.zip
    sha256sum -c "$profiles/package.sha256"
    sha256sum -c "$controller/package.sha256"

    bash "$profiles/build.sh"
    python3 "$profiles/prepare.py"

    first_job=$(sbatch --parsable --array=0 "$profiles/submit.sh")
    first_job="${first_job%%;*}"
    echo "Initial profile job: $first_job"
    nohup python3 -u "$controller/run_remaining.py" \
        --first-job "$first_job" \
        > "$controller/controller.out" 2>&1 < /dev/null &
)
```

The controller watches the submitted initial job whether pending, running or
already finished. It validates the actual outputs before submitting index1,
then proceeds one case at a time through index11. It retains each accepted job
ID, checks provenance before future submissions, and does not infer success
from disappearance from squeue. Failures stop further submissions. Only one
controller may run. Existing data and old results are preserved.

The prior paired runs took roughly seven minutes per case; case0 takes about
fourteen minutes because it includes the additional regression pair. Budget
roughly 1.5–2 hours of execution for the twelve cases, plus queue waits. Each
job uses partition/qos test and a 30-minute limit. No new eigenvectors are
generated; this campaign produces small logs and CSVs.

Monitor from Masterarbeit:

```bash
tail -n 25 laplace_static/validation/hyp_projection_configs_qx_160_190_queue/controller.out
squeue -u "$USER" -o '%.18i %.14j %.2t %.10M %.6D %R'
```

After all twelve cases pass, the controller rechecks the outputs and creates
hyp_projection_configs_160_190_review_mixed_FIRST_JOB.zip directly inside
Masterarbeit. Expected aggregate coverage is 7968 raw profile rows, 3984 paired
rows, 96 SU(3) time rows and 24 SU(3) window rows, excluding the separate n100
regression. Upload that ZIP for assessment and combination with the earlier
five-configuration pilot.

If the controller stops, inspect its error and RUNS.json before retrying. After
resolving the issue, resume the recorded campaign on one login node:

```bash
controller=laplace_static/validation/hyp_projection_configs_qx_160_190_queue
nohup python3 -u "$controller/run_remaining.py" --resume \
    >> "$controller/controller.out" 2>&1 < /dev/null &
```

Ambiguous sbatch acceptance is retained for inspection and is not blindly
resubmitted. Failed measurement outputs are preserved. Do not start a second
campaign to bypass a failure or edit checksum files to bypass changed inputs.

## Local validation

All fifteen actual completed profile cases were rechecked using the inherited
scientific checkers. Eleven controller tests cover sequential admission,
pending/running initial jobs, failure stops, submission timeout handling,
recorded acceptance recovery, provenance and identity mismatches, complete
coverage, partial-result labels and resume without duplicate submission.
Scientific-check stubs in controller fixtures test orchestration only, not new
physics. Python and Bash syntax, byte identity of numerical sources, package
checksums and ZIP integrity are verified. Compilation/linking and the new
measurements still run on SuperMUC through the commands above.

adaptation.diff records changes from the completed campaign and its controller.
SAVED_MODE_REVIEW.json retains the successful current audit summary and scope.
