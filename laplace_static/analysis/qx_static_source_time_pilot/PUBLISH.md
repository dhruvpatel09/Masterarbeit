# Publish the completed pilot, then pause

These commands publish the existing VS Code spell-checker edit in a dedicated
commit, followed by a separate analysis commit. No cluster jobs are submitted.
All staging archives and bundles stay inside the respective repository roots.

## 1. SuperMUC: commit the existing settings edit

The reported base is `69368a6`. The explicit base/staging checks stop before
committing if another change has already advanced that state. Review the
displayed settings diff; the commit stages only this already tracked file.

```bash
(
    set -e
    cd /hppfs/work/pn29se/di24gus/Masterarbeit
    test "$(git branch --show-current)" = main
    test "$(git rev-parse HEAD)" = 69368a64e4f9a4982641c255ae11a2b713e7145f
    git diff --cached --quiet
    git --no-pager diff -- .vscode/settings.json
    git add -- .vscode/settings.json
    git commit -m "Update VS Code spell-checker dictionary"
    git bundle create qx_source_time_settings_sync.bundle main \
        ^69368a64e4f9a4982641c255ae11a2b713e7145f
    git bundle verify qx_source_time_settings_sync.bundle
    git status
)
```

## 2. Stromboli: import settings, verify analysis, publish both commits

Download `qx_static_source_time_pilot_results.zip` into the Stromboli repository
root before running this block. The original native review is already inside
the analysis package; no second review transfer is needed.

```bash
(
    set -e
    cd /home/m2130292/Masterarbeit
    test "$(git branch --show-current)" = main
    test -z "$(git status --porcelain --untracked-files=no)"
    git pull --ff-only origin main

    scp di24gus@skx.supermuc.lrz.de:/hppfs/work/pn29se/di24gus/Masterarbeit/qx_source_time_settings_sync.bundle .
    git bundle verify qx_source_time_settings_sync.bundle
    git fetch ./qx_source_time_settings_sync.bundle main
    git merge --ff-only FETCH_HEAD

    unzip -n qx_static_source_time_pilot_results.zip
    qx_analysis=laplace_static/analysis/qx_static_source_time_pilot
    python3 "$qx_analysis/verify_files.py"
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python3 "$qx_analysis/test_analysis.py"
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python3 "$qx_analysis/reproduce.py" --verify

    git add -- "$qx_analysis"
    git add -f -- "$qx_analysis/inputs/job6036903_review.zip" \
        "$qx_analysis/inputs/published/laplace_static/analysis/qx_static_matrix_71configs/inputs/primary_matrices.npz"
    git ls-files --error-unmatch -- "$qx_analysis/inputs/job6036903_review.zip" \
        "$qx_analysis/inputs/published/laplace_static/analysis/qx_static_matrix_71configs/inputs/primary_matrices.npz"
    git commit -m "Add qX source-time pilot results and position controls"
    git push origin main

    git bundle create qx_source_time_results_sync.bundle main \
        ^69368a64e4f9a4982641c255ae11a2b713e7145f
    git bundle verify qx_source_time_results_sync.bundle
    scp qx_source_time_results_sync.bundle \
        di24gus@skx.supermuc.lrz.de:/hppfs/work/pn29se/di24gus/Masterarbeit/

    rm -- qx_source_time_settings_sync.bundle \
        qx_source_time_results_sync.bundle qx_static_source_time_pilot_results.zip
    git log -2 --oneline
    git status
)
```

The returned bundle contains both commits. SuperMUC already has the settings
commit, so the final import fast-forwards it to the analysis commit.

## 3. SuperMUC: import and verify, then stop

```bash
(
    set -e
    cd /hppfs/work/pn29se/di24gus/Masterarbeit
    git bundle verify qx_source_time_results_sync.bundle
    git fetch ./qx_source_time_results_sync.bundle main
    git merge --ff-only FETCH_HEAD

    qx_analysis=laplace_static/analysis/qx_static_source_time_pilot
    python3 "$qx_analysis/verify_files.py"
    git ls-files --error-unmatch -- "$qx_analysis/inputs/job6036903_review.zip" \
        "$qx_analysis/inputs/published/laplace_static/analysis/qx_static_matrix_71configs/inputs/primary_matrices.npz"
    qx_previous_origin=$(git rev-parse refs/remotes/origin/main)
    git update-ref refs/remotes/origin/main "$(git rev-parse HEAD)" "$qx_previous_origin"

    rm -- qx_source_time_results_sync.bundle qx_source_time_settings_sync.bundle
    git log -2 --oneline
    git status
)
```

Keep the native campaign directories, builds and saved modes. The root review
archive can remain until the published analysis and its tracked exact input
copy have been verified on both clusters. This workflow ends here; the next
measurement campaign is deferred.
