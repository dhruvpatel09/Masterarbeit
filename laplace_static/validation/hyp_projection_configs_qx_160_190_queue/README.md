# Sequential controller for n160–n190 paired profiles

This controller watches an already submitted index0 job, then submits indices
1–11 one at a time. The scientific campaign is in the sibling directory
hyp_projection_configs_qx_160_190; its README contains the complete build,
preparation, start and resume instructions.

The initial job may still be pending or running when --first-job is supplied.
The controller records its ID and waits for successful actual outputs,
including the n100 regression, before admitting any later case. RUNS.json
persists submission state. A lock prevents concurrent controllers; accepted
IDs are never blindly resubmitted. Failures stop the campaign with records
retained. Collection requires all twelve cases and rechecks their outputs.

Collection creates a review ZIP directly in Masterarbeit. This controller
does not alter gauge fields, saved eigenvectors or solver libraries.
