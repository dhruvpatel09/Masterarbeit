# Historical four-iteration HYP2 projection on qX

Reference: Donnellan, Knechtli, Leder and Sommer, arXiv:1012.3037, section 2.1.
They use alpha=(1,1,0.5), the normalization in hep-lat/0506008 Eq.(2.24),
then four updates of Eq.(2.25). The finite analytic map is part of the historical
definition. Four projection updates do not mean four complete HYP steps.

This diagnostic uses one temporal HYP step and independently resets the
original qX_p2n100 field before piter4, repeat4 and piter40. It retains the
original field and uses solver_mpiInit with the verified 6x4x4x4 solver build.
No eigenvectors or HDF5 outputs are written.

From the Masterarbeit root:

    bash laplace_static/validation/hyp_projection_qx/build.sh
    sbatch laplace_static/validation/hyp_projection_qx/submit.sh

Structural checks cover import/export, zero alpha, unchanged spatial/absent
upper-boundary temporal links, temporal effect, finite components, repeatability
and the number of physical temporal links. The expected count is
191*64^3 = 50069504. Gauge allocation/init and MPI calls follow the prior test.

For each physical temporal link define E=U^dagger U-I. The output reports:

- UNIT_MAX_ENTRY: maximum complex magnitude of an entry of E, matching the
  earlier temporal-unitarity diagnostic.
- UNIT_MAX_FROB and UNIT_RMS_FROB: maximum and RMS of the Frobenius norm of E
  across links (RMS is sqrt(mean(norm(E)^2))).
- DET_MAX and DET_RMS: maximum and RMS of abs(det(U)-1).
- SU3_EXCEED: counts/fractions above 1e-12, 1e-8, 1e-4 and 1e-2; these are
  diagnostic thresholds, not advisor-approved acceptance bounds.
- PITER_4_VS_40: maximum real/imaginary link-component differences.

The converged piter40 reference retains a 1e-12 SU(3) tolerance. No SU(3)
acceptance threshold or agreement-to-piter40 threshold is imposed on piter4.
PASS means structural checks and repeatability passed; it does not approve the
piter4 prescription for production, quantify observable effects, or independently
validate all HYP staples. Review defects and then compare actual L, LS and S
measurements using the original clover convention. Existing piter40 results
remain useful as diagnostics, not substitutes for the historical definition.

Local verification: C syntax with the reviewed qcd/solver headers and local MPI
declarations; maximum/RMS defects and exceedance fractions checked against
independent NumPy calculations for 63 known/random matrices, including SU(3)
matrices and perturbed matrices. An invalid absent upper-boundary link was
excluded as intended. This checks the diagnostic mathematics, not MPI reduction
behavior or HYP execution on the cluster; those checks remain pending.
