Stage 1: Laplace static trial-state correlator

Goal:
Implement Höllwieser et al. Eq. 7 and Eq. 8 using already computed Nv=10 eigenvectors.

Inputs:
- Gauge configs: /home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n*
- Eigenvectors: /home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n*/eigenvectors
- Eigenvalues: /home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n*/eigenvalues

First implementation:
- Nv = 10
- rho_i = 1
- on-axis separations first: R/a = 1,2,3,4
- T/a = 1,...,12
- average over spatial source x and source time t0
- output L(R,T)

Target formulas:
tau_ij(y,t0,t1) = v_i^\dagger(y,t0) U_t(y,t0,t1) v_j(y,t1)

L(R,T) = < sum_ij tau_ij(y,t0,t1) tau_ji(x,t1,t0) >
