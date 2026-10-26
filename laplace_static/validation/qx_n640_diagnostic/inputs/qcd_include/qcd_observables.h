/** @file qcd_observables.h
 *  @brief header file for qcd_observables.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_OBSERVABLES
#define H_QCD_OBSERVABLES 1  

/* prototypes */
/**
 * @fn double qcd_calculatePlaquette(qcd_gaugeField *u)
 * @brief Computes the average plaquette
 * @details The real part of the trace of a plaquette is computed and averaged
 *          over the whole global lattice. All processes must call this 
 *          function simultaneously.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return The average plaquette value
 */
double qcd_calculatePlaquette(qcd_gaugeField *u);

/**
 * @fn int qcd_calculateCloverField(qcd_cloverField *q, qcd_gaugeField *u)
 * @brief Calculates the clover field.
 * @details Given a gauge-field u, the clover field q is computed:\n
 *          q->D[x][munu][c1][c2] contains the c1/c2 color component of the 
 *          clover field \f$ Q_{\mu\nu}(x)\f$. munu=0..5 corresponds to planes
 *          mu/nu=0/1,0/2,0/3,1/2,1/3,2/3. The definition of the clover field is
 *\f$ Q_{\mu\nu}(x) = \ldots \f$
 *          All processes must call this function simultaneously.
 * @param[in,out] q pointer to an initialized qcd_cloverField structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_calculateCloverField(qcd_cloverField *q, qcd_gaugeField *u);

/**
 * @fn double qcd_calculateWeightedPlaquette(qcd_gaugeField *u, double *w0, double *wk, double normfac)
 * @brief Calculates the average plaquette, with different weights for spatial and
 *        temporal plaquettes on different time-slices.
 * @details Given the gauge field u, two arrays w0 and wk of length geo->lT,
 *          and an overall normalization factor normfac, the average plaquette 
 *          is computed. Temporal plaquettes with positive corner (local corrdinates)
 *          (lt,lx,ly,lz) are weighted with w0[lt] and spatial ones with wk[lt].
 *          The final result is normalized (i.e. divided by) normfac and returned.
 *          All processes must call this function simultaneously.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] w0 array of weights for temporal plaquettes
 * @param[in] wk array of weights for spatial plaquettes
 * @param[in] normfac overall normalization
 * @return the weighted average plaquette
 */
double qcd_calculateWeightedPlaquette(qcd_gaugeField *u, double *w0, double *wk, double normfac);

/**
 * @fn int qcd_isSFGaugeField(qcd_gaugeField *u)
 * @brief Tests whether a gauge field is a Schroedinger Functional one.
 * @details In the Schroedinger functional, gauge field components at global
 *          time coordinate x0=T-1, in the mu=0 direction, are zero. 
 *          At global times x0=0 the spatial gauge fields are constant. 
 *          At global times x0=T-1 the spatial gauge fields are constant (but
 *          possibly different than at x0=0).
 *          A 1 is returned if u has these properties, otherwise a 0 is returned.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return 0 or 1
 */
int qcd_isSFGaugeField(qcd_gaugeField *u);


/**
 * @fn int qcd_getSFBoundaries(double *phi, double *phiP, qcd_gaugeField *u);
 * @brief Extracts the Abelian phases of Schroedinger Functional boundary conditions.
 * @details It is assumed that u has SF boundaries (test with #qcd_isSFGaugeField).
 *          Then from the values of u at global
 *          time x0=0 phi is extracted assuming that u has the form \n
 *          \f$U_k(x) =\begin{pmatrix}e^{i{\tt phi[0]}} & 0                 & 0 \\
                                         0              & e^{i{\tt phi[1]}} & 0 \\
                                         0              & 0                 & e^{i{\tt phi[2]}} 
                       \end{pmatrix}\f$\n
 *          and at global time T-1 phip is extracted assuming\n
 *          \f$U_k(x) =\begin{pmatrix}e^{i{\tt phip[0]}} & 0                  & 0 \\
                                         0               & e^{i{\tt phip[1]}} & 0 \\
                                         0               & 0                  & e^{i{\tt phip[2]}} 
                       \end{pmatrix}\f$ \n
 *          A 0 is returned on success, a 1 otherwise.
 *          All processes must to call this function simultaneously.
 * @param[out] phi Abelian phase for the boundaries at x0=0
 * @param[out] phiP Abelian phase for the boundaries at x0=T-1
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return 0 or 1
 */
int qcd_getSFBoundaries(double *phi, double *phiP, qcd_gaugeField *u);

/**
 * @fn int qcd_hasSFAngles(qcd_gaugeField *u, double *phi, double *phip)
 * @brief Tests whether a gaguge field has Abelian Schroedinger Functional boundaries
 *        with the given angles.
 * @details It is assumed that u has SF boundaries (test with #qcd_isSFGaugeField).
 *          For these, it is tested whether spatial gauge-fields at global
 *          time x0=0 have the form \n
 *          \f$U_k(x) =\begin{pmatrix}e^{i{\tt phi[0]}} & 0                 & 0 \\
                                         0              & e^{i{\tt phi[1]}} & 0 \\
                                         0              & 0                 & e^{i{\tt phi[2]}} 
                       \end{pmatrix}\f$\n
 *          and at global time T-1 have the form\n
 *          \f$U_k(x) =\begin{pmatrix}e^{i{\tt phip[0]}} & 0                  & 0 \\
                                         0               & e^{i{\tt phip[1]}} & 0 \\
                                         0               & 0                  & e^{i{\tt phip[2]}} 
                       \end{pmatrix}\f$ \n
 *          If this is the case, a
 *          1 is returned, otherwise a 0.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] phi Abelian phase for the boundaries at x0=0
 * @param[in] phip Abelian phase for the boundaries at x0=T-1
 * @return 0 or 1
 * @todo change parameter order to match getSFAngles
 */
int qcd_hasSFAngles(qcd_gaugeField *u, double *phi, double *phip);

/**
 * @fn int qcd_isOpenGaugeField(qcd_gaugeField *u)
 * @brief Tests whether a gauge field has open boundary conditions
 * @details A gauge field u with open boundary conditions is zero at
 *          global time x0=T-1 in direction mu=0. A 1 is returned if this
 *          is the case, a 1 otherwise. All processes must call this function
 *          simultaneously.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return 0 or 1
 * @bug This assumes that u=0 <=> norm(u)=0, which is not the case with the present norm
 */
int qcd_isOpenGaugeField(qcd_gaugeField *u);

#endif
