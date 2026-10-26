/** @file qcd_smearing.h
 *  @brief header file for qcd_smearing.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_SMEARING
#define H_QCD_SMEARING 1

/* prototypes */
/**
 * @fn int qcd_gaussIteration3d(qcd_vector *v, qcd_gaugeField *u, double alpha, int t)
 * @brief  Performs 1 iteration of Gaussian smearing on a time-slice of a spinor field.
 * @details The smearing is performed on time-slice t of the spinor field v. The
 *          gauge field u is used to construct the 3d covariant Hopping matrix H[U].\n
 *          \f$ v = \frac{1}{1+6\alpha}[v+\alpha\,H[U]\, v)  \f$,\n
 *          where H[U] is the 3d Hopping matrix.
 *          Only the necessary communications are performed, based on the flags.
 *          All processes must call this function simultaneously.
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] alpha parameter for the Gauss smearing (positive, real)
 * @param[in] t the time-slice which to Gauss-smear
 * @return 0 on success, 1 on failure
 * @todo communication: do not necessarily communicate U. Use 3d communication.
 */
int qcd_gaussIteration3d(qcd_vector *v, qcd_gaugeField *u, double alpha, int t);

/**
 * @fn int qcd_gaussIteration3dAll(qcd_vector *v, qcd_gaugeField *u, double alpha)
 * @brief Perform one iteration of 3d Gaussian smearing on all time-slices of a spinor field.
 * @details The smearing is performed on all time slices of the spinor field v. The
 *          gauge field u is used to construct the 3d covariant Hopping H[U] matrix.\n
 *          \f$ v = \frac{1}{1+6\alpha}[v+\alpha\,H[U]\, v)  \f$,\n
 *          Only the necessary communications are performed, based on the flags.
 *          All processes must call this function simultaneously.
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] alpha parameter for the Gauss smearing (positive, real)
 * @return 0 on success, 1 on failure
 */
int qcd_gaussIteration3dAll(qcd_vector *v, qcd_gaugeField *u, double alpha);

/**
 * @fn int qcd_apeSmear3d(qcd_gaugeField *apeu, qcd_gaugeField *u, double alpha)
 * @brief Performs one iteration of 3d APE smearing on a gauge field.
 * @details Gauge field u is 3d APE smeared with parameter alpha. The result 
 *          is stored in apeu.
 *          In 3d APE smearing, the temporal links stay untouched. The spatial
 *          links are replaced by the sum of the old link with alpha times 
 *          the sum of spatial staples. The result is projected back to SU(3).
 *          All processes must call this function simultaneously.
 * @param[in,out] apeu pointer to an initialized qcd_gaugeField structure
 * @param[in] u  pointer to an initialized qcd_gaugeField structure
 * @param[in] alpha parameter for the APE smearing (positive, real)
 * @return 0 on success, 1 on failure
 */
int qcd_apeSmear3d(qcd_gaugeField *apeu, qcd_gaugeField *u, double alpha);

#endif 
