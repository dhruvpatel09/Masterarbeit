/** @file qcd_gaugeFixing.h
 *  @brief header file for qcd_gaugeFixing.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_GAUGEFIXING
#define H_QCD_GAUGEFIXING 1

/**
 * @def qcd_LANDAU_GAUGE_FIXING_PRECISION
 * @brief Stopping criterion for Landau gauge fixing.
 */
#define qcd_LANDAU_GAUGE_FIXING_PRECISION 1.0e-12

/**
 * @def qcd_LANDAU_GAUGE_FIXING_MAX_ITER
 * @brief Maximal iteration count for Landau gauge fixing.
 */
#define qcd_LANDAU_GAUGE_FIXING_MAX_ITER 10000

/**
 * @fn void qcd_gluonField(qcd_gaugeField *a, qcd_gaugeField *u)
 * @brief Computes the lattice gluon field from a gauge field.
 * @details The lattice gluon field a is calculated from the 
 * SU(3) gauge-field u:\n
 * a = 1/(2i) [u-u<sup>+</sup>] - trace
 * @param[in,out] a pointer to an initialized qcd_gaugeField structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @todo This function computes a also on the outer boundary. Maybe: make it
 *       optional.
 */
void qcd_gluonField(qcd_gaugeField *a, qcd_gaugeField *u);

/**
 * @fn void qcd_localGluonField(qcd_complex_16 a[3][3], qcd_gaugeField *u, long x, int mu)
 * @brief Computes the lattice gluon field on a single link.
 * @details The lattice gluon field a is computed from a single link of the
 *          SU(3) gauge field u:\n
 *          a = 1/(2i) [u(x,mu)-u<sup>+</sup>(x,mu)] - trace
 * @param[out] a pointer to a complex 3x3 matrix
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] x the gauge field on link (x,mu) is used
 * @param[in] mu the gauge field on link (x,mu) is used
 */
void qcd_localGluonField(qcd_complex_16 a[3][3], qcd_gaugeField *u, long x, int mu);


/**
 * @fn int qcd_gaugeTransform(qcd_gaugeField *u, qcd_gaugeTransformation *g)
 * @brief Applies local gauge transformation g to gaugeField u
 * @details The qcd_gaugeField u gets replaced by\n
 *          U<sub>&mu;</sub>(x) -> g<sup>+</sup>(x) U<sub>&mu;</sub>(x) g(x+&mu;)
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @param[in] g pointer to an initialized qcd_gaugeTransformation structure
 * @return 0 on success, 1 on failure
 */
int qcd_gaugeTransform(qcd_gaugeField *u, qcd_gaugeTransformation *g);

/**
 * @fn int qcd_gaugeTransformVector(qcd_vector *v, qcd_gaugeTransformation *g)
 * @brief Applies local gauge transformation g to spinor v
 * @details The spinor v gets replaced by\n
 *          v(x) -> g(x) v(x)\n
 *          If the outer boundaries of g and v are up to date, the transformation
 *          is also carried out on those.
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in] g pointer to an initialized qcd_gaugeTransformation structure
 * @return 0 on success, 1 on failure
 */
int qcd_gaugeTransformVector(qcd_vector *v, qcd_gaugeTransformation *g);

/**
 * @fn int qcd_landauGauge(qcd_gaugeField *landauu, qcd_gaugeField *u, double overparam)
 * @brief Landau gauge fixing procedure.
 * @details The stochastic over-relaxation algorithm to fix u to landau-gauge.
 *          The result is stored in landauu.\n
 *          A local maximum of the functional
 *          &sum;<sub>&mu;,x</sub> Re Tr [ U<sub>&mu;</sub> (x) ]
 *          is returned. The iteration is stopped when the precision 
 *          qcd_LANDAU_GAUGE_FIXING_PRECISION is reached
 *          for the quality-measure 
 *          1/V &sum;<sub>x</sub> Tr(A<sub>&mu;</sub>(x)-A<sub>&mu;</sub>(x-&mu;))<sup>2</sup>
 *          or the iteration number surpasses qcd_LANDAU_GAUGE_FIXING_MAX_ITER.
 *          This is partially based on a fortran routine by Constantia Alexandrou.
 * @param[in,out] landauu pointer to an initialized qcd_gaugeField structure
 * @param u pointer to an initialized qcd_gaugeField structure
 * @param overparam overrelaxation parameter 1&le;overpararm&ge;2
 * @return 0 on success, 1 on failure
 * @warning Low quality random numbers are used 
 */
int qcd_landauGauge(qcd_gaugeField *landauu, qcd_gaugeField *u, double overparam);


int qcd_maxTreeGauge(qcd_gaugeField *u);

#endif 
