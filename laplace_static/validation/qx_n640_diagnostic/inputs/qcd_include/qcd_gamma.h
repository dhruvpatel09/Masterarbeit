/** @file qcd_gamma.h
 *  @brief header file for qcd_gamma.c
 *  @details Details & license. ALPHA conventions
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2012-2015
 *  @copyright GNU Public License.
 */


#ifndef H_QCD_GAMMA
#define H_QCD_GAMMA 1

/** 
 * @var qcd_EPS[6][3]
 * @brief Indices of non-vanishing components of the 3d epsilon tensor.
 * @details &epsilon;<sub>ijk</sub> has only six non-vanishing entries:\n
 *          &epsilon;<sub>012</sub> = &epsilon;<sub>201</sub> = &epsilon;<sub>120</sub> = +1 and\n
 *          &epsilon;<sub>210</sub> = &epsilon;<sub>021</sub> = &epsilon;<sub>102</sub> = -1.\n
 *          #qcd_EPS stores the indices that belong to these six entries.
 *          The actual value +1 or -1 at these indices is stored in #qcd_SGN_EPS. \n
 *          e.g. qcd_EPS[0]= {0,1,2}, qcd_SGN_EPS[0] = +1
 */
extern const int qcd_EPS[6][3];

/** 
 * @var qcd_SGN_EPS[6][3]
 * @brief Values of the 3d epsilon tensor.
 * @details &epsilon;<sub>ijk</sub> has only six non-vanishing entries:\n
 *          &epsilon;<sub>012</sub> = &epsilon;<sub>201</sub> = &epsilon;<sub>120</sub> = +1 and\n
 *          &epsilon;<sub>210</sub> = &epsilon;<sub>021</sub> = &epsilon;<sub>102</sub> = -1.\n
 *          #qcd_EPS stores the indices that belong to these six entries.
 *          The actual value +1 or -1 at these indices is stored in #qcd_SGN_EPS. \n
 *          e.g. qcd_EPS[0]= {0,1,2}, qcd_SGN_EPS[0] = +1
 */
extern const int qcd_SGN_EPS[6];

/**
 * @var qcd_complex_16 qcd_ONE[4][4]
 * @brief 4x4 unit matrix
 */
extern const qcd_complex_16 qcd_ONE[4][4];

/**
 * @var qcd_complex_16 qcd_GAMMA[8][4][4]
 * @brief Dirac matrices with ALPHA conventions.
 * @details qcd_Gamma[0]: &gamma;<sub>0</sub> = 
 \verbatim
  [ 0  0 -1  0]
  [ 0  0  0 -1]
  [-1  0  0  0]
  [ 0 -1  0  0]
 \endverbatim
 ** qcd_Gamma[1]: &gamma;<sub>1</sub> =
 \verbatim
  [ 0  0  0 -i]
  [ 0  0 -i  0]
  [ 0  i  0  0]
  [ i  0  0  0]
 \endverbatim
 * qcd_Gamma[2]: &gamma;<sub>2</sub> =
 \verbatim
  [ 0  0  0 -1]
  [ 0  0  1  0]
  [ 0  1  0  0]
  [-1  0  0  0]
 \endverbatim
 * qcd_Gamma[3]: &gamma;<sub>3</sub> =
 \verbatim
  [ 0  0 -i  0]
  [ 0  0  0  i]
  [ i  0  0  0]
  [ 0 -i  0  0]
 \endverbatim
 * qcd_Gamma[4] = -qcd_Gamma[0] \n
 * qcd_Gamma[5]: &gamma;<sub>5</sub> = 
 *                &gamma;<sub>0</sub> &gamma;<sub>1</sub> &gamma;<sub>2</sub> &gamma;<sub>3</sub> 
 *              = &gamma;<sub>1</sub> &gamma;<sub>2</sub> &gamma;<sub>3</sub> &gamma;<sub>4</sub> =
 \verbatim
  [ 1  0  0  0]
  [ 0  1  0  0]
  [ 0  0 -1  0]
  [ 0  0  0 -1]
 \endverbatim
 * qcd_Gamma[6]: &gamma;<sub>+</sub> = (&gamma;<sub>1</sub> + i&gamma;<sub>2</sub>)/2 \n
 * qcd_Gamma[7]: &gamma;<sub>-</sub> = (&gamma;<sub>1</sub> - i&gamma;<sub>2</sub>)/2 \n
 */
extern const qcd_complex_16 qcd_GAMMA[8][4][4];

/**
 * @var qcd_complex_16 qcd_CGAMMA[6][4][4]
 * @brief Dirac matrices multiplied by C.
 * @details The Dirac matrices are in the ALPHA conventions
 *          (see #qcd_GAMMA). The charge conjugation matrix
 *          is given by\n
 *          C = -&gamma;<sub>0</sub>  &gamma;<sub>2</sub>\n
 *          #qcd_CGAMMA contains the matrices \n
 *          C &gamma;<sub>i</sub> with i=0,1,2,3,4,5.
 */
extern const qcd_complex_16 qcd_CGAMMA[6][4][4];

/**
 * @var qcd_complex_16 qcd_BAR_CGAMMA[6][4][4]
 * @brief Matrices <span style="text-decoration:overline">C &gamma;<sub>i</sub></span> = 
 *        &gamma;<sub>0</sub> [C &gamma;<sub>i</sub>]<sup>+</sup> &gamma;<sub>0</sub>
 * @details The Dirac matrices are in the ALPHA conventions
 *          (see #qcd_GAMMA). The charge conjugation matrix
 *          is given by\n
 *          C = -&gamma;<sub>0</sub>  &gamma;<sub>2</sub>\n
 *          #qcd_BAR_CGAMMA contains the matrices \n
 *          &gamma;<sub>0</sub> [C &gamma;<sub>i</sub>]<sup>+</sup> &gamma;<sub>0</sub> with i=0,1,2,3,4,5.
 */
extern const qcd_complex_16 qcd_BAR_CGAMMA[6][4][4];

/**
 * @var qcd_complex_16 qcd_G5GAMMA[8][4][4]
 * @brief Matrices &gamma;<sub>5</sub> &gamma;<sub>i</sub>
 * @details #qcd_G5GAMMA contains the matrices \n
 *          &gamma;<sub>5</sub> &gamma;<sub>i</sub>, i=0,1,2,3,4,5,6,7\n
 *          where &gamma;<sub>i</sub> are defined as in
 *          #qcd_GAMMA.
 */
extern const qcd_complex_16 qcd_G5GAMMA[8][4][4];

/**
 * @var qcd_complex_16 qcd_BAR_G5GAMMA[8][4][4]
 * @brief Matrices <span style="text-decoration:overline">&gamma;<sub>5</sub> &gamma;<sub>i</sub></span> =
 *  &gamma;<sub>0</sub> [&gamma;<sub>5</sub> &gamma;<sub>i</sub>]<sup>+</sup> &gamma;<sub>0</sub>
 * @details The Dirac matrices are in the ALPHA conventions
 *          (see #qcd_GAMMA).
 *          #qcd_BAR_G5GAMMA contains the matrices \n
 *          &gamma;<sub>0</sub> [&gamma;<sub>5</sub> &gamma;<sub>i</sub>]<sup>+</sup> &gamma;<sub>0</sub>
 *           with i=0,1,2,3,4,5,6,7.
 */
extern const qcd_complex_16 qcd_BAR_G5GAMMA[8][4][4];

/**
 * @var qcd_complex_16 qcd_ONE_PLUS_GAMMA[6][4][4]
 * @brief Matrices 1+&gamma;<sub>i</sub>
 * @details The Dirac matrices are in the ALPHA conventions
 *          (see #qcd_GAMMA).
 *          #qcd_ONE_PLUS_GAMMA contains the matrices \n
 *          1+&gamma;<sub>i</sub> with i=0,1,2,3,4,5.
 */
extern const qcd_complex_16 qcd_ONE_PLUS_GAMMA[6][4][4];

/**
 * @var qcd_complex_16 qcd_ONE_MINUS_GAMMA[6][4][4]
 * @brief Matrices 1-&gamma;<sub>i</sub>
 * @details The Dirac matrices are in the ALPHA conventions
 *          (see #qcd_GAMMA).
 *          #qcd_ONE_MINUS_GAMMA contains the matrices \n
 *          1-&gamma;<sub>i</sub> with i=0,1,2,3,4,5.
 */
extern const qcd_complex_16 qcd_ONE_MINUS_GAMMA[6][4][4];

/* prototypes */
void qcd_gamma0Vector(qcd_vector *v);
void qcd_gamma1Vector(qcd_vector *v);
void qcd_gamma2Vector(qcd_vector *v);
void qcd_gamma3Vector(qcd_vector *v);
void qcd_gamma0gamma1Vector(qcd_vector *v);
void qcd_gamma0gamma2Vector(qcd_vector *v);
void qcd_gamma0gamma3Vector(qcd_vector *v);
void qcd_gamma0gamma5Vector(qcd_vector *v);
void qcd_gamma1gamma2Vector(qcd_vector *v);
void qcd_gamma1gamma3Vector(qcd_vector *v);
void qcd_gamma1gamma5Vector(qcd_vector *v);
void qcd_gamma2gamma3Vector(qcd_vector *v);
void qcd_gamma2gamma5Vector(qcd_vector *v);
void qcd_gamma3gamma5Vector(qcd_vector *v);

/**
 * @fn void qcd_gamma5Vector(qcd_vector *v)
 * @brief Applies the &gamma;<sub>5</sub> matrix to a spinor field.
 * @details The &gamma;<sub>5</sub> matrix (see #qcd_GAMMA) is applied
 *          to the spinor field v on the local lattice. The result is stored
 *          again in v.\n
 *          v.D[x][d1][c] = &sum;<sub>d2</sub> qcd_GAMMA5[5][d1][d2] v.D[x][d2][c]
 * @param[in,out] v pointer to an initialized qcd_vector structure
 */
void qcd_gamma5Vector(qcd_vector *v);

/**
 * @fn void qcd_gamma5Propagator(qcd_propagator *p);
 * @brief Applies the &gamma;<sub>5</sub> matrix to a propagator field.
 * @details The &gamma;<sub>5</sub> matrix (see #qcd_GAMMA) is applied
 *          to the propagator  p on the local lattice. The result is stored
 *          again in p.\n
 *          p.D[x][d1][d2][c1][c2] = &sum;<sub>d3</sub> qcd_GAMMA5[5][d1][d3] p.D[x][d3][d2][c1][c2]
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 */
void qcd_gamma5Propagator(qcd_propagator *p);

/**
 * @fn void qcd_propagatorGamma5(qcd_propagator *p);
 * @brief Applies the &gamma;<sub>5</sub> matrix to a propagator field from the right.
 * @details The propagator p is multiplied with the 
 *          &gamma;<sub>5</sub> matrix (see #qcd_GAMMA) on the local lattice. The result is stored
 *          again in p.\n
 *          p.D[x][d1][d2][c1][c2] = &sum;<sub>d3</sub>  p.D[x][d1][d3][c1][c2] qcd_GAMMA5[5][d3][d2]
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 */
void qcd_propagatorGamma5(qcd_propagator *p);

/**
 * @fn void qcd_gamma5PropagatorGamma5(qcd_propagator *p);
 * @brief Applies the &gamma;<sub>5</sub> matrix to a propagator field.
 * @details The propagator p is multiplied from both sides with the &gamma;<sub>5</sub> matrix 
 *          (see #qcd_GAMMA) on the local lattice. The result is stored
 *          again in p.\n
 *          p.D[x][d1][d2][c1][c2] = 
 *          &sum;<sub>d3,d4</sub> qcd_GAMMA5[5][d1][d3] p.D[x][d3][d4][c1][c2] qcd_GAMMA5[5][d4][d2]
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 */
void qcd_gamma5PropagatorGamma5(qcd_propagator *p);

#endif
