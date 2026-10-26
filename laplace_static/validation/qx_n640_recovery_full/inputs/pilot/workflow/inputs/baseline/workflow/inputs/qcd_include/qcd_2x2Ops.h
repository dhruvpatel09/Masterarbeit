/** @file qcd_2x2Ops.h
 *  @brief header file for qcd_2x2Ops.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2022
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_2x2
#define H_QCD_2x2 1
/** @def qcd_MUL2x2(a,b,c)
 *  @brief a = b*c, where a,b,c are (complex) 2x2 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\tt b[i][k]}\, {\tt c[k][j]}\f$
 */
#define qcd_MUL2x2(a,b,c)  { a[0][0].re =  b[0][0].re*c[0][0].re - b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[1][0].re - b[0][1].im*c[1][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          +b[0][1].re*c[1][0].im + b[0][1].im*c[1][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re - b[0][0].im*c[0][1].im \
                                          +b[0][1].re*c[1][1].re - b[0][1].im*c[1][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im + b[0][0].im*c[0][1].re \
                                          +b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re - b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re - b[1][1].im*c[1][0].im;\
                             a[1][0].im =  b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im + b[1][1].im*c[1][0].re;\
                             a[1][1].re =  b[1][0].re*c[0][1].re - b[1][0].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re - b[1][1].im*c[1][1].im;\
                             a[1][1].im =  b[1][0].re*c[0][1].im + b[1][0].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re;\
}
 

/** @def qcd_MULADJOINT2x2(a,b,c)
 *  @brief a = b* adjoint(c), where a,b,c are (complex) 2x2 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\tt b[i][k]}\, {\rm conj}({\tt c[j][k]}) \f$
 */
#define qcd_MULADJOINT2x2(a,b,c) { \
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im;\
                             a[0][0].im = -b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          -b[0][1].re*c[0][1].im + b[0][1].im*c[0][1].re;\
                             a[0][1].re =  b[0][0].re*c[1][0].re + b[0][0].im*c[1][0].im \
                                          +b[0][1].re*c[1][1].re + b[0][1].im*c[1][1].im;\
                             a[0][1].im = -b[0][0].re*c[1][0].im + b[0][0].im*c[1][0].re \
                                          -b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re + b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[0][1].re + b[1][1].im*c[0][1].im;\
                             a[1][0].im = -b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          -b[1][1].re*c[0][1].im + b[1][1].im*c[0][1].re;\
                             a[1][1].re =  b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im;\
                             a[1][1].im = -b[1][0].re*c[1][0].im + b[1][0].im*c[1][0].re \
                                          -b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re;\
}


/** @def qcd_ADJOINTMUL2x2(a,b,c)
 *  @brief a = adjoint(b)*c, where a,b,c are (complex) 2x2 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\rm conj}({\tt b[k][i]})\,{\tt c[k][j]} \f$
 */
#define qcd_ADJOINTMUL2x2(a,b,c)  {\
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im - b[0][0].im*c[0][0].re \
                                          +b[1][0].re*c[1][0].im - b[1][0].im*c[1][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re + b[0][0].im*c[0][1].im \
                                          +b[1][0].re*c[1][1].re + b[1][0].im*c[1][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im - b[0][0].im*c[0][1].re \
                                          +b[1][0].re*c[1][1].im - b[1][0].im*c[1][1].re;\
                                                                                         \
                             a[1][0].re =  b[0][1].re*c[0][0].re + b[0][1].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re + b[1][1].im*c[1][0].im;\
                             a[1][0].im =  b[0][1].re*c[0][0].im - b[0][1].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im - b[1][1].im*c[1][0].re;\
                             a[1][1].re =  b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im,\
                             a[1][1].im =  b[0][1].re*c[0][1].im - b[0][1].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im - b[1][1].im*c[1][1].re;\
}


/** @def qcd_ADJOINTMULADJOINT2x2(a,b,c)
 *  @brief a = adjoint(b)*adjoint(c), where a,b,c are (complex) 2x2 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\rm conj}({\tt b[k][i]})\, 
 *                                             {\rm conj}({\tt c[j][k]}) \f$
 */
#define qcd_ADJOINTMULADJOINT2x2(a,b,c)  { \
                             a[0][0].re =  b[0][0].re*c[0][0].re - b[0][0].im*c[0][0].im \
                                          +b[1][0].re*c[0][1].re - b[1][0].im*c[0][1].im;\
                             a[0][0].im = -b[0][0].re*c[0][0].im - b[0][0].im*c[0][0].re \
                                          -b[1][0].re*c[0][1].im - b[1][0].im*c[0][1].re;\
                             a[0][1].re =  b[0][0].re*c[1][0].re - b[0][0].im*c[1][0].im \
                                          +b[1][0].re*c[1][1].re - b[1][0].im*c[1][1].im;\
                             a[0][1].im = -b[0][0].re*c[1][0].im - b[0][0].im*c[1][0].re \
                                          -b[1][0].re*c[1][1].im - b[1][0].im*c[1][1].re;\
                                                                                         \
                             a[1][0].re =  b[0][1].re*c[0][0].re - b[0][1].im*c[0][0].im \
                                          +b[1][1].re*c[0][1].re - b[1][1].im*c[0][1].im;\
                             a[1][0].im = -b[0][1].re*c[0][0].im - b[0][1].im*c[0][0].re \
                                          -b[1][1].re*c[0][1].im - b[1][1].im*c[0][1].re;\
                             a[1][1].re =  b[0][1].re*c[1][0].re - b[0][1].im*c[1][0].im \
                                          +b[1][1].re*c[1][1].re - b[1][1].im*c[1][1].im;\
                             a[1][1].im = -b[0][1].re*c[1][0].im - b[0][1].im*c[1][0].re \
                                          -b[1][1].re*c[1][1].im - b[1][1].im*c[1][1].re;\
}

/**
 * @def qcd_SU2TRACER(u)
 * @brief Real part of the trace of u.
 * @details \f$ {\tt result} = \sum\limits_k {\rm Re}({\tt u[k][k]}) \f$
 */
#define qcd_SU2TRACER(u) ( (double) (u[0][0].re + u[1][1].re))

/******************************************************************************
 * prototypes 
 ******************************************************************************/

/**
 * @fn void qcd_zero2x2(qcd_complex_16 a[2][2])
 * @brief Sets the complex 2x2 matrix a to zero.
 * @param[in,out] a pointer to a complex 2x2 matrix
 */
void qcd_zero2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_unit2x2(qcd_complex_16 a[2][2])
 * @brief Sets the complex 2x2 matrix to a unit matrix.
 * @param[in,out] a pointer to a complex 2x2 matrix
 */
void qcd_unit2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_minusUnit2x2(qcd_complex_16 a[2][2])
 * @brief Sets the complex 2x2 matrix to minus 1.
 * @details Where 1 is a complex unit matrix.
 * @param[in,out] a pointer to a complex 2x2 matrix
 */
void qcd_minusUnit2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_sub2x2(qcd_complex_16 difference[2][2], qcd_complex_16 minuend[2][2], qcd_complex_16 subtrahend[2][2]);
 * @brief Subtracts two complex 2x2 matrices from each other.
 * @details difference = minuend - subtrahend
 * @param[out] difference pointer to a complex 2x2 matrix
 * @param[in] minuend pointer to a complex 2x2 matrix
 * @param[in] subtrahend pointer to a complex 2x2 matrix
 */
void qcd_sub2x2(qcd_complex_16 difference[2][2], qcd_complex_16 minuend[2][2], qcd_complex_16 subtrahend[2][2]);

/**
 * @fn void qcd_add2x2(qcd_complex_16 sum[2][2], qcd_complex_16 summand1[2][2], qcd_complex_16 summand2[2][2]);
 * @brief Sums two complex 2x2 matrices.
 * @details sum = summand1 + summand2\n
 * Two or three of the parameters may be identical.
 * @param[out] sum pointer to a complex 2x2 matrix
 * @param[in] summand1 pointer to a complex 2x2 matrix
 * @param[in] summand2 pointer to a complex 2x2 matrix
 */
void qcd_add2x2(qcd_complex_16 sum[2][2], qcd_complex_16 summand1[2][2], qcd_complex_16 summand2[2][2]);

/**
 * @fn void qcd_addAdjoint2x2(qcd_complex_16 sum[2][2], qcd_complex_16 summand1[2][2], qcd_complex_16 summand2[2][2]);
 * @brief Adds a 2x2 matrix to the adjoint of another.
 * @details sum = summand1 + adj(summand2)\n
 * Two or three of the parameters may be identical.
 * @param[out] sum pointer to a complex 2x2 matrix
 * @param[in] summand1 pointer to a complex 2x2 matrix
 * @param[in] summand2 pointer to a complex 2x2 matrix
 */
void qcd_addAdjoint2x2(qcd_complex_16 sum[2][2], qcd_complex_16 summand1[2][2], qcd_complex_16 summand2[2][2]);

/**
 * @fn void qcd_mul2x2(qcd_complex_16 product[2][2], qcd_complex_16 factor1[2][2], qcd_complex_16 factor2[2][2])
 * @brief Multiplies two complex 2x2 matrices
 * @details product = factor1 * factor2\n
 * factor1 and factor2 may be identical.
 * @param[out] product pointer to a complex 2x2 matrix
 * @param[in] factor1 pointer to a complex 2x2 matrix
 * @param[in] factor2 pointer to a complex 2x2 matrix
 */
void qcd_mul2x2(qcd_complex_16 product[2][2], qcd_complex_16 factor1[2][2], qcd_complex_16 factor2[2][2]);

/**
 * @fn void qcd_mulAdjoint2x2(qcd_complex_16 product[2][2], qcd_complex_16 factor1[2][2], qcd_complex_16 factor2[2][2])
 * @brief Multiply a 2x2 matrix with the adjoint of another
 * @details product = factor1 * adj(factor2)
 * @param[out] product pointer to a complex 2x2 matrix
 * @param[in] factor1 pointer to a complex 2x2 matrix
 * @param[in] factor2 pointer to a complex 2x2 matrix
 */
void qcd_mulAdjoint2x2(qcd_complex_16 product[2][2], qcd_complex_16 factor1[2][2], qcd_complex_16 factor2[2][2]);

/**
 * @fn void qcd_scale2x2(qcd_complex_16 a[2][2], double r)
 * @brief Scales a complex 2x2 matrix by a real number
 * @details a = a*r
 * @param[in,out] a pointer to a complex 2x2 matrix
 * @param[in] r real rescaling factor
 */
void qcd_scale2x2(qcd_complex_16 a[2][2], double r);

/**
 * @fn void qcd_cScale2x2(qcd_complex_16 a[2][2], qcd_complex_16 c)
 * @brief Scales a complex 2x2 matrix by a complex number
 * @details a = a*c
 * @param[in,out] a pointer to a complex 2x2 matrix
 * @param[in] c complex rescaling factor
 */
void qcd_cScale2x2(qcd_complex_16 a[2][2], qcd_complex_16 c);

/**
 * @fn qcd_complex_16 qcd_trace2x2(qcd_complex_16 a[2][2])
 * @brief Computes the trace of a complex 2x2 matrix.
 * @details trace = a[0][0]+a[1][1]+a[2][2]
 * @param[in] a pointer to a complex 2x2 matrix
 * @return the (complex) trace of the 2x2 matrix
 */
qcd_complex_16 qcd_trace2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_copy2x2(qcd_complex_16 dest[2][2], qcd_complex_16 src[2][2])
 * @brief Copies a complex 2x2 matrix
 * @details dest[i][j] = src[i][j] for i,j=0,1,2
 * @param[out] dest a pointer to a complex 2x2 matrix
 * @param[in] src pointer to a complex 2x2 matrix
 */
void qcd_copy2x2(qcd_complex_16 dest[2][2], qcd_complex_16 src[2][2]);

/**
 * @fn double qcd_norm2x2(qcd_complex_16 a[2][2]);
 * @brief Computes the Froebenius norm of a complex 2x2 matrix
 * @details nrm = sqrt(a^dag*a)
 * @param[in] a pointer to a complex 2x2 matrix
 * @return the real, positive norm
 */
double qcd_norm2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_adjoint2x2(qcd_complex_16 a[2][2])
 * @brief Computes the adjoint of a complex 2x2 matrix
 * @details a[i][j] = conj(a[j][i]) for i,j=0,1
 * @param[in,out] a pointer to a complex 2x2 matrix
 */
void qcd_adjoint2x2(qcd_complex_16 a[2][2]);

/**
 * @fn qcd_complex_16 qcd_det2x2(qcd_complex_16 a[2][2])
 * @brief Computes the determinant of a complex 2x2 matrix
 * @details a[i][j] = a[0][0]*a[1][1]-a[0][1]*a[1][0]
 * @param[in] a pointer to a complex 2x2 matrix
 * @return the determinant of the matrix
 */
qcd_complex_16 qcd_det2x2(qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_eig2x2(qcd_complex_16 l[2], qcd_complex_16 a[2][2])
 * @brief Computes the eigenvalues of a complex 2x2 matrix
 * @details l[0] and l[1] are the two complex eigenvalues
 * @param[in] a pointer to a complex 2x2 matrix
 * @param[out] l pointer to a complex 2x1 array to store the eigenvalues
 */
void qcd_eig2x2(qcd_complex_16 l[2], qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_inv2x2(qcd_complex_16 ai[2][2], qcd_complex_16 a[2][2])
 * @brief Computes the inverse of a complex 2x2 matrix
 * @details ai is the inverse of the complex 2x2 matrix a. If a is singular,
 *          the inverse will contain infinite entries caused by divisions by zero.
 * @param[in] a pointer to a complex 2x2 matrix
 * @param[out] ai pointer to a complex 2x2 matrix
 */
void qcd_inv2x2(qcd_complex_16 ai[2][2], qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_2x2ToQuat(qcd_complex_16 y[4], qcd_complex_16 a[2][2])
 * @brief Decomposition of a 2x2 Matrix into a linear combination of Pauli matrices and 1
 * @details Every complex 2x2 matrix can be written as\n
 *          \f$ A = {\tt y[0]} {\bf 1} + {\tt y[1]} \sigma_1 + {\tt y[2]} \sigma_2 + {\tt y[3]} \sigma_3 \f$ \n
 *          with complex coefficients y. This functions computes these coefficients.
 * @param[in] a pointer to a complex 2x2 matrix
 * @param[in,out] y pointer to a 4-component array with complex coefficients
 */
void qcd_2x2ToQuat(qcd_complex_16 y[4], qcd_complex_16 a[2][2]);

/**
 * @fn void qcd_SU2ToQuat(double y[4], qcd_complex_16 a[2][2])
 * @brief Decomposition of a 2x2 Matrix into a linear combination of Pauli matrices and 1
 * @details \f$a=\sqrt{2}\f$\n
 *         Every SU(2) matrix can be written as\n
 *           \n
 *          with real coefficients y. This functions computes these coefficients.
 *          In the case of a non-SU(2) matrix the coefficients will not be normalized and
 *          only their real part is computed.
 * @param[in] a pointer to a complex 2x2 matrix
 * @param[in,out] y pointer to a 4-component array with real coefficients
 */
void qcd_SU2ToQuat(double y[4], qcd_complex_16 a[2][2]);

void qcd_quatTo2x2(qcd_complex_16 a[2][2], qcd_complex_16 y[4]);
void qcd_quatToSU2(qcd_complex_16 a[2][2], double y[4]);
void qcd_3x3To2x2(qcd_complex_16 W2x2[2][2], qcd_complex_16 W[3][3], int k);
void qcd_2x2To3x3(qcd_complex_16 W[3][3], qcd_complex_16 W2x2[2][2], int k);
void qcd_mul2x23x3(qcd_complex_16 P[3][3], qcd_complex_16 F2x2[2][2], qcd_complex_16 F[3][3], int k);
#endif
