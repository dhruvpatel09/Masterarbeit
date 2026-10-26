/** @file qcd_3x3Ops.h
 *  @brief header file for qcd_3x3Ops.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2022
 *  @copyright GNU Public License.
 *  @todo 3x3: eigenvalues
 */
 
#ifndef H_QCD_3x3
#define H_QCD_3x3 1
/** @def qcd_MUL3x3(a,b,c)
 *  @brief a = b*c, where a,b,c are (complex) 3x3 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\tt b[i][k]}\, {\tt c[k][j]}\f$
 */
#define qcd_MUL3x3(a,b,c)  { a[0][0].re =  b[0][0].re*c[0][0].re - b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[1][0].re - b[0][1].im*c[1][0].im \
                                          +b[0][2].re*c[2][0].re - b[0][2].im*c[2][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          +b[0][1].re*c[1][0].im + b[0][1].im*c[1][0].re \
                                          +b[0][2].re*c[2][0].im + b[0][2].im*c[2][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re - b[0][0].im*c[0][1].im \
                                          +b[0][1].re*c[1][1].re - b[0][1].im*c[1][1].im \
                                          +b[0][2].re*c[2][1].re - b[0][2].im*c[2][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im + b[0][0].im*c[0][1].re \
                                          +b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re \
                                          +b[0][2].re*c[2][1].im + b[0][2].im*c[2][1].re;\
                             a[0][2].re =  b[0][0].re*c[0][2].re - b[0][0].im*c[0][2].im \
                                          +b[0][1].re*c[1][2].re - b[0][1].im*c[1][2].im \
                                          +b[0][2].re*c[2][2].re - b[0][2].im*c[2][2].im;\
                             a[0][2].im =  b[0][0].re*c[0][2].im + b[0][0].im*c[0][2].re \
                                          +b[0][1].re*c[1][2].im + b[0][1].im*c[1][2].re \
                                          +b[0][2].re*c[2][2].im + b[0][2].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re - b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re - b[1][1].im*c[1][0].im \
                                          +b[1][2].re*c[2][0].re - b[1][2].im*c[2][0].im;\
                             a[1][0].im =  b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im + b[1][1].im*c[1][0].re \
                                          +b[1][2].re*c[2][0].im + b[1][2].im*c[2][0].re;\
                             a[1][1].re =  b[1][0].re*c[0][1].re - b[1][0].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re - b[1][1].im*c[1][1].im \
                                          +b[1][2].re*c[2][1].re - b[1][2].im*c[2][1].im;\
                             a[1][1].im =  b[1][0].re*c[0][1].im + b[1][0].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re \
                                          +b[1][2].re*c[2][1].im + b[1][2].im*c[2][1].re;\
                             a[1][2].re =  b[1][0].re*c[0][2].re - b[1][0].im*c[0][2].im \
                                          +b[1][1].re*c[1][2].re - b[1][1].im*c[1][2].im \
                                          +b[1][2].re*c[2][2].re - b[1][2].im*c[2][2].im;\
                             a[1][2].im =  b[1][0].re*c[0][2].im + b[1][0].im*c[0][2].re \
                                          +b[1][1].re*c[1][2].im + b[1][1].im*c[1][2].re \
                                          +b[1][2].re*c[2][2].im + b[1][2].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[2][0].re*c[0][0].re - b[2][0].im*c[0][0].im \
                                          +b[2][1].re*c[1][0].re - b[2][1].im*c[1][0].im \
                                          +b[2][2].re*c[2][0].re - b[2][2].im*c[2][0].im;\
                             a[2][0].im =  b[2][0].re*c[0][0].im + b[2][0].im*c[0][0].re \
                                          +b[2][1].re*c[1][0].im + b[2][1].im*c[1][0].re \
                                          +b[2][2].re*c[2][0].im + b[2][2].im*c[2][0].re;\
                             a[2][1].re =  b[2][0].re*c[0][1].re - b[2][0].im*c[0][1].im \
                                          +b[2][1].re*c[1][1].re - b[2][1].im*c[1][1].im \
                                          +b[2][2].re*c[2][1].re - b[2][2].im*c[2][1].im;\
                             a[2][1].im =  b[2][0].re*c[0][1].im + b[2][0].im*c[0][1].re \
                                          +b[2][1].re*c[1][1].im + b[2][1].im*c[1][1].re \
                                          +b[2][2].re*c[2][1].im + b[2][2].im*c[2][1].re;\
                             a[2][2].re =  b[2][0].re*c[0][2].re - b[2][0].im*c[0][2].im \
                                          +b[2][1].re*c[1][2].re - b[2][1].im*c[1][2].im \
                                          +b[2][2].re*c[2][2].re - b[2][2].im*c[2][2].im;\
                             a[2][2].im =  b[2][0].re*c[0][2].im + b[2][0].im*c[0][2].re \
                                          +b[2][1].re*c[1][2].im + b[2][1].im*c[1][2].re \
                                          +b[2][2].re*c[2][2].im + b[2][2].im*c[2][2].re;\
}
 

/** @def qcd_MULADJOINT3x3(a,b,c)
 *  @brief a = b* adjoint(c), where a,b,c are (complex) 3x3 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\tt b[i][k]}\, {\rm conj}({\tt c[j][k]}) \f$
 */
#define qcd_MULADJOINT3x3(a,b,c) {\
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im \
                                          +b[0][2].re*c[0][2].re + b[0][2].im*c[0][2].im;\
                             a[0][0].im = -b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          -b[0][1].re*c[0][1].im + b[0][1].im*c[0][1].re \
                                          -b[0][2].re*c[0][2].im + b[0][2].im*c[0][2].re;\
                             a[0][1].re =  b[0][0].re*c[1][0].re + b[0][0].im*c[1][0].im \
                                          +b[0][1].re*c[1][1].re + b[0][1].im*c[1][1].im \
                                          +b[0][2].re*c[1][2].re + b[0][2].im*c[1][2].im;\
                             a[0][1].im = -b[0][0].re*c[1][0].im + b[0][0].im*c[1][0].re \
                                          -b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re \
                                          -b[0][2].re*c[1][2].im + b[0][2].im*c[1][2].re;\
                             a[0][2].re =  b[0][0].re*c[2][0].re + b[0][0].im*c[2][0].im \
                                          +b[0][1].re*c[2][1].re + b[0][1].im*c[2][1].im \
                                          +b[0][2].re*c[2][2].re + b[0][2].im*c[2][2].im;\
                             a[0][2].im = -b[0][0].re*c[2][0].im + b[0][0].im*c[2][0].re \
                                          -b[0][1].re*c[2][1].im + b[0][1].im*c[2][1].re \
                                          -b[0][2].re*c[2][2].im + b[0][2].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re + b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[0][1].re + b[1][1].im*c[0][1].im \
                                          +b[1][2].re*c[0][2].re + b[1][2].im*c[0][2].im;\
                             a[1][0].im = -b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          -b[1][1].re*c[0][1].im + b[1][1].im*c[0][1].re \
                                          -b[1][2].re*c[0][2].im + b[1][2].im*c[0][2].re;\
                             a[1][1].re =  b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im \
                                          +b[1][2].re*c[1][2].re + b[1][2].im*c[1][2].im;\
                             a[1][1].im = -b[1][0].re*c[1][0].im + b[1][0].im*c[1][0].re \
                                          -b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re \
                                          -b[1][2].re*c[1][2].im + b[1][2].im*c[1][2].re;\
                             a[1][2].re =  b[1][0].re*c[2][0].re + b[1][0].im*c[2][0].im \
                                          +b[1][1].re*c[2][1].re + b[1][1].im*c[2][1].im \
                                          +b[1][2].re*c[2][2].re + b[1][2].im*c[2][2].im;\
                             a[1][2].im = -b[1][0].re*c[2][0].im + b[1][0].im*c[2][0].re \
                                          -b[1][1].re*c[2][1].im + b[1][1].im*c[2][1].re \
                                          -b[1][2].re*c[2][2].im + b[1][2].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[2][0].re*c[0][0].re + b[2][0].im*c[0][0].im \
                                          +b[2][1].re*c[0][1].re + b[2][1].im*c[0][1].im \
                                          +b[2][2].re*c[0][2].re + b[2][2].im*c[0][2].im;\
                             a[2][0].im = -b[2][0].re*c[0][0].im + b[2][0].im*c[0][0].re \
                                          -b[2][1].re*c[0][1].im + b[2][1].im*c[0][1].re \
                                          -b[2][2].re*c[0][2].im + b[2][2].im*c[0][2].re;\
                             a[2][1].re =  b[2][0].re*c[1][0].re + b[2][0].im*c[1][0].im \
                                          +b[2][1].re*c[1][1].re + b[2][1].im*c[1][1].im \
                                          +b[2][2].re*c[1][2].re + b[2][2].im*c[1][2].im;\
                             a[2][1].im = -b[2][0].re*c[1][0].im + b[2][0].im*c[1][0].re \
                                          -b[2][1].re*c[1][1].im + b[2][1].im*c[1][1].re \
                                          -b[2][2].re*c[1][2].im + b[2][2].im*c[1][2].re;\
                             a[2][2].re =  b[2][0].re*c[2][0].re + b[2][0].im*c[2][0].im \
                                          +b[2][1].re*c[2][1].re + b[2][1].im*c[2][1].im \
                                          +b[2][2].re*c[2][2].re + b[2][2].im*c[2][2].im;\
                             a[2][2].im = -b[2][0].re*c[2][0].im + b[2][0].im*c[2][0].re \
                                          -b[2][1].re*c[2][1].im + b[2][1].im*c[2][1].re \
                                          -b[2][2].re*c[2][2].im + b[2][2].im*c[2][2].re;\
}


/** @def qcd_ADJOINTMUL3x3(a,b,c)
 *  @brief a = adjoint(b)*c, where a,b,c are (complex) 3x3 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\rm conj}({\tt b[k][i]})\,{\tt c[k][j]} \f$
 */
#define qcd_ADJOINTMUL3x3(a,b,c)  {\
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im \
                                          +b[2][0].re*c[2][0].re + b[2][0].im*c[2][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im - b[0][0].im*c[0][0].re \
                                          +b[1][0].re*c[1][0].im - b[1][0].im*c[1][0].re \
                                          +b[2][0].re*c[2][0].im - b[2][0].im*c[2][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re + b[0][0].im*c[0][1].im \
                                          +b[1][0].re*c[1][1].re + b[1][0].im*c[1][1].im \
                                          +b[2][0].re*c[2][1].re + b[2][0].im*c[2][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im - b[0][0].im*c[0][1].re \
                                          +b[1][0].re*c[1][1].im - b[1][0].im*c[1][1].re \
                                          +b[2][0].re*c[2][1].im - b[2][0].im*c[2][1].re;\
                             a[0][2].re =  b[0][0].re*c[0][2].re + b[0][0].im*c[0][2].im \
                                          +b[1][0].re*c[1][2].re + b[1][0].im*c[1][2].im \
                                          +b[2][0].re*c[2][2].re + b[2][0].im*c[2][2].im;\
                             a[0][2].im =  b[0][0].re*c[0][2].im - b[0][0].im*c[0][2].re \
                                          +b[1][0].re*c[1][2].im - b[1][0].im*c[1][2].re \
                                          +b[2][0].re*c[2][2].im - b[2][0].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[0][1].re*c[0][0].re + b[0][1].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re + b[1][1].im*c[1][0].im \
                                          +b[2][1].re*c[2][0].re + b[2][1].im*c[2][0].im;\
                             a[1][0].im =  b[0][1].re*c[0][0].im - b[0][1].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im - b[1][1].im*c[1][0].re \
                                          +b[2][1].re*c[2][0].im - b[2][1].im*c[2][0].re;\
                             a[1][1].re =  b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im \
                                          +b[2][1].re*c[2][1].re + b[2][1].im*c[2][1].im;\
                             a[1][1].im =  b[0][1].re*c[0][1].im - b[0][1].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im - b[1][1].im*c[1][1].re \
                                          +b[2][1].re*c[2][1].im - b[2][1].im*c[2][1].re;\
                             a[1][2].re =  b[0][1].re*c[0][2].re + b[0][1].im*c[0][2].im \
                                          +b[1][1].re*c[1][2].re + b[1][1].im*c[1][2].im \
                                          +b[2][1].re*c[2][2].re + b[2][1].im*c[2][2].im;\
                             a[1][2].im =  b[0][1].re*c[0][2].im - b[0][1].im*c[0][2].re \
                                          +b[1][1].re*c[1][2].im - b[1][1].im*c[1][2].re \
                                          +b[2][1].re*c[2][2].im - b[2][1].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[0][2].re*c[0][0].re + b[0][2].im*c[0][0].im \
                                          +b[1][2].re*c[1][0].re + b[1][2].im*c[1][0].im \
                                          +b[2][2].re*c[2][0].re + b[2][2].im*c[2][0].im;\
                             a[2][0].im =  b[0][2].re*c[0][0].im - b[0][2].im*c[0][0].re \
                                          +b[1][2].re*c[1][0].im - b[1][2].im*c[1][0].re \
                                          +b[2][2].re*c[2][0].im - b[2][2].im*c[2][0].re;\
                             a[2][1].re =  b[0][2].re*c[0][1].re + b[0][2].im*c[0][1].im \
                                          +b[1][2].re*c[1][1].re + b[1][2].im*c[1][1].im \
                                          +b[2][2].re*c[2][1].re + b[2][2].im*c[2][1].im;\
                             a[2][1].im =  b[0][2].re*c[0][1].im - b[0][2].im*c[0][1].re \
                                          +b[1][2].re*c[1][1].im - b[1][2].im*c[1][1].re \
                                          +b[2][2].re*c[2][1].im - b[2][2].im*c[2][1].re;\
                             a[2][2].re =  b[0][2].re*c[0][2].re + b[0][2].im*c[0][2].im \
                                          +b[1][2].re*c[1][2].re + b[1][2].im*c[1][2].im \
                                          +b[2][2].re*c[2][2].re + b[2][2].im*c[2][2].im;\
                             a[2][2].im =  b[0][2].re*c[0][2].im - b[0][2].im*c[0][2].re \
                                          +b[1][2].re*c[1][2].im - b[1][2].im*c[1][2].re \
                                          +b[2][2].re*c[2][2].im - b[2][2].im*c[2][2].re;\
}


/** @def qcd_ADJOINTMULADJOINT3x3(a,b,c)
 *  @brief a = adjoint(b)*adjoint(c), where a,b,c are (complex) 3x3 matrices.
 *  @details \f$ {\tt a[i][j]} = \sum\limits_k {\rm conj}({\tt b[k][i]})\, 
 *                                             {\rm conj}({\tt c[j][k]}) \f$
 */
#define qcd_ADJOINTMULADJOINT3x3(a,b,c)  { \
                             a[0][0].re =  b[0][0].re*c[0][0].re - b[0][0].im*c[0][0].im \
                                          +b[1][0].re*c[0][1].re - b[1][0].im*c[0][1].im \
                                          +b[2][0].re*c[0][2].re - b[2][0].im*c[0][2].im;\
                             a[0][0].im = -b[0][0].re*c[0][0].im - b[0][0].im*c[0][0].re \
                                          -b[1][0].re*c[0][1].im - b[1][0].im*c[0][1].re \
                                          -b[2][0].re*c[0][2].im - b[2][0].im*c[0][2].re;\
                             a[0][1].re =  b[0][0].re*c[1][0].re - b[0][0].im*c[1][0].im \
                                          +b[1][0].re*c[1][1].re - b[1][0].im*c[1][1].im \
                                          +b[2][0].re*c[1][2].re - b[2][0].im*c[1][2].im;\
                             a[0][1].im = -b[0][0].re*c[1][0].im - b[0][0].im*c[1][0].re \
                                          -b[1][0].re*c[1][1].im - b[1][0].im*c[1][1].re \
                                          -b[2][0].re*c[1][2].im - b[2][0].im*c[1][2].re;\
                             a[0][2].re =  b[0][0].re*c[2][0].re - b[0][0].im*c[2][0].im \
                                          +b[1][0].re*c[2][1].re - b[1][0].im*c[2][1].im \
                                          +b[2][0].re*c[2][2].re - b[2][0].im*c[2][2].im;\
                             a[0][2].im = -b[0][0].re*c[2][0].im - b[0][0].im*c[2][0].re \
                                          -b[1][0].re*c[2][1].im - b[1][0].im*c[2][1].re \
                                          -b[2][0].re*c[2][2].im - b[2][0].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[0][1].re*c[0][0].re - b[0][1].im*c[0][0].im \
                                          +b[1][1].re*c[0][1].re - b[1][1].im*c[0][1].im \
                                          +b[2][1].re*c[0][2].re - b[2][1].im*c[0][2].im;\
                             a[1][0].im = -b[0][1].re*c[0][0].im - b[0][1].im*c[0][0].re \
                                          -b[1][1].re*c[0][1].im - b[1][1].im*c[0][1].re \
                                          -b[2][1].re*c[0][2].im - b[2][1].im*c[0][2].re;\
                             a[1][1].re =  b[0][1].re*c[1][0].re - b[0][1].im*c[1][0].im \
                                          +b[1][1].re*c[1][1].re - b[1][1].im*c[1][1].im \
                                          +b[2][1].re*c[1][2].re - b[2][1].im*c[1][2].im;\
                             a[1][1].im = -b[0][1].re*c[1][0].im - b[0][1].im*c[1][0].re \
                                          -b[1][1].re*c[1][1].im - b[1][1].im*c[1][1].re \
                                          -b[2][1].re*c[1][2].im - b[2][1].im*c[1][2].re;\
                             a[1][2].re =  b[0][1].re*c[2][0].re - b[0][1].im*c[2][0].im \
                                          +b[1][1].re*c[2][1].re - b[1][1].im*c[2][1].im \
                                          +b[2][1].re*c[2][2].re - b[2][1].im*c[2][2].im;\
                             a[1][2].im = -b[0][1].re*c[2][0].im - b[0][1].im*c[2][0].re \
                                          -b[1][1].re*c[2][1].im - b[1][1].im*c[2][1].re \
                                          -b[2][1].re*c[2][2].im - b[2][1].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[0][2].re*c[0][0].re - b[0][2].im*c[0][0].im \
                                          +b[1][2].re*c[0][1].re - b[1][2].im*c[0][1].im \
                                          +b[2][2].re*c[0][2].re - b[2][2].im*c[0][2].im;\
                             a[2][0].im = -b[0][2].re*c[0][0].im - b[0][2].im*c[0][0].re \
                                          -b[1][2].re*c[0][1].im - b[1][2].im*c[0][1].re \
                                          -b[2][2].re*c[0][2].im - b[2][2].im*c[0][2].re;\
                             a[2][1].re =  b[0][2].re*c[1][0].re - b[0][2].im*c[1][0].im \
                                          +b[1][2].re*c[1][1].re - b[1][2].im*c[1][1].im \
                                          +b[2][2].re*c[1][2].re - b[2][2].im*c[1][2].im;\
                             a[2][1].im = -b[0][2].re*c[1][0].im - b[0][2].im*c[1][0].re \
                                          -b[1][2].re*c[1][1].im - b[1][2].im*c[1][1].re \
                                          -b[2][2].re*c[1][2].im - b[2][2].im*c[1][2].re;\
                             a[2][2].re =  b[0][2].re*c[2][0].re - b[0][2].im*c[2][0].im \
                                          +b[1][2].re*c[2][1].re - b[1][2].im*c[2][1].im \
                                          +b[2][2].re*c[2][2].re - b[2][2].im*c[2][2].im;\
                             a[2][2].im = -b[0][2].re*c[2][0].im - b[0][2].im*c[2][0].re \
                                          -b[1][2].re*c[2][1].im - b[1][2].im*c[2][1].re \
                                          -b[2][2].re*c[2][2].im - b[2][2].im*c[2][2].re;\
}

/** @def qcd_PROJ_ANTIHERMITIAN_TRACELESS3x3(b,a)
 *  @brief Projects a onto the space of antihermitian traceless 3x3 matrices
 *  @details \f$ {\tt tmp} = \frac{1}{2}[a-a^\dagger]\quad  \f$ (antihermitian)\n
 *           \f$ b = {\tt tmp} - \frac{{\rm tr}({\tt tmp})}{3}  {\mathds 1}\quad \f$
 *           (antihermitian, traceless)
 *  @todo untested
 */
#define qcd_PROJ_ANTIHERMITIAN_TRACELESS3x3(b,a) \
   { b[2][2].im = (a[0][0].im+a[1][1].im+a[2][2].im)*0.3333333333333333; \
     b[0][0].re = 0;                              b[0][0].im = a[0][0].im-b[2][2].im; \
     b[0][1].re = 0.5*(a[0][1].re - a[1][0].re);  b[0][1].im = 0.5*(a[0][1].im + a[1][0].im); \
     b[0][2].re = 0.5*(a[0][2].re - a[2][0].re);  b[0][2].im = 0.5*(a[0][2].im + a[2][0].im); \
     b[1][0].re = -b[0][1].re;                    b[1][0].im = b[0][1].im; \
     b[1][1].re = 0;                              b[1][1].im = a[1][1].im-b[2][2].im; \
     b[1][2].re = 0.5*(a[1][2].re - a[2][1].re);  b[1][2].im = 0.5*(a[1][2].im + a[2][1].im); \
     b[2][0].re = -b[0][2].re;                    b[2][0].im = b[0][2].im; \
     b[2][1].re = -b[1][2].re;                    b[2][1].im = b[1][2].im; \
     b[2][2].re = 0;                              b[2][2].im = a[2][2].im-b[2][2].im; \
  }


/** @def qcd_PROJ_ANTIHERMITIAN_TRACELESS_ECO3x3(b,a)
 *  @brief Projects a onto the space of antihermitian traceless 3x3 matrices
 *  @details \f$ {\tt tmp} = \frac{1}{2}[a-a^\dagger]\quad  \f$ (antihermitian)\n
 *           \f$ b = {\tt tmp} - \frac{{\rm tr}({\tt tmp})}{3}  {\mathds 1}\quad \f$
 *           (antihermitian, traceless)\n
 *  only elements Im(00) 01 02 Im(11) 12 Im(22) are overwritten
 *  @todo untested
 */
#define qcd_PROJ_ANTIHERMITIAN_TRACELESS_ECO3x3(b,a) \
  {  b[2][2].im = (a[0][0].im+a[1][1].im+a[2][2].im)*0.3333333333333333; \
                                                  b[0][0].im = a[0][0].im-b[2][2].im; \
     b[0][1].re = 0.5*(a[0][1].re - a[1][0].re);  b[0][1].im = 0.5*(a[0][1].im + a[1][0].im); \
     b[0][2].re = 0.5*(a[0][2].re - a[2][0].re);  b[0][2].im = 0.5*(a[0][2].im + a[2][0].im); \
                                                  b[1][1].im = a[1][1].im-b[2][2].im; \
     b[1][2].re = 0.5*(a[1][2].re - a[2][1].re);  b[1][2].im = 0.5*(a[1][2].im + a[2][1].im); \
                                                  b[2][2].im = a[2][2].im-b[2][2].im; \
  }

/**
 * @def qcd_SU3TRACER(u)
 * @brief Real part of the trace of u.
 * @details \f$ {\tt result} = \sum\limits_k {\rm Re}({\tt u[k][k]}) \f$
 */
#define qcd_SU3TRACER(u) ( (double) (u[0][0].re + u[1][1].re + u[2][2].re))

/**
 * @def qcd_APPLY_U(u,x,y)
 * @brief Applies gauge field u to spinor y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[dirac][color] and Y[dirac][color] are stored such that
 \verbatim
           x[0] = real(X[0][0])
           x[1] = imag(X[0][0])
           x[2] = real(X[0][1])
           ...
          x[23] = imag(X[3][2])
 \endverbatim
 */
#define qcd_APPLY_U(u,x,y)  \
    { x[0] = u[0] *y[0]  - u[1] *y[1]  + u[2] *y[2]  - u[3] *y[3]  + u[4] *y[4]  - u[5] *y[5];     \
      x[1] = u[0] *y[1]  + u[1] *y[0]  + u[2] *y[3]  + u[3] *y[2]  + u[4] *y[5]  + u[5] *y[4];     \
      x[2] = u[6] *y[0]  - u[7] *y[1]  + u[8] *y[2]  - u[9] *y[3]  + u[10]*y[4]  - u[11]*y[5];     \
      x[3] = u[6] *y[1]  + u[7] *y[0]  + u[8] *y[3]  + u[9] *y[2]  + u[10]*y[5]  + u[11]*y[4];     \
      x[4] = u[12]*y[0]  - u[13]*y[1]  + u[14]*y[2]  - u[15]*y[3]  + u[16]*y[4]  - u[17]*y[5];     \
      x[5] = u[12]*y[1]  + u[13]*y[0]  + u[14]*y[3]  + u[15]*y[2]  + u[16]*y[5]  + u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  - u[1] *y[7]  + u[2] *y[8]  - u[3] *y[9]  + u[4] *y[10] - u[5] *y[11];    \
      x[7] = u[0] *y[7]  + u[1] *y[6]  + u[2] *y[9]  + u[3] *y[8]  + u[4] *y[11] + u[5] *y[10];    \
      x[8] = u[6] *y[6]  - u[7] *y[7]  + u[8] *y[8]  - u[9] *y[9]  + u[10]*y[10] - u[11]*y[11];    \
      x[9] = u[6] *y[7]  + u[7] *y[6]  + u[8] *y[9]  + u[9] *y[8]  + u[10]*y[11] + u[11]*y[10];    \
      x[10]= u[12]*y[6]  - u[13]*y[7]  + u[14]*y[8]  - u[15]*y[9]  + u[16]*y[10] - u[17]*y[11];    \
      x[11]= u[12]*y[7]  + u[13]*y[6]  + u[14]*y[9]  + u[15]*y[8]  + u[16]*y[11] + u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] - u[1] *y[13] + u[2] *y[14] - u[3] *y[15] + u[4] *y[16] - u[5] *y[17];    \
      x[13]= u[0] *y[13] + u[1] *y[12] + u[2] *y[15] + u[3] *y[14] + u[4] *y[17] + u[5] *y[16];    \
      x[14]= u[6] *y[12] - u[7] *y[13] + u[8] *y[14] - u[9] *y[15] + u[10]*y[16] - u[11]*y[17];    \
      x[15]= u[6] *y[13] + u[7] *y[12] + u[8] *y[15] + u[9] *y[14] + u[10]*y[17] + u[11]*y[16];    \
      x[16]= u[12]*y[12] - u[13]*y[13] + u[14]*y[14] - u[15]*y[15] + u[16]*y[16] - u[17]*y[17];    \
      x[17]= u[12]*y[13] + u[13]*y[12] + u[14]*y[15] + u[15]*y[14] + u[16]*y[17] + u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] - u[1] *y[19] + u[2] *y[20] - u[3] *y[21] + u[4] *y[22] - u[5] *y[23];    \
      x[19]= u[0] *y[19] + u[1] *y[18] + u[2] *y[21] + u[3] *y[20] + u[4] *y[23] + u[5] *y[22];    \
      x[20]= u[6] *y[18] - u[7] *y[19] + u[8] *y[20] - u[9] *y[21] + u[10]*y[22] - u[11]*y[23];    \
      x[21]= u[6] *y[19] + u[7] *y[18] + u[8] *y[21] + u[9] *y[20] + u[10]*y[23] + u[11]*y[22];    \
      x[22]= u[12]*y[18] - u[13]*y[19] + u[14]*y[20] - u[15]*y[21] + u[16]*y[22] - u[17]*y[23];    \
      x[23]= u[12]*y[19] + u[13]*y[18] + u[14]*y[21] + u[15]*y[20] + u[16]*y[23] + u[17]*y[22];    \
   }

/**
 * @def qcd_APPLY_U_HALFSPINOR(u,x,y)
 * @brief Applies gauge field u to half-spinor y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[dirac][color] and Y[dirac][color] are stored such that
 \verbatim
           x[0] = real(X[0][0])
           x[1] = imag(X[0][0])
           x[2] = real(X[0][1])
           ...
          x[11] = imag(X[1][2])
 \endverbatim
 */
#define qcd_APPLY_U_HALFSPINOR(u,x,y)  \
    { x[0] = u[0] *y[0]  - u[1] *y[1]  + u[2] *y[2]  - u[3] *y[3]  + u[4] *y[4]  - u[5] *y[5];     \
      x[1] = u[0] *y[1]  + u[1] *y[0]  + u[2] *y[3]  + u[3] *y[2]  + u[4] *y[5]  + u[5] *y[4];     \
      x[2] = u[6] *y[0]  - u[7] *y[1]  + u[8] *y[2]  - u[9] *y[3]  + u[10]*y[4]  - u[11]*y[5];     \
      x[3] = u[6] *y[1]  + u[7] *y[0]  + u[8] *y[3]  + u[9] *y[2]  + u[10]*y[5]  + u[11]*y[4];     \
      x[4] = u[12]*y[0]  - u[13]*y[1]  + u[14]*y[2]  - u[15]*y[3]  + u[16]*y[4]  - u[17]*y[5];     \
      x[5] = u[12]*y[1]  + u[13]*y[0]  + u[14]*y[3]  + u[15]*y[2]  + u[16]*y[5]  + u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  - u[1] *y[7]  + u[2] *y[8]  - u[3] *y[9]  + u[4] *y[10] - u[5] *y[11];    \
      x[7] = u[0] *y[7]  + u[1] *y[6]  + u[2] *y[9]  + u[3] *y[8]  + u[4] *y[11] + u[5] *y[10];    \
      x[8] = u[6] *y[6]  - u[7] *y[7]  + u[8] *y[8]  - u[9] *y[9]  + u[10]*y[10] - u[11]*y[11];    \
      x[9] = u[6] *y[7]  + u[7] *y[6]  + u[8] *y[9]  + u[9] *y[8]  + u[10]*y[11] + u[11]*y[10];    \
      x[10]= u[12]*y[6]  - u[13]*y[7]  + u[14]*y[8]  - u[15]*y[9]  + u[16]*y[10] - u[17]*y[11];    \
      x[11]= u[12]*y[7]  + u[13]*y[6]  + u[14]*y[9]  + u[15]*y[8]  + u[16]*y[11] + u[17]*y[10];    \
   }

/**
 * @def qcd_APPLY_U_SPINORCOMPONENT(u,x,y)
 * @brief Applies gauge field u to color-vector y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[color] and Y[color] are stored such that
 \verbatim
           x[0] = real(X[0])
           x[1] = imag(X[0])
           x[2] = real(X[1])
           ...
           x[5] = imag(X[2])
 \endverbatim
 */
#define qcd_APPLY_U_SPINORCOMPONENT(u,x,y)  \
    { x[0] = u[0] *y[0]  - u[1] *y[1]  + u[2] *y[2]  - u[3] *y[3]  + u[4] *y[4]  - u[5] *y[5];     \
      x[1] = u[0] *y[1]  + u[1] *y[0]  + u[2] *y[3]  + u[3] *y[2]  + u[4] *y[5]  + u[5] *y[4];     \
      x[2] = u[6] *y[0]  - u[7] *y[1]  + u[8] *y[2]  - u[9] *y[3]  + u[10]*y[4]  - u[11]*y[5];     \
      x[3] = u[6] *y[1]  + u[7] *y[0]  + u[8] *y[3]  + u[9] *y[2]  + u[10]*y[5]  + u[11]*y[4];     \
      x[4] = u[12]*y[0]  - u[13]*y[1]  + u[14]*y[2]  - u[15]*y[3]  + u[16]*y[4]  - u[17]*y[5];     \
      x[5] = u[12]*y[1]  + u[13]*y[0]  + u[14]*y[3]  + u[15]*y[2]  + u[16]*y[5]  + u[17]*y[4];     \
   }

/**
 * @def qcd_APPLY_U_DAGGER(u,x,y)
 * @brief Applies gauge field adj(u) to spinor y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[dirac][color] and Y[dirac][color] are stored such that
 \verbatim
           x[0] = real(X[0][0])
           x[1] = imag(X[0][0])
           x[2] = real(X[0][1])
           ...
          x[23] = imag(X[3][2])
 \endverbatim
 */
#define qcd_APPLY_U_DAGGER(u,x,y)  \
    { x[0] = u[0] *y[0]  + u[1] *y[1]  + u[6] *y[2]  + u[7] *y[3]  + u[12]*y[4]  + u[13]*y[5];     \
      x[1] = u[0] *y[1]  - u[1] *y[0]  + u[6] *y[3]  - u[7] *y[2]  + u[12]*y[5]  - u[13]*y[4];     \
      x[2] = u[2] *y[0]  + u[3] *y[1]  + u[8] *y[2]  + u[9] *y[3]  + u[14]*y[4]  + u[15]*y[5];     \
      x[3] = u[2] *y[1]  - u[3] *y[0]  + u[8] *y[3]  - u[9] *y[2]  + u[14]*y[5]  - u[15]*y[4];     \
      x[4] = u[4] *y[0]  + u[5] *y[1]  + u[10]*y[2]  + u[11]*y[3]  + u[16]*y[4]  + u[17]*y[5];     \
      x[5] = u[4] *y[1]  - u[5] *y[0]  + u[10]*y[3]  - u[11]*y[2]  + u[16]*y[5]  - u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  + u[1] *y[7]  + u[6] *y[8]  + u[7] *y[9]  + u[12]*y[10] + u[13]*y[11];    \
      x[7] = u[0] *y[7]  - u[1] *y[6]  + u[6] *y[9]  - u[7] *y[8]  + u[12]*y[11] - u[13]*y[10];    \
      x[8] = u[2] *y[6]  + u[3] *y[7]  + u[8] *y[8]  + u[9] *y[9]  + u[14]*y[10] + u[15]*y[11];    \
      x[9] = u[2] *y[7]  - u[3] *y[6]  + u[8] *y[9]  - u[9] *y[8]  + u[14]*y[11] - u[15]*y[10];    \
      x[10]= u[4] *y[6]  + u[5] *y[7]  + u[10]*y[8]  + u[11]*y[9]  + u[16]*y[10] + u[17]*y[11];    \
      x[11]= u[4] *y[7]  - u[5] *y[6]  + u[10]*y[9]  - u[11]*y[8]  + u[16]*y[11] - u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] + u[1] *y[13] + u[6] *y[14] + u[7] *y[15] + u[12]*y[16] + u[13]*y[17];    \
      x[13]= u[0] *y[13] - u[1] *y[12] + u[6] *y[15] - u[7] *y[14] + u[12]*y[17] - u[13]*y[16];    \
      x[14]= u[2] *y[12] + u[3] *y[13] + u[8] *y[14] + u[9] *y[15] + u[14]*y[16] + u[15]*y[17];    \
      x[15]= u[2] *y[13] - u[3] *y[12] + u[8] *y[15] - u[9] *y[14] + u[14]*y[17] - u[15]*y[16];    \
      x[16]= u[4] *y[12] + u[5] *y[13] + u[10]*y[14] + u[11]*y[15] + u[16]*y[16] + u[17]*y[17];    \
      x[17]= u[4] *y[13] - u[5] *y[12] + u[10]*y[15] - u[11]*y[14] + u[16]*y[17] - u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] + u[1] *y[19] + u[6] *y[20] + u[7] *y[21] + u[12]*y[22] + u[13]*y[23];    \
      x[19]= u[0] *y[19] - u[1] *y[18] + u[6] *y[21] - u[7] *y[20] + u[12]*y[23] - u[13]*y[22];    \
      x[20]= u[2] *y[18] + u[3] *y[19] + u[8] *y[20] + u[9] *y[21] + u[14]*y[22] + u[15]*y[23];    \
      x[21]= u[2] *y[19] - u[3] *y[18] + u[8] *y[21] - u[9] *y[20] + u[14]*y[23] - u[15]*y[22];    \
      x[22]= u[4] *y[18] + u[5] *y[19] + u[10]*y[20] + u[11]*y[21] + u[16]*y[22] + u[17]*y[23];    \
      x[23]= u[4] *y[19] - u[5] *y[18] + u[10]*y[21] - u[11]*y[20] + u[16]*y[23] - u[17]*y[22];    \
   }

/**
 * @def qcd_APPLY_U_DAGGER_HALFSPINOR(u,x,y)
 * @brief Applies gauge field adj(u) to half-spinor y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[dirac][color] and Y[dirac][color] are stored such that
 \verbatim
           x[0] = real(X[0][0])
           x[1] = imag(X[0][0])
           x[2] = real(X[0][1])
           ...
          x[11] = imag(X[1][2])
 \endverbatim
 */
#define qcd_APPLY_U_DAGGER_HALFSPINOR(u,x,y)  \
    { x[0] = u[0] *y[0]  + u[1] *y[1]  + u[6] *y[2]  + u[7] *y[3]  + u[12]*y[4]  + u[13]*y[5];     \
      x[1] = u[0] *y[1]  - u[1] *y[0]  + u[6] *y[3]  - u[7] *y[2]  + u[12]*y[5]  - u[13]*y[4];     \
      x[2] = u[2] *y[0]  + u[3] *y[1]  + u[8] *y[2]  + u[9] *y[3]  + u[14]*y[4]  + u[15]*y[5];     \
      x[3] = u[2] *y[1]  - u[3] *y[0]  + u[8] *y[3]  - u[9] *y[2]  + u[14]*y[5]  - u[15]*y[4];     \
      x[4] = u[4] *y[0]  + u[5] *y[1]  + u[10]*y[2]  + u[11]*y[3]  + u[16]*y[4]  + u[17]*y[5];     \
      x[5] = u[4] *y[1]  - u[5] *y[0]  + u[10]*y[3]  - u[11]*y[2]  + u[16]*y[5]  - u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  + u[1] *y[7]  + u[6] *y[8]  + u[7] *y[9]  + u[12]*y[10] + u[13]*y[11];    \
      x[7] = u[0] *y[7]  - u[1] *y[6]  + u[6] *y[9]  - u[7] *y[8]  + u[12]*y[11] - u[13]*y[10];    \
      x[8] = u[2] *y[6]  + u[3] *y[7]  + u[8] *y[8]  + u[9] *y[9]  + u[14]*y[10] + u[15]*y[11];    \
      x[9] = u[2] *y[7]  - u[3] *y[6]  + u[8] *y[9]  - u[9] *y[8]  + u[14]*y[11] - u[15]*y[10];    \
      x[10]= u[4] *y[6]  + u[5] *y[7]  + u[10]*y[8]  + u[11]*y[9]  + u[16]*y[10] + u[17]*y[11];    \
      x[11]= u[4] *y[7]  - u[5] *y[6]  + u[10]*y[9]  - u[11]*y[8]  + u[16]*y[11] - u[17]*y[10];    \
   }

/**
 * @def qcd_APPLY_U_DAGGER_SPINORCOMPONENT(u,x,y)
 * @brief Applies gauge field adj(u) to color-vector y and stores the result in x.
 * @details The complex 3x3 matrix U is stored such that\n
 \verbatim
           u[0] = real(U[0][0]),
           u[1] = imag(U[0][0]),
           u[2] = real(U[0][1]),
           ...
           u[17]= imag(U[2][2]),
 \endverbatim
 * and X[color] and Y[color] are stored such that
 \verbatim
           x[0] = real(X[0])
           x[1] = imag(X[0])
           x[2] = real(X[1])
           ...
           x[5] = imag(X[2])
 \endverbatim
 */
#define qcd_APPLY_U_DAGGER_SPINORCOMPONENT(u,x,y)  \
    { x[0] = u[0] *y[0]  + u[1] *y[1]  + u[6] *y[2]  + u[7] *y[3]  + u[12]*y[4]  + u[13]*y[5];     \
      x[1] = u[0] *y[1]  - u[1] *y[0]  + u[6] *y[3]  - u[7] *y[2]  + u[12]*y[5]  - u[13]*y[4];     \
      x[2] = u[2] *y[0]  + u[3] *y[1]  + u[8] *y[2]  + u[9] *y[3]  + u[14]*y[4]  + u[15]*y[5];     \
      x[3] = u[2] *y[1]  - u[3] *y[0]  + u[8] *y[3]  - u[9] *y[2]  + u[14]*y[5]  - u[15]*y[4];     \
      x[4] = u[4] *y[0]  + u[5] *y[1]  + u[10]*y[2]  + u[11]*y[3]  + u[16]*y[4]  + u[17]*y[5];     \
      x[5] = u[4] *y[1]  - u[5] *y[0]  + u[10]*y[3]  - u[11]*y[2]  + u[16]*y[5]  - u[17]*y[4];     \
   }

/******************************************************************************
 * prototypes 
 ******************************************************************************/

/**
 * @fn void qcd_zero3x3(qcd_complex_16 a[3][3])
 * @brief Sets the complex 3x3 matrix a to zero.
 * @param[in,out] a pointer to a complex 3x3 matrix
 */
void qcd_zero3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_unit3x3(qcd_complex_16 a[3][3])
 * @brief Sets the complex 3x3 matrix to a unit matrix.
 * @param[in,out] a pointer to a complex 3x3 matrix
 */
void qcd_unit3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_minusUnit3x3(qcd_complex_16 a[3][3])
 * @brief Sets the complex 3x3 matrix to minus 1.
 * @details Where 1 is a complex unit matrix.
 * @param[in,out] a pointer to a complex 3x3 matrix
 */
void qcd_minusUnit3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_sub3x3(qcd_complex_16 difference[3][3], qcd_complex_16 minuend[3][3], qcd_complex_16 subtrahend[3][3]);
 * @brief Subtracts two complex 3x3 matrices from each other.
 * @details difference = minuend - subtrahend
 * @param[out] difference pointer to a complex 3x3 matrix
 * @param[in] minuend pointer to a complex 3x3 matrix
 * @param[in] subtrahend pointer to a complex 3x3 matrix
 */
void qcd_sub3x3(qcd_complex_16 difference[3][3], qcd_complex_16 minuend[3][3], qcd_complex_16 subtrahend[3][3]);

/**
 * @fn void qcd_add3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);
 * @brief Sums two complex 3x3 matrices.
 * @details sum = summand1 + summand2\n
 * Two or three of the parameters may be identical.
 * @param[out] sum pointer to a complex 3x3 matrix
 * @param[in] summand1 pointer to a complex 3x3 matrix
 * @param[in] summand2 pointer to a complex 3x3 matrix
 */
void qcd_add3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);

/**
 * @fn void qcd_addAdjoint3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);
 * @brief Adds a 3x3 matrix to the adjoint of another.
 * @details sum = summand1 + adj(summand2)\n
 * Two or three of the parameters may be identical.
 * @param[out] sum pointer to a complex 3x3 matrix
 * @param[in] summand1 pointer to a complex 3x3 matrix
 * @param[in] summand2 pointer to a complex 3x3 matrix
 */
void qcd_addAdjoint3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);

/**
 * @fn void qcd_mul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
 * @brief Multiplies two complex 3x3 matrices
 * @details product = factor1 * factor2\n
 * factor1 and factor2 may be identical.
 * @param[out] product pointer to a complex 3x3 matrix
 * @param[in] factor1 pointer to a complex 3x3 matrix
 * @param[in] factor2 pointer to a complex 3x3 matrix
 */
void qcd_mul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);

/**
 * @fn void qcd_mulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
 * @brief Multiply a 3x3 matrix with the adjoint of another
 * @details product = factor1 * adj(factor2)
 * @param[out] product pointer to a complex 3x3 matrix
 * @param[in] factor1 pointer to a complex 3x3 matrix
 * @param[in] factor2 pointer to a complex 3x3 matrix
 */
void qcd_mulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);

/**
 * @fn void qcd_adjointMul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
 * @brief Multiply the adjoint of a 3x3 matrix with another 3x3 matrix
 * @details product = adj(factor1) * factor2
 * @param[out] product pointer to a complex 3x3 matrix
 * @param[in] factor1 pointer to a complex 3x3 matrix
 * @param[in] factor2 pointer to a complex 3x3 matrix
 */
void qcd_adjointMul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);

/**
 * @fn void qcd_adjointMulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
 * @brief Multiply the adjoint of a 3x3 matrix with the adjoint of another
 * @details product = adj(factor1) * adj(factor2)
 * @param[out] product pointer to a complex 3x3 matrix
 * @param[in] factor1 pointer to a complex 3x3 matrix
 * @param[in] factor2 pointer to a complex 3x3 matrix
 */
void qcd_adjointMulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);


/**
 * @fn void qcd_scale3x3(qcd_complex_16 a[3][3], double r)
 * @brief Scales a complex 3x3 matrix by a real number
 * @details a = a*r
 * @param[in,out] a pointer to a complex 3x3 matrix
 * @param[in] r real rescaling factor
 */
void qcd_scale3x3(qcd_complex_16 a[3][3], double r);

/**
 * @fn void qcd_cScale3x3(qcd_complex_16 a[3][3], qcd_complex_16 c)
 * @brief Scales a complex 3x3 matrix by a complex number
 * @details a = a*c
 * @param[in,out] a pointer to a complex 3x3 matrix
 * @param[in] c complex rescaling factor
 */
void qcd_cScale3x3(qcd_complex_16 a[3][3], qcd_complex_16 c);

/**
 * @fn qcd_complex_16 qcd_trace3x3(qcd_complex_16 a[3][3])
 * @brief Computes the trace of a complex 3x3 matrix.
 * @details trace = a[0][0]+a[1][1]+a[2][2]
 * @param[in] a pointer to a complex 3x3 matrix
 * @return the (complex) trace of the 3x3 matrix
 */
qcd_complex_16 qcd_trace3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_copy3x3(qcd_complex_16 dest[3][3], qcd_complex_16 src[3][3])
 * @brief Copies a complex 3x3 matrix
 * @details dest[i][j] = src[i][j] for i,j=0,1,2
 * @param[out] dest a pointer to a complex 3x3 matrix
 * @param[in] src pointer to a complex 3x3 matrix
 */
void qcd_copy3x3(qcd_complex_16 dest[3][3], qcd_complex_16 src[3][3]);

/**
 * @fn double qcd_norm3x3(qcd_complex_16 a[3][3]);
 * @brief Computes the Froebenius norm of a complex 3x3 matrix
 * @details nrm = sqrt(a^dag * a)
 * @param[in] a pointer to a complex 3x3 matrix
 * @return the real, positive norm
 */
double qcd_norm3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_adjoint3x3(qcd_complex_16 a[3][3])
 * @brief Computes the adjoint of a complex 3x3 matrix
 * @details a[i][j] = conj(a[j][i]) for i,j=0,1,2
 * @param[in,out] a pointer to a complex 3x3 matrix
 */
void qcd_adjoint3x3(qcd_complex_16 a[3][3]);

/**
 * @fn qcd_complex_16 qcd_det3x3(qcd_complex_16 a[3][3])
 * @brief Computes the determinant of a complex 3x3 matrix
 * @details The (complex) determinant is computed for a general complex 
 *          3x3 matrix a.
 * @param[in] a pointer to a complex 3x3 matrix
 * @return complex determinant of a
 */
qcd_complex_16 qcd_det3x3(qcd_complex_16 a[3][3]);

/**
 * @fn void qcd_inv3x3(qcd_complex_16 ai[3][3], qcd_complex_16 a[3][3])
 * @brief Computes the inverse of a complex 3x3 matrix
 * @details ai is the inverse of the complex 3x3 matrix a. If a is singular,
 *          the inverse will contain infinite entries caused by divisions by zero.
 * @param[in] a pointer to a complex 3x3 matrix
 * @param[out] ai pointer to a complex 3x3 matrix
 */
void qcd_inv3x3(qcd_complex_16 ai[3][3], qcd_complex_16 a[3][3]);
/**
 * @fn void qcd_projectSU33x3(qcd_complex_16 g[3][3])
 * @brief Projects a 3x3 matrix to SU(3).
 * @details If M is an arbitrary 3x3 matrix it can be factorized as the product
 *          of a hermitian and a unitary matrix, M = U H. We will calculate U 
 *          as the projection of M to SU(3). Now M^dag M = H^2. By
 *          going to diagonal basis for H^2 we calculate \n
 *          H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^dag \n
 *          where v is the 3x3 matrix whos columns are the eigenvectors of H^2 
 *          and L is the diagonal representation of H^2. Thus: U = M H^-1 \n
 *          [1] Phys. Lett. B307, 375-382, (1993)
 * @author Giannis Koutsou
 * @param[in,out] g  pointer to a complex 3x3 matrix
 */
void qcd_projectSU33x3(qcd_complex_16 g[3][3]);

#endif
