/** @file qcd_COps.h
 *  @brief header file for qcd_COps.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2022
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_C
#define H_QCD_C 1
/* preprocessor macros, use carefully */
/** @def qcd_CONJ(x)
 *  @brief Complex conjugate of x.
 *  @details \f$ {\tt result} = {\rm Re}(x) - i\,{\rm Im}(x) \f$
 */
#define qcd_CONJ(x)    ( (qcd_complex_16) {x.re,-x.im} )

/** @def qcd_CMUL(x,y)
 *  @brief Complex multiplication of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)-{\rm Im}(x){\rm Im}(y) 
 *                           +i[{\rm Re}(x){\rm Im}(y)+{\rm Im}(x){\rm Re}(y)] \f$
 */
#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )

/** @def qcd_CMULR(x,y)
 *  @brief Real part of the complex multiplication of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)-{\rm Im}(x){\rm Im}(y) \f$
 */
#define qcd_CMULR(x,y) ( (double) (x.re * y.re - x.im * y.im) )

/** @def qcd_CMULI(x,y)
 *  @brief Imaginary part of the complex multiplication of x and y.
 *  @details \f$ {\tt result} = [{\rm Re}(x){\rm Im}(y)+{\rm Im}(x){\rm Re}(y)] \f$
 */
#define qcd_CMULI(x,y) ( (double) (x.re * y.im + x.im * y.re) )

/** @def qcd_CCONJMUL(x,y)
 *  @brief Complex multiplication of conj(x) and y.
 *  @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y) 
 *                           +i[{\rm Re}(x){\rm Im}(y)-{\rm Im}(x){\rm Re}(y)]\f$
 */
#define qcd_CCONJMUL(x,y)  ( (qcd_complex_16) {x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re } )

/** @def qcd_CADD(x,y)
 *  @brief Complex sum of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x)+{\rm Re}(y) + i[{\rm Im}(x)+{\rm Im}(y)] \f$
 */
#define qcd_CADD(x,y)  ( (qcd_complex_16) {x.re+y.re, x.im+y.im})

/** @def qcd_CADDR(x,y)
 *  @brief Real part of complex sum of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x)+{\rm Re}(y) \f$
 */
#define qcd_CADDR(x,y) ( (double) (x.re+y.re))

/** @def qcd_CADDI(x,y)
 *  @brief Imaginary part of complex sum of x and y.
 *  @details \f$ {\tt result} = {\rm Im}(x)+{\rm Im}(y) \f$
 */
#define qcd_CADDI(x,y) ( (double) (x.im+y.im))

/** @def qcd_CSUB(x,y)
 *  @brief Complex difference of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x)-{\rm Re}(y) + i[{\rm Im}(x)-{\rm Im}(y)] \f$
 */
#define qcd_CSUB(x,y)  ( (qcd_complex_16) {x.re-y.re, x.im-y.im})

/** @def qcd_CSUBR(x,y)
 *  @brief Real part of complex difference of x and y.
 *  @details \f$ {\tt result} = {\rm Re}(x)-{\rm Re}(y) \f$
 */
#define qcd_CSUBR(x,y) ( (double) (x.re-y.re))

/** @def qcd_CSUBI(x,y)
 *  @brief Imaginary part of complex difference of x and y.
 *  @details \f$ {\tt result} = {\rm Im}(x)-{\rm Im}(y) \f$
 */
#define qcd_CSUBI(x,y) ( (double) (x.im-y.im))

/** @def qcd_CSCALE(x,a)
 *  @brief Multiplication of complex x by real a.
 *  @details \f$ {\tt result} = a\,{\rm Re}(x) + i\,a\,{\rm Im}(x) \f$
 */
#define qcd_CSCALE(x,a)  ( (qcd_complex_16) {x.re*(a), x.im*(a)})

/** @def qcd_CSCALER(x,a)
 *  @brief Real part of multiplication of complex x by real a.
 *  @details \f$ {\tt result} = a\, {\rm Re}(x) \f$
 */
#define qcd_CSCALER(x,a) ( (double) x.re*(a))

/** @def qcd_CSCALEI(x,a)
 *  @brief Imaginary part of multiplication of complex x by real a.
 *  @details \f$ {\tt result} = a\,{\rm Im}(x) \f$
 */
#define qcd_CSCALEI(x,a) ( (double) x.im*(a))

/** @def qcd_CDIV(x,y)
 *  @brief Complex division of x by y.
 *  @details \f$ {\tt result} = \frac{  {\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y) 
 *                                   +i[{\rm Im}(x){\rm Re}(y)-{\rm Re}(x){\rm Im}(y)]}
 *                                   {{\rm Re}(y)^2 + {\rm Im}(y)^2)} \f$
 */
#define qcd_CDIV(x,y)  ( (qcd_complex_16) {(x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im), (x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) } )

/** @def qcd_CDIVR(x,y)
 *  @brief Real part of complex division of x by y.
 *  @details \f$ {\tt result} = \frac{{\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y)}{{\rm Re}(y)^2 + {\rm Im}(y)^2} \f$
 */
#define qcd_CDIVR(x,y) ( (double)     ((x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im) ))

/** @def qcd_CDIVI(x,y)
 *  @brief Imaginary part of complex division of x by y.
 *  @details \f$ {\tt result} = \frac{[{\rm Im}(x){\rm Re}(y)-{\rm Re}(x){\rm Im}(y)]}
 *                                   {{\rm Re}(y)^2 + {\rm Im}(y)^2} \f$
 */
#define qcd_CDIVI(x,y) ( (double)     ((x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) ))

/** @def qcd_CNORM(x)
 *  @brief Complex norm of x
 *  @details \f$ {\tt result} = \sqrt{ {\rm Re}(x)^2 + {\rm Im}(x)^2 } \f$
 */
#define qcd_CNORM(x)    ( (double) sqrt(x.re * x.re + x.im * x.im))

/** @def qcd_CNORMSQUARED(x)
 *  @brief Square of complex norm of x
 *  @details \f$ {\tt result} = {\rm Re}(x)^2 + {\rm Im}(x)^2 \f$
 */
#define qcd_CNORMSQUARED(x) ( (double) (x.re * x.re + x.im * x.im))

/** @def qcd_CARG(x)
 *  @brief The angle &phi; in the polar representation of the complex number x.
 *  @details If \f$ x = |x| \, e^{i\,\phi}\f$, then \n
 *           \f$ \phi = {\rm atan2}({\rm Im}(x),{\rm Re}(x)) \f$.\n
 *           (see <a href=https://en.wikipedia.org/wiki/Atan2>atan2</a>)
 */
#define qcd_CARG(x)     ( (double) atan2(x.im,x.re))

/** @def qcd_CPOW(x,a)
 *  @brief Complex number x to the power of real a
 *  @details If \f$ x   = |x|   \, e^{i\, \phi}\f$, then \n
 *              \f$ x^a = |x|^a \, e^{i\,a\,\phi}\f$
 */
#define qcd_CPOW(x,a)  ( (qcd_complex_16) {pow(qcd_CNORM(x),(a))*cos(qcd_CARG(x)*(a)), pow(qcd_CNORM(x),(a))*sin(qcd_CARG(x)*(a))})

/** @def qcd_CPOWR(x,a)
 *  @brief Real part of the complex number x raised to the power of real a
 *  @details If \f$ x   = |x| \, e^{i\, \phi}\f$, then \n
 *              \f${\rm Re}(x^a) = |x|^a \, \cos(a\, \phi)\f$
 */
#define qcd_CPOWR(x,a) ( (double) pow(qcd_CNORM(x),(a))*cos(qcd_CARG(x)*(a)))

/** @def qcd_CPOWI(x,a)
 *  @brief Imaginary part of the complex number x raised to the power of real a
 *  @details If \f$ x   = |x| \, e^{i\, \phi}\f$, then \n
 *              \f${\rm Im}(x^a) = |x|^a \, \sin(a\, \phi)\f$
 */
#define qcd_CPOWI(x,a) ( (double) pow(qcd_CNORM(x),(a))*sin(qcd_CARG(x)*(a))) 


/******************************************************************************
 * prototypes 
 ******************************************************************************/

/**
 * @fn qcd_complex_16 qcd_CConj(qcd_complex_16 c)
 * @brief Complex conjugate of c.
 * @details \f$ {\tt result} = {\rm Re}(c) - i\,{\rm Im}(c) \f$
 * @param[in] c complex number
 * @return a complex number which is the complex-conjugate of c
 */
qcd_complex_16 qcd_CConj(qcd_complex_16 c);

/**
 * @fn qcd_complex_16 qcd_CMul(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Complex multiplication of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)-{\rm Im}(x){\rm Im}(y) 
 *                           +i[{\rm Re}(x){\rm Im}(y)+{\rm Im}(x){\rm Re}(y)] \f$
 * @param[in] x complex factor 1
 * @param[in] y complex factor 2
 * @return complex product of x and y
 */
qcd_complex_16 qcd_CMul(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CMulR(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Real part of complex multiplication of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)-{\rm Im}(x){\rm Im}(y) \f$
 * @param[in] x complex factor 1
 * @param[in] y complex factor 2
 * @return real part of complex product of x and y
 */
double qcd_CMulR(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CMulI(qcd_complex_16 x,qcd_complex_16 y)
 * @brief Imaginary part of complex multiplication of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x){\rm Im}(y)+{\rm Im}(x){\rm Re}(y) \f$
 * @param[in] x complex factor 1
 * @param[in] y complex factor 2
 * @return imaginary part of complex product of x and y
 */
double qcd_CMulI(qcd_complex_16 x,qcd_complex_16 y);

/**
 * @fn qcd_complex_16 qcd_CConjMul(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Complex multiplication of conj(x) and y.
 * @details \f$ {\tt result} = {\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y) 
 *                           +i[{\rm Re}(x){\rm Im}(y)-{\rm Im}(x){\rm Re}(y)]\f$
 * @param[in] x complex factor 1
 * @param[in] y complex factor 2
 * @return complex product of conj(x) and y
 */
qcd_complex_16 qcd_CConjMul(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn qcd_complex_16 qcd_CAdd(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Complex sum of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x)+{\rm Re}(y) + i[{\rm Im}(x)+{\rm Im}(y)] \f$
 * @param[in] x complex summand 1
 * @param[in] y complex summand 2
 * @return complex sum of x and y
 */
qcd_complex_16 qcd_CAdd(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CAddR(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Real part of complex sum of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x)+{\rm Re}(y) \f$
 * @param[in] x complex summand 1
 * @param[in] y complex summand 2
 * @return sum of real parts of x and y
 */
double qcd_CAddR(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CAddI(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Imaginary part of complex sum of x and y.
 * @details \f$ {\tt result} = {\rm Im}(x)+{\rm Im}(y) \f$
 * @param[in] x complex summand 1
 * @param[in] y complex summand 2
 * @return sum of imaginary parts of x and y
 */
double qcd_CAddI(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn qcd_complex_16 qcd_CSub(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Complex difference of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x)-{\rm Re}(y) + i[{\rm Im}(x)-{\rm Im}(y)] \f$
 * @param[in] x complex minuend
 * @param[in] y complex subtrahend
 * @return complex difference of x and y
 */
qcd_complex_16 qcd_CSub(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CSubR(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Real part of complex difference of x and y.
 * @details \f$ {\tt result} = {\rm Re}(x)-{\rm Re}(y) \f$
 * @param[in] x complex minuend
 * @param[in] y complex subtrahend
 * @return difference of real parts of x and y
 */
double qcd_CSubR(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CSubI(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Imaginary part of complex difference of x and y.
 * @details \f$ {\tt result} = {\rm Im}(x)-{\rm Im}(y) \f$
 * @param[in] x complex minuend
 * @param[in] y complex subtrahend
 * @return difference of imaginary parts of x and y
 */
double qcd_CSubI(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn qcd_complex_16 qcd_CScale(qcd_complex_16 x, double a)
 * @brief product of complex x with real a.
 * @details \f$ {\tt result} = a\,{\rm Re}(x) + i\,a\,{\rm Im}(x) \f$
 * @param[in] x complex factor 1
 * @param[in] a real factor 2
 * @return complex product of x and a
 */
qcd_complex_16 qcd_CScale(qcd_complex_16 x, double a);

/**
 * @fn double qcd_CScaleR(qcd_complex_16 x, double a)
 * @brief Real part of the product of complex x with real a.
 * @details \f$ {\tt result} = a\,{\rm Re}(x) \f$
 * @param[in] x complex factor 1
 * @param[in] a real factor 2
 * @return product of Re(x) and a
 */
double qcd_CScaleR(qcd_complex_16 x, double a);

/**
 * @fn double qcd_CscaleI(qcd_complex_16 x, double a)
 * @brief Imaginary part of the product of complex x with real a.
 * @details \f$ {\tt result} = a\,{\rm Im}(x) \f$
 * @param[in] x complex factor 1
 * @param[in] a real factor 2
 * @return product of Im(x) and a
 */
double qcd_CScaleI(qcd_complex_16 x, double a);

/**
 * @fn qcd_complex_16 qcd_CDiv(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Complex division of x by y.
 * @details \f$ {\tt result} = \frac{  {\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y) 
 *                                   +i[{\rm Im}(x){\rm Re}(y)-{\rm Re}(x){\rm Im}(y)]}
 *                                   {{\rm Re}(y)^2 + {\rm Im}(y)^2)} \f$
 * @param[in] x complex dividend
 * @param[in] y complex divisor
 * @return complex division of x by y
 */
qcd_complex_16 qcd_CDiv(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CDivR(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Real part of complex division of x by y.
 * @details \f$ {\tt result} = \frac{  {\rm Re}(x){\rm Re}(y)+{\rm Im}(x){\rm Im}(y)}
 *                                   {{\rm Re}(y)^2 + {\rm Im}(y)^2)} \f$
 * @param[in] x complex dividend
 * @param[in] y complex divisor
 * @return real part of the division of x by y
 */
double qcd_CDivR(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CDivI(qcd_complex_16 x, qcd_complex_16 y)
 * @brief Imaginary part of complex division of x by y.
 * @details \f$ {\tt result} = \frac{{\rm Im}(x){\rm Re}(y)-{\rm Re}(x){\rm Im}(y)}
 *                                   {{\rm Re}(y)^2 + {\rm Im}(y)^2)} \f$
 * @param[in] x complex dividend
 * @param[in] y complex divisor
 * @return imaginary part of the division of x by y
 */
double qcd_CDivI(qcd_complex_16 x, qcd_complex_16 y);

/**
 * @fn double qcd_CNorm(qcd_complex_16 x)
 * @brief Norm of the complex number x
 * @details \f$ {\tt result} = \sqrt{ {\rm Re}(x)^2 + {\rm Im}(x)^2 } \f$
 * @param[in] x complex number
 * @return norm of x
 */
double qcd_CNorm(qcd_complex_16 x);

/**
 * @fn double qcd_CNormSquared(qcd_complex_16 x)
 * @brief Squared norm of the complex number x
 * @details \f$ {\tt result} =  {\rm Re}(x)^2 + {\rm Im}(x)^2  \f$
 * @param[in] x complex number
 * @return squared norm of x
 */
double qcd_CNormSquared(qcd_complex_16 x);

/**
 * @fn double qcd_CArg(qcd_complex_16 x)
 * @brief The angle &phi; in the polar representation of the complex number x.
 * @details If \f$ x = |x| \, e^{i\,\phi}\f$, then \n
 *           \f$ \phi = {\rm atan2}({\rm Im}(x),{\rm Re}(x)) \f$.\n
 *           (see <a href=https://en.wikipedia.org/wiki/Atan2>atan2</a>)
 * @param[in] x complex number
 * @return argument of x
 */
double qcd_CArg(qcd_complex_16 x);

/**
 * @fn qcd_complex_16 qcd_CPow(qcd_complex_16 x,double a)
 * @brief x to the power of a with complex x and real a
 * @details If \f$ x   = |x|   \, e^{i\, \phi}\f$, then \n
 *             \f$ x^a = |x|^a \, e^{i\,a\,\phi}\f$
 * @param[in] x complex base
 * @param[in] a real exponent
 * @return x to the power a
 */
qcd_complex_16 qcd_CPow(qcd_complex_16 x,double a);

/**
 * @fn double qcd_CPowR(qcd_complex_16 x, double a)
 * @brief Real part of x to the power of a with complex x and real a
 * @details If \f$ x   = |x|   \, e^{i\, \phi}\f$, then \n
 *             \f$ {\rm Re}(x^a) = |x|^a \, \cos{a\,\phi}\f$
 * @param[in] x complex base
 * @param[in] a real exponent
 * @return real part of x to the power a
 */
double qcd_CPowR(qcd_complex_16 x, double a);

/**
 * @fn double qcd_CPowI(qcd_complex_16 x, double a)
 * @brief Imaginary part of x to the power of a with complex x and real a
 * @details If \f$ x   = |x|   \, e^{i\, \phi}\f$, then \n
 *             \f$ {\rm Im}(x^a) = |x|^a \, \sin{a\,\phi}\f$
 * @param[in] x complex base
 * @param[in] a real exponent
 * @return imaginary part of x to the power a
 */
double qcd_CPowI(qcd_complex_16 x, double a);


#endif
