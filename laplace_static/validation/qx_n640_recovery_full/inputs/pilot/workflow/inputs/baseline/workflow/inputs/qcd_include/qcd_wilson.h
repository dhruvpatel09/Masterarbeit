/** @file qcd_wilson.h
 *  @brief header file for qcd_wilson.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_WILSON
#define H_QCD_WILSON 1

/* some macros 
 * these are valid in the TM chiral representation
 * and in the alpha-collab chiral representation
 * assuming that x and y are 4*3*2 component vectors
 * with index i= (is_imaginary)+2*color+2*3*dirac 
 */

/**
 * @def qcd_ONE_PLUS_GAMMA0(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_0\f$.
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_0\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA0(x,y) \
  {  x[0] = y[0] - y[12]; \
     x[1] = y[1] - y[13]; \
     x[2] = y[2] - y[14]; \
     x[3] = y[3] - y[15]; \
     x[4] = y[4] - y[16]; \
     x[5] = y[5] - y[17]; \
                          \
     x[6] = y[6] - y[18]; \
     x[7] = y[7] - y[19]; \
     x[8] = y[8] - y[20]; \
     x[9] = y[9] - y[21]; \
     x[10]= y[10]- y[22]; \
     x[11]= y[11]- y[23]; \
                          \
     x[12]= -x[0]; \
     x[13]= -x[1]; \
     x[14]= -x[2]; \
     x[15]= -x[3]; \
     x[16]= -x[4]; \
     x[17]= -x[5]; \
                   \
     x[18]= -x[6]; \
     x[19]= -x[7]; \
     x[20]= -x[8]; \
     x[21]= -x[9]; \
     x[22]= -x[10];\
     x[23]= -x[11];\
   }

/**
 * @def qcd_ONE_PLUS_GAMMA0_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_0\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_0\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA0_HALFSPINOR(x,y) \
  {  x[0] = y[0] - y[12]; \
     x[1] = y[1] - y[13]; \
     x[2] = y[2] - y[14]; \
     x[3] = y[3] - y[15]; \
     x[4] = y[4] - y[16]; \
     x[5] = y[5] - y[17]; \
                          \
     x[6] = y[6] - y[18]; \
     x[7] = y[7] - y[19]; \
     x[8] = y[8] - y[20]; \
     x[9] = y[9] - y[21]; \
     x[10]= y[10]- y[22]; \
     x[11]= y[11]- y[23]; \
   }

/**
 * @def qcd_ONE_PLUS_GAMMA1(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_1\f$
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_1\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA1(x,y) \
  {  x[0] = y[0] + y[19]; \
     x[1] = y[1] - y[18]; \
     x[2] = y[2] + y[21]; \
     x[3] = y[3] - y[20]; \
     x[4] = y[4] + y[23]; \
     x[5] = y[5] - y[22]; \
                          \
     x[6] = y[6] + y[13]; \
     x[7] = y[7] - y[12]; \
     x[8] = y[8] + y[15]; \
     x[9] = y[9] - y[14]; \
     x[10]= y[10]+ y[17]; \
     x[11]= y[11]- y[16]; \
                          \
     x[12]= -x[7]; \
     x[13]=  x[6]; \
     x[14]= -x[9]; \
     x[15]=  x[8]; \
     x[16]= -x[11];\
     x[17]=  x[10];\
                   \
     x[18]= -x[1]; \
     x[19]=  x[0]; \
     x[20]= -x[3]; \
     x[21]=  x[2]; \
     x[22]= -x[5]; \
     x[23]=  x[4]; \
   }

/**
 * @def qcd_ONE_PLUS_GAMMA1_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_1\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_1\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA1_HALFSPINOR(x,y) \
  {  x[0] = y[0] + y[19]; \
     x[1] = y[1] - y[18]; \
     x[2] = y[2] + y[21]; \
     x[3] = y[3] - y[20]; \
     x[4] = y[4] + y[23]; \
     x[5] = y[5] - y[22]; \
                          \
     x[6] = y[6] + y[13]; \
     x[7] = y[7] - y[12]; \
     x[8] = y[8] + y[15]; \
     x[9] = y[9] - y[14]; \
     x[10]= y[10]+ y[17]; \
     x[11]= y[11]- y[16]; \
   }

/**
 * @def qcd_ONE_PLUS_GAMMA2(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_2\f$
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_2\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA2(x,y) \
  {  x[0] = y[0] - y[18]; \
     x[1] = y[1] - y[19]; \
     x[2] = y[2] - y[20]; \
     x[3] = y[3] - y[21]; \
     x[4] = y[4] - y[22]; \
     x[5] = y[5] - y[23]; \
                          \
     x[6] = y[6] + y[12]; \
     x[7] = y[7] + y[13]; \
     x[8] = y[8] + y[14]; \
     x[9] = y[9] + y[15]; \
     x[10]= y[10]+ y[16]; \
     x[11]= y[11]+ y[17]; \
                          \
     x[12]=   x[6]; \
     x[13]=   x[7]; \
     x[14]=   x[8]; \
     x[15]=   x[9]; \
     x[16]=   x[10];\
     x[17]=   x[11];\
                    \
     x[18]=  -x[0]; \
     x[19]=  -x[1]; \
     x[20]=  -x[2]; \
     x[21]=  -x[3]; \
     x[22]=  -x[4]; \
     x[23]=  -x[5]; \
   }


/**
 * @def qcd_ONE_PLUS_GAMMA2_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_2\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_2\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA2_HALFSPINOR(x,y) \
  {  x[0] = y[0] - y[18]; \
     x[1] = y[1] - y[19]; \
     x[2] = y[2] - y[20]; \
     x[3] = y[3] - y[21]; \
     x[4] = y[4] - y[22]; \
     x[5] = y[5] - y[23]; \
                          \
     x[6] = y[6] + y[12]; \
     x[7] = y[7] + y[13]; \
     x[8] = y[8] + y[14]; \
     x[9] = y[9] + y[15]; \
     x[10]= y[10]+ y[16]; \
     x[11]= y[11]+ y[17]; \
   }


/**
 * @def qcd_ONE_PLUS_GAMMA3(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_3\f$
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_3\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA3(x,y) \
  {  x[0] = y[0] + y[13]; \
     x[1] = y[1] - y[12]; \
     x[2] = y[2] + y[15]; \
     x[3] = y[3] - y[14]; \
     x[4] = y[4] + y[17]; \
     x[5] = y[5] - y[16]; \
                          \
     x[6] = y[6] - y[19]; \
     x[7] = y[7] + y[18]; \
     x[8] = y[8] - y[21]; \
     x[9] = y[9] + y[20]; \
     x[10]= y[10]- y[23]; \
     x[11]= y[11]+ y[22]; \
                          \
     x[12]=  -x[1]; \
     x[13]=   x[0]; \
     x[14]=  -x[3]; \
     x[15]=   x[2]; \
     x[16]=  -x[5]; \
     x[17]=   x[4]; \
                    \
     x[18]=   x[7]; \
     x[19]=  -x[6]; \
     x[20]=   x[9]; \
     x[21]=  -x[8]; \
     x[22]=   x[11];\
     x[23]=  -x[10];\
   }

/**
 * @def qcd_ONE_PLUS_GAMMA3_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}+\gamma_3\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}+\gamma_3\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_PLUS_GAMMA3_HALFSPINOR(x,y) \
  {  x[0] = y[0] + y[13]; \
     x[1] = y[1] - y[12]; \
     x[2] = y[2] + y[15]; \
     x[3] = y[3] - y[14]; \
     x[4] = y[4] + y[17]; \
     x[5] = y[5] - y[16]; \
                          \
     x[6] = y[6] - y[19]; \
     x[7] = y[7] + y[18]; \
     x[8] = y[8] - y[21]; \
     x[9] = y[9] + y[20]; \
     x[10]= y[10]- y[23]; \
     x[11]= y[11]+ y[22]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA0(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_0\f$
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_0\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA0(x,y) \
  {  x[0] = y[0] + y[12]; \
     x[1] = y[1] + y[13]; \
     x[2] = y[2] + y[14]; \
     x[3] = y[3] + y[15]; \
     x[4] = y[4] + y[16]; \
     x[5] = y[5] + y[17]; \
                          \
     x[6] = y[6] + y[18]; \
     x[7] = y[7] + y[19]; \
     x[8] = y[8] + y[20]; \
     x[9] = y[9] + y[21]; \
     x[10]= y[10]+ y[22]; \
     x[11]= y[11]+ y[23]; \
                          \
     x[12]=  x[0]; \
     x[13]=  x[1]; \
     x[14]=  x[2]; \
     x[15]=  x[3]; \
     x[16]=  x[4]; \
     x[17]=  x[5]; \
                   \
     x[18]=  x[6]; \
     x[19]=  x[7]; \
     x[20]=  x[8]; \
     x[21]=  x[9]; \
     x[22]=  x[10];\
     x[23]=  x[11];\
   }

/**
 * @def qcd_ONE_MINUS_GAMMA0_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_0\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_0\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA0_HALFSPINOR(x,y) \
  {  x[0] = y[0] + y[12]; \
     x[1] = y[1] + y[13]; \
     x[2] = y[2] + y[14]; \
     x[3] = y[3] + y[15]; \
     x[4] = y[4] + y[16]; \
     x[5] = y[5] + y[17]; \
                          \
     x[6] = y[6] + y[18]; \
     x[7] = y[7] + y[19]; \
     x[8] = y[8] + y[20]; \
     x[9] = y[9] + y[21]; \
     x[10]= y[10]+ y[22]; \
     x[11]= y[11]+ y[23]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA1(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_1\f$
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_1\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA1(x,y) \
  {  x[0] = y[0] - y[19]; \
     x[1] = y[1] + y[18]; \
     x[2] = y[2] - y[21]; \
     x[3] = y[3] + y[20]; \
     x[4] = y[4] - y[23]; \
     x[5] = y[5] + y[22]; \
                          \
     x[6] = y[6] - y[13]; \
     x[7] = y[7] + y[12]; \
     x[8] = y[8] - y[15]; \
     x[9] = y[9] + y[14]; \
     x[10]= y[10]- y[17]; \
     x[11]= y[11]+ y[16]; \
                          \
     x[12]=  x[7]; \
     x[13]= -x[6]; \
     x[14]=  x[9]; \
     x[15]= -x[8]; \
     x[16]=  x[11];\
     x[17]= -x[10];\
                   \
     x[18]=  x[1]; \
     x[19]= -x[0]; \
     x[20]=  x[3]; \
     x[21]= -x[2]; \
     x[22]=  x[5]; \
     x[23]= -x[4]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA1_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_1\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_1\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA1_HALFSPINOR(x,y) \
  {  x[0] = y[0] - y[19]; \
     x[1] = y[1] + y[18]; \
     x[2] = y[2] - y[21]; \
     x[3] = y[3] + y[20]; \
     x[4] = y[4] - y[23]; \
     x[5] = y[5] + y[22]; \
                          \
     x[6] = y[6] - y[13]; \
     x[7] = y[7] + y[12]; \
     x[8] = y[8] - y[15]; \
     x[9] = y[9] + y[14]; \
     x[10]= y[10]- y[17]; \
     x[11]= y[11]+ y[16]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA2(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_2\f$.
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_2\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA2(x,y) \
  {  x[0] = y[0] + y[18]; \
     x[1] = y[1] + y[19]; \
     x[2] = y[2] + y[20]; \
     x[3] = y[3] + y[21]; \
     x[4] = y[4] + y[22]; \
     x[5] = y[5] + y[23]; \
                          \
     x[6] = y[6] - y[12]; \
     x[7] = y[7] - y[13]; \
     x[8] = y[8] - y[14]; \
     x[9] = y[9] - y[15]; \
     x[10]= y[10]- y[16]; \
     x[11]= y[11]- y[17]; \
                          \
     x[12]=  -x[6]; \
     x[13]=  -x[7]; \
     x[14]=  -x[8]; \
     x[15]=  -x[9]; \
     x[16]=  -x[10];\
     x[17]=  -x[11];\
                   \
     x[18]=  x[0]; \
     x[19]=  x[1]; \
     x[20]=  x[2]; \
     x[21]=  x[3]; \
     x[22]=  x[4]; \
     x[23]=  x[5]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA2_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_2\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_2\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA2_HALFSPINOR(x,y) \
  {  x[0] = y[0] + y[18]; \
     x[1] = y[1] + y[19]; \
     x[2] = y[2] + y[20]; \
     x[3] = y[3] + y[21]; \
     x[4] = y[4] + y[22]; \
     x[5] = y[5] + y[23]; \
                          \
     x[6] = y[6] - y[12]; \
     x[7] = y[7] - y[13]; \
     x[8] = y[8] - y[14]; \
     x[9] = y[9] - y[15]; \
     x[10]= y[10]- y[16]; \
     x[11]= y[11]- y[17]; \
   }

/**
 * @def qcd_ONE_MINUS_GAMMA3(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_3\f$.
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_3\f$ is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA3(x,y) \
  {  x[0] = y[0] - y[13]; \
     x[1] = y[1] + y[12]; \
     x[2] = y[2] - y[15]; \
     x[3] = y[3] + y[14]; \
     x[4] = y[4] - y[17]; \
     x[5] = y[5] + y[16]; \
                          \
     x[6] = y[6] + y[19]; \
     x[7] = y[7] - y[18]; \
     x[8] = y[8] + y[21]; \
     x[9] = y[9] - y[20]; \
     x[10]= y[10]+ y[23]; \
     x[11]= y[11]- y[22]; \
                          \
     x[12]=   x[1]; \
     x[13]=  -x[0]; \
     x[14]=   x[3]; \
     x[15]=  -x[2]; \
     x[16]=   x[5]; \
     x[17]=  -x[4]; \
                    \
     x[18]=  -x[7]; \
     x[19]=   x[6]; \
     x[20]=  -x[9]; \
     x[21]=   x[8]; \
     x[22]=  -x[11];\
     x[23]=   x[10];\
   }

/**
 * @def qcd_ONE_MINUS_GAMMA3_HALFSPINOR(x,y)
 * @brief Multiply a spinor by \f${\mathds 1}-\gamma_3\f$, compute only half the components
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f${\mathds 1}-\gamma_3\f$ is applied to y. Only the first 12
 *          components of x are computed. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra.
 */
#define qcd_ONE_MINUS_GAMMA3_HALFSPINOR(x,y) \
  {  x[0] = y[0] - y[13]; \
     x[1] = y[1] + y[12]; \
     x[2] = y[2] - y[15]; \
     x[3] = y[3] + y[14]; \
     x[4] = y[4] - y[17]; \
     x[5] = y[5] + y[16]; \
                          \
     x[6] = y[6] + y[19]; \
     x[7] = y[7] - y[18]; \
     x[8] = y[8] + y[21]; \
     x[9] = y[9] - y[20]; \
     x[10]= y[10]+ y[23]; \
     x[11]= y[11]- y[22]; \
   }


/**
 * @def qcd_MASSTERMS(x,y,mp4,mu)
 * @brief Multiply a spinor by \f$({\tt mp4}{\mathds 1}+i*\gamma_5 {\tt mu})\f$.
 * @details It is assumed that the 24 components of y contain
 *          a spinor, such that the real part of a color/dirac component
 *          is at position 2*color+6*dirac and the imaginary at
 *          2*color+6*dirac+1. The x is a spinor in the same ordering
 *          that is obtained, if \f$({\tt mp4}{\mathds 1}+i*\gamma_5 {\tt mu})\f$
 *          is applied to y. This works
 *          only int the TM chiral representation or the alpha-collaboration
 *          chiral representation of the clifford algebra. In Dirac operators
 *          mp4 is usually m0 + 4 (hence the name).
 */
#define qcd_MASSTERMS(x,y,mp4,mu) \
  {  x[0] = mp4*y[0]-mu*y[1];  \
     x[1] = mp4*y[1]+mu*y[0];  \
     x[2] = mp4*y[2]-mu*y[3];  \
     x[3] = mp4*y[3]+mu*y[2];  \
     x[4] = mp4*y[4]-mu*y[5];  \
     x[5] = mp4*y[5]+mu*y[4];  \
                               \
     x[6] = mp4*y[6] -mu*y[7]; \
     x[7] = mp4*y[7] +mu*y[6]; \
     x[8] = mp4*y[8] -mu*y[9]; \
     x[9] = mp4*y[9] +mu*y[8]; \
     x[10]= mp4*y[10]-mu*y[11];\
     x[11]= mp4*y[11]+mu*y[10];\
                               \
     x[12]= mp4*y[12]+mu*y[13];\
     x[13]= mp4*y[13]-mu*y[12];\
     x[14]= mp4*y[14]+mu*y[15];\
     x[15]= mp4*y[15]-mu*y[14];\
     x[16]= mp4*y[16]+mu*y[17];\
     x[17]= mp4*y[17]-mu*y[16];\
                               \
     x[18]= mp4*y[18]+mu*y[19];\
     x[19]= mp4*y[19]-mu*y[18];\
     x[20]= mp4*y[20]+mu*y[21];\
     x[21]= mp4*y[21]-mu*y[20];\
     x[22]= mp4*y[22]+mu*y[23];\
     x[23]= mp4*y[23]-mu*y[22];\
   }



/**
 * @def qcd_SUM_OP(x,a,b,c,d,e,f,g,h,i)
 * @brief sum x[j] = a[j] -(b[j]+c[j]+d[j]+e[j]+f[j]+g[j]+h[j]+i[j])/2 for j=0..23
 */
#define qcd_SUM_OP(x,a,b,c,d,e,f,g,h,i)  \
   {  x[0] = a[0] -0.5*(b[0] +c[0] +d[0] +e[0] +f[0] +g[0] +h[0] +i[0]);        \
      x[1] = a[1] -0.5*(b[1] +c[1] +d[1] +e[1] +f[1] +g[1] +h[1] +i[1]);        \
      x[2] = a[2] -0.5*(b[2] +c[2] +d[2] +e[2] +f[2] +g[2] +h[2] +i[2]);        \
      x[3] = a[3] -0.5*(b[3] +c[3] +d[3] +e[3] +f[3] +g[3] +h[3] +i[3]);        \
      x[4] = a[4] -0.5*(b[4] +c[4] +d[4] +e[4] +f[4] +g[4] +h[4] +i[4]);        \
      x[5] = a[5] -0.5*(b[5] +c[5] +d[5] +e[5] +f[5] +g[5] +h[5] +i[5]);        \
      x[6] = a[6] -0.5*(b[6] +c[6] +d[6] +e[6] +f[6] +g[6] +h[6] +i[6]);        \
      x[7] = a[7] -0.5*(b[7] +c[7] +d[7] +e[7] +f[7] +g[7] +h[7] +i[7]);        \
      x[8] = a[8] -0.5*(b[8] +c[8] +d[8] +e[8] +f[8] +g[8] +h[8] +i[8]);        \
      x[9] = a[9] -0.5*(b[9] +c[9] +d[9] +e[9] +f[9] +g[9] +h[9] +i[9]);        \
      x[10]= a[10]-0.5*(b[10]+c[10]+d[10]+e[10]+f[10]+g[10]+h[10]+i[10]);       \
      x[11]= a[11]-0.5*(b[11]+c[11]+d[11]+e[11]+f[11]+g[11]+h[11]+i[11]);       \
      x[12]= a[12]-0.5*(b[12]+c[12]+d[12]+e[12]+f[12]+g[12]+h[12]+i[12]);       \
      x[13]= a[13]-0.5*(b[13]+c[13]+d[13]+e[13]+f[13]+g[13]+h[13]+i[13]);       \
      x[14]= a[14]-0.5*(b[14]+c[14]+d[14]+e[14]+f[14]+g[14]+h[14]+i[14]);       \
      x[15]= a[15]-0.5*(b[15]+c[15]+d[15]+e[15]+f[15]+g[15]+h[15]+i[15]);       \
      x[16]= a[16]-0.5*(b[16]+c[16]+d[16]+e[16]+f[16]+g[16]+h[16]+i[16]);       \
      x[17]= a[17]-0.5*(b[17]+c[17]+d[17]+e[17]+f[17]+g[17]+h[17]+i[17]);       \
      x[18]= a[18]-0.5*(b[18]+c[18]+d[18]+e[18]+f[18]+g[18]+h[18]+i[18]);       \
      x[19]= a[19]-0.5*(b[19]+c[19]+d[19]+e[19]+f[19]+g[19]+h[19]+i[19]);       \
      x[20]= a[20]-0.5*(b[20]+c[20]+d[20]+e[20]+f[20]+g[20]+h[20]+i[20]);       \
      x[21]= a[21]-0.5*(b[21]+c[21]+d[21]+e[21]+f[21]+g[21]+h[21]+i[21]);       \
      x[22]= a[22]-0.5*(b[22]+c[22]+d[22]+e[22]+f[22]+g[22]+h[22]+i[22]);       \
      x[23]= a[23]-0.5*(b[23]+c[23]+d[23]+e[23]+f[23]+g[23]+h[23]+i[23]);       \
   }

/**
 * @def qcd_SUM_OP_HALFSPINORS(x,a,b,c,d,e,f,g,h,i)
 * @brief sum x[j] = a[j] -(b[j]+c[j]+d[j]+e[j]+f[j]+g[j]+h[j]+i[j])/2 for j=0..23
 * @details Only components j=0..11 of b..i are used. It is assumed that b..i 
 *          are \f$ [1\pm \gamma_\mu]\psi \f$ so that the other half of the
 *          components can be obtained through symmetries.
 */
#define qcd_SUM_OP_HALFSPINORS(x,a,b,c,d,e,f,g,h,i)  \
   {  x[0] = a[0] -0.5*(b[0] +c[0] +d[0] +e[0] +f[0] +g[0] +h[0] +i[0]);        \
      x[1] = a[1] -0.5*(b[1] +c[1] +d[1] +e[1] +f[1] +g[1] +h[1] +i[1]);        \
      x[2] = a[2] -0.5*(b[2] +c[2] +d[2] +e[2] +f[2] +g[2] +h[2] +i[2]);        \
      x[3] = a[3] -0.5*(b[3] +c[3] +d[3] +e[3] +f[3] +g[3] +h[3] +i[3]);        \
      x[4] = a[4] -0.5*(b[4] +c[4] +d[4] +e[4] +f[4] +g[4] +h[4] +i[4]);        \
      x[5] = a[5] -0.5*(b[5] +c[5] +d[5] +e[5] +f[5] +g[5] +h[5] +i[5]);        \
      x[6] = a[6] -0.5*(b[6] +c[6] +d[6] +e[6] +f[6] +g[6] +h[6] +i[6]);        \
      x[7] = a[7] -0.5*(b[7] +c[7] +d[7] +e[7] +f[7] +g[7] +h[7] +i[7]);        \
      x[8] = a[8] -0.5*(b[8] +c[8] +d[8] +e[8] +f[8] +g[8] +h[8] +i[8]);        \
      x[9] = a[9] -0.5*(b[9] +c[9] +d[9] +e[9] +f[9] +g[9] +h[9] +i[9]);        \
      x[10]= a[10]-0.5*(b[10]+c[10]+d[10]+e[10]+f[10]+g[10]+h[10]+i[10]);       \
      x[11]= a[11]-0.5*(b[11]+c[11]+d[11]+e[11]+f[11]+g[11]+h[11]+i[11]);       \
      x[12]= a[12]-0.5*(-b[0] +c[0] -d[7] +e[7] +f[6] -g[6] -h[1] +i[1] );      \
      x[13]= a[13]-0.5*(-b[1] +c[1] +d[6] -e[6] +f[7] -g[7] +h[0] -i[0] );      \
      x[14]= a[14]-0.5*(-b[2] +c[2] -d[9] +e[9] +f[8] -g[8] -h[3] +i[3] );      \
      x[15]= a[15]-0.5*(-b[3] +c[3] +d[8] -e[8] +f[9] -g[9] +h[2] -i[2] );      \
      x[16]= a[16]-0.5*(-b[4] +c[4] -d[11]+e[11]+f[10]-g[10]-h[5] +i[5] );      \
      x[17]= a[17]-0.5*(-b[5] +c[5] +d[10]-e[10]+f[11]-g[11]+h[4] -i[4] );      \
      x[18]= a[18]-0.5*(-b[6] +c[6] -d[1] +e[1] -f[0] +g[0] +h[7] -i[7] );      \
      x[19]= a[19]-0.5*(-b[7] +c[7] +d[0] -e[0] -f[1] +g[1] -h[6] +i[6] );      \
      x[20]= a[20]-0.5*(-b[8] +c[8] -d[3] +e[3] -f[2] +g[2] +h[9] -i[9] );      \
      x[21]= a[21]-0.5*(-b[9] +c[9] +d[2] -e[2] -f[3] +g[3] -h[8] +i[8] );      \
      x[22]= a[22]-0.5*(-b[10]+c[10]-d[5] +e[5] -f[4] +g[4] +h[11]-i[11]);      \
      x[23]= a[23]-0.5*(-b[11]+c[11]+d[4] -e[4] -f[5] +g[5] -h[10]+i[10]);      \
   }

/**
 * @def qcd_PHASE_MUL(x,t,pr,pi)
 * @brief Multiplies a spinor by a complex phase.
 * @details A spinor is multiplied by the phase
 *          \f$ e^{i\phi} = {\tt pr} +i {\tt pi} \f$.
 *          The spinor components are ordered such, that even indices are 
 *          real parts, and odd indices are imaginary parts.
 */
#define qcd_PHASE_MUL(x,t,pr,pi)   \
   {  t    =x[0] *pr -x[1] *pi;    \
      x[1] =x[0] *pi +x[1] *pr;    \
      x[0] =t;                     \
      t    =x[2] *pr -x[3] *pi;    \
      x[3] =x[2] *pi +x[3] *pr;    \
      x[2] =t;                     \
      t    =x[4] *pr -x[5] *pi;    \
      x[5] =x[4] *pi +x[5] *pr;    \
      x[4] =t;                     \
      t    =x[6] *pr -x[7] *pi;    \
      x[7] =x[6] *pi +x[7] *pr;    \
      x[6] =t;                     \
      t    =x[8] *pr -x[9] *pi;    \
      x[9] =x[8] *pi +x[9] *pr;    \
      x[8] =t;                     \
      t    =x[10]*pr -x[11]*pi;    \
      x[11]=x[10]*pi +x[11]*pr;    \
      x[10]=t;                     \
      t    =x[12]*pr -x[13]*pi;    \
      x[13]=x[12]*pi +x[13]*pr;    \
      x[12]=t;                     \
      t    =x[14]*pr -x[15]*pi;    \
      x[15]=x[14]*pi +x[15]*pr;    \
      x[14]=t;                     \
      t    =x[16]*pr -x[17]*pi;    \
      x[17]=x[16]*pi +x[17]*pr;    \
      x[16]=t;                     \
      t    =x[18]*pr -x[19]*pi;    \
      x[19]=x[18]*pi +x[19]*pr;    \
      x[18]=t;                     \
      t    =x[20]*pr -x[21]*pi;    \
      x[21]=x[20]*pi +x[21]*pr;    \
      x[20]=t;                     \
      t    =x[22]*pr -x[23]*pi;    \
      x[23]=x[22]*pi +x[23]*pr;    \
      x[22]=t;                     \
   }

/**
 * @def qcd_PHASE_MUL_HALFSPINOR(x,t,pr,pi)
 * @brief Multiplies a half-spinor by a complex phase.
 * @details A spinor is multiplied by the phase
 *          \f$ e^{i\phi} = {\tt pr} +i {\tt pi} \f$.
 *          The spinor components are ordered such, that even indices are 
 *          real parts, and odd indices are imaginary parts. The operation
 *          is performed on the first half of the spinor's components.
 */
#define qcd_PHASE_MUL_HALFSPINOR(x,t,pr,pi)   \
   {  t    =x[0] *pr -x[1] *pi;    \
      x[1] =x[0] *pi +x[1] *pr;    \
      x[0] =t;                     \
      t    =x[2] *pr -x[3] *pi;    \
      x[3] =x[2] *pi +x[3] *pr;    \
      x[2] =t;                     \
      t    =x[4] *pr -x[5] *pi;    \
      x[5] =x[4] *pi +x[5] *pr;    \
      x[4] =t;                     \
      t    =x[6] *pr -x[7] *pi;    \
      x[7] =x[6] *pi +x[7] *pr;    \
      x[6] =t;                     \
      t    =x[8] *pr -x[9] *pi;    \
      x[9] =x[8] *pi +x[9] *pr;    \
      x[8] =t;                     \
      t    =x[10]*pr -x[11]*pi;    \
      x[11]=x[10]*pi +x[11]*pr;    \
      x[10]=t;                     \
   }


 /* prototypes */
/**
 * @fn void qcd_onePlusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], int mu)
 * @brief Applies \f$ 1+\gamma_\mu \f$ to a spinor.
 * @param[out] result pointer to a spinor with 4x3 complex components
 * @param[in] vec pointer to a spinor with 4x3 complex components
 * @param[in] mu specifies which projector to apply mu=0,1,2,3
 */
 void qcd_onePlusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], int mu);
 
 /**
 * @fn void qcd_oneMinusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], int mu)
 * @brief Applies \f$ 1-\gamma_\mu \f$ to a spinor.
 * @param[out] result pointer to a spinor with 4x3 complex components
 * @param[in] vec pointer to a spinor with 4x3 complex components
 * @param[in] mu specifies which projector to apply mu=0,1,2,3
 */
 void qcd_oneMinusGamma(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], int mu);
 
 /**
  * @fn void qcd_applyU(qcd_complex_16 result[][3],qcd_complex_16 u[][3],qcd_complex_16 vec[][3])
  * @brief Multiplies a spinor with a gauge field.
  * @details The spinor vec is multiplied with the 3x3 matrix u. The result is
  *          stored in result.
  * @param[out] result pointer to a 4x3 complex spinor
  * @param[in] u pointer to a 3x3 complex matrix
  * @param[in] vec pointer to a 4x3 complex spinor
  */
 void qcd_applyU(qcd_complex_16 result[][3],qcd_complex_16 u[][3],qcd_complex_16 vec[][3]);

/**
  * @fn void qcd_applyAdjointU(qcd_complex_16 result[][3],qcd_complex_16 u[][3],qcd_complex_16 vec[][3])
  * @brief Multiplies a spinor with a adjoint gauge field.
  * @details The spinor vec is multiplied with the adjoint of the 3x3 matrix u. 
  *          The result is stored in result.
  * @param[out] result pointer to a 4x3 complex spinor
  * @param[in] u pointer to a 3x3 complex matrix
  * @param[in] vec pointer to a 4x3 complex spinor
  */
 void qcd_applyAdjointU(qcd_complex_16 result[][3], qcd_complex_16 u[][3],qcd_complex_16 vec[][3]);
 
 /**
  * @fn void qcd_massTerms(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], double mp4, double tm)
  * @brief Applies \f$({\mathds 1}mp4 + i\gamma_5{\tt tm}) \f$ to a spinor.
  * @details The spinor vec is multiplied by the 4x4 matrix 
  *          \f$({\mathds 1}mp4 + i\gamma_5{\tt tm}) \f$ and the result is 
  *          stored in the spinor result.
  * @param[out] result pointer to a 4x3 complex spinor
  * @param[in] vec pointer to a 4x3 complex spinor
  * @param[in] mp4 bare mass plus 4
  * @param[in] tm bare twisted mass
  */
 void qcd_massTerms(qcd_complex_16 result[][3], qcd_complex_16 vec[][3], double mp4, double tm);
 
 /**
  * @fn void qcd_phaseMul(qcd_complex_16 vec[][3], qcd_complex_16 phi)
  * @brief Multiplies a spinor by a complex number.
  * @details the 4x3 component spinor vec is multiplied (in place) by the
  *          complex number phi.
  * @param[in,out] vec pointer to a 4x3 complex spinor
  * @param[in] phi a complex number
  */
 void qcd_phaseMul(qcd_complex_16 vec[][3], qcd_complex_16 phi);
 
 /**
  * @fn void qcd_phaseMul3(qcd_complex_16 vec[3], qcd_complex_16 phi)
  * @brief Multiplies a spinor component by a complex number.
  * @details the 3-component color vector  vec is multiplied (in place) by the
  *          complex number phi.
  * @param[in,out] vec pointer to a 3 component complex spinor
  * @param[in] phi a complex number
  */
 void qcd_phaseMul3(qcd_complex_16 ve[3], qcd_complex_16 phi);



  /**
  * @fn int qcd_DVector(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, int mu)
  * @brief Applies the covariant forward derivative to a spinor field
  * @details Dv(x) = U_mu(x) v(x+mu) - v(x)
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] mu direction of the derivative (0...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DVector(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, int mu);
 
 /**
  * @fn int qcd_DBackVector(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, int mu)
  * @brief Applies the backward covariant forward derivative to a spinor field
  * @details Dv(x) = v(x) - U^dag_mu(x-mu) v(x-mu)
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] mu direction of the derivative (0...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DBackVector(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, int mu);
 
 /**
  * @fn int qcd_applyWilsonTMOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm)
  * @brief Applies the Wilson twisted mass operator on a spinor field.
  * @details The (uni-improved) Wilson operator with bare mass m and twisted mass
  *           tm, with a given gauge field u, is applied to the spinor field v.
  *           The result is stored in Dv.
  *           All processes must call this function simultaneously.
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] m bare mass parameter
  * @param[in] tm bare twisted mass parameter
  * @return 0 on success, 1 on failure
  * @todo Make gauge field communication optional.
  * @todo communicate only half-spinors
  * @todo write a fast version with vector instructions and/or openMP
  */
 int qcd_applyWilsonTMOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm);
 
 /**
  * @fn int qcd_applyWilsonTMOpNoScripts(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm)
  * @brief Applies the Wilson twisted mass operator on a spinor field.
  * @details The (uni-improved) Wilson operator with bare mass m and twisted mass
  *           tm, with a given gauge field u, is applied to the spinor field v.
  *           The result is stored in Dv.
  *           All processes must call this function simultaneously.
  *           This version does not rely on preprocessor scripts for the
  *           \f$ {\mathds 1} \pm \gamma_\mu \f$ projectors.
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] m bare mass parameter
  * @param[in] tm bare twisted mass parameter
  * @return 0 on success, 1 on failure
  */
 int qcd_applyWilsonTMOpNoScripts(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm);
 
 /**
  * @fn int qcd_applyWilsonTMOpNoOptim(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm)
  * @brief Applies the Wilson twisted mass operator on a spinor field.
  * @details The (uni-improved) Wilson operator with bare mass m and twisted mass
  *           tm, with a given gauge field u, is applied to the spinor field v.
  *           The result is stored in Dv.
  *           All processes must call this function simultaneously.
  *           This version does not use half-spinors.
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] m bare mass parameter
  * @param[in] tm bare twisted mass parameter
  * @return 0 on success, 1 on failure
  */
 int qcd_applyWilsonTMOpNoOptim(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, double m, double tm);
 
 
 /**
  * @fn int qcd_makeCloverDiagonal(qcd_cloverDiagonal *diag, qcd_cloverField *q, double m0, double tm, double csw)
  * @brief Create the spacetime-diagonal part of the O(a) improved Wilson operator.
  * @details Given the bare mass m0, the bare twisted mass tm, the improvement parameter
  *          csw and the clover field \f$ Q_{\mu \nu}\f$ given by q (calculated
  *          using #qcd_calculateCloverField), the space-time
  *          diagonal part of the improved Wilson operator is computed and stored
  *          in diag.
  * @param[in,out] diag pointer to an initialized qcd_cloverDiagonal structure
  * @param[in] q pointer to an already computed clover field in a qcd_cloverField structure
  * @param[in] m0 bare mass parameter
  * @param[in] tm bare twisted mass parameter
  * @param[in] csw improvement coefficient
  * @return 0 on success, 1 on failure
  * @bug This is not working correctly.
  */
 int qcd_makeCloverDiagonal(qcd_cloverDiagonal *diag, qcd_cloverField *q, double m0, double tm, double csw);
 
 /**
  * @fn int qcd_applyWilsonTMCloverOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_cloverDiagonal *diagD)
  * @brief Apply the O(a) improved Wilson operator with twisted mass term.
  * @details Given the gauge field u  and the clover-diagonal field diagD, 
  *          the clover Wilson
  *          operator is applied to spinor field v. The result is stored in 
  *          Dv. Bare mass, bare twisted mass and the improvement coefficient
  *          are all inside diagD, which should be created with qcd_CloverDiagonal
  *          whenever u has changed.
  *          All processes must call this function simultaneously.
  * @param[in,out] Dv pointer to an initialized qcd_vector structure
  * @param[in] v pointer to an initialized qcd_vector structure
  * @param[in] u pointer to an initialized qcd_gaugeField structure
  * @param[in] diagD pointer to an already computed qcd_cloverDiagonal structure
  * @return 0 on success, 1 on failure
  * @bug This is not working correctly yet (unless csw=0).
  */
 int qcd_applyWilsonTMCloverOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_cloverDiagonal *diagD);
 
 /**
  * @fn int qcd_isSFVector(qcd_vector *v)
  * @brief Tests whether a spinor field satisfies SF boundary conditions
  * @details Given the spinor field v, it is tested whetehr its values at 
  *          global coordinates x0=0 and x0=T-1 are all zero. 
  *          All processes must call this function simultaneously.
  * @param v pointer to an initialized qcd_vector structure
  * @return 1 if v satisfies SF b.c., 0 otherwise
  */
 int qcd_isSFVector(qcd_vector *v);
 
 /**
  * @fn void qcd_forceSFVector(qcd_vector *v)
  * @brief Sets the boundary values of v to zero.
  * @details Given the spinor field v, its values at 
  *          global coordinates x0=0 and x0=T-1 are all set to zero. 
  *          All processes must call this function simultaneously.
  * @param v pointer to an initialized qcd_vector structure
  */
 void qcd_forceSFVector(qcd_vector *v);
#endif
