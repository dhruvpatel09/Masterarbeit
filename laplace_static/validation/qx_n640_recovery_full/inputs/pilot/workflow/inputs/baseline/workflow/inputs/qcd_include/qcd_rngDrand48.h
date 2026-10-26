/** @file qcd_rngDrand48.h
 *  @brief header file for qcd_rngDrand48.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2016
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_RNGDRAND48
#define H_QCD_RNGDRAND48 1

/**
 * @fn int qcd_initDrand48(int seed)
 * @brief Initializes the drand48 random number generator.
 * @details drand48 is a linear congruential pseudo random number generator
 *          that creates the sequence \n
 *          \f$ x_{n+1} = (a\, x_n + c)\ \mod n \f$ \n
 *          with
 *          * \f$ a = \mathtt{0x5deece66d} \f$
 *          * \f$ c = \mathtt{0xb} \f$
 *          * \f$ n = 2^{48} \f$
 *          This function initializes \f$ x_0\f$. The low 16 bits are 
 *          set to 0x330E+process-ID, the high 32 bits are set to the value 
 *          seed + process-ID*137137.
 * @param[in] seed determines the initial condition of drand48
 *          
 */
int qcd_initDrand48(int seed);

/**
 * @fn void qcd_randDrand48(double *v, int N)
 * @brief generates double precision random numbers, flatly distributed between 0 and 1
 * @details The the N next 48-bit integers of the drand48 sequence are generated.
 *          their values, divided by 2^48, are stored in v.
 * @param[out] v pointer to allocated memory with space for N*sizeof(double)
 * @param[in] N number of random numbers to generate
 * @warning drand48 is a low quality RNG
 */
void qcd_randDrand48(double *v, int N);

/**
 * @fn void qcd_getRngTypeDrand48(int *rngType, int *n)
 * @brief Returns type and size of the drand48 RNG
 * @param[out] rngType #qcd_RNG_DRAND48
 * @param[out] n 48bits = 6 byte
 */
void qcd_getRngTypeDrand48(int *rngType, int *n);

/**
 * @fn void qcd_getRngStateArrayDrand48(uint8_t *v)
 * @brief Returns a byte-array with the state of the RNG on the current process
 * @details The state of the drand48 random number generator is given by the 
 *          low 48 bits of the 64-bit state integer. These are written in 
 *          little-endian order into the byte-array v.
 * @param v pointer to allocated space for at least 6 bytes
 */
void qcd_getRngStateArrayDrand48(uint8_t *v);

/**
 * @fn void qcd_setRngStateArrayDrand48(uint8_t *v)
 * @brief initializes the drand48 RNG to the state encoded by a given byte-array.
 * @details The state of the drand48 random number generator is given by the 
 *          low 48 bits of the 64-bit state integer. These are obtained from the
 *          6 bytes in the array v, which are interpreted as bytes of the 
 *          state vector in little-endian order.
 * @param v array with 6 bytes
 */
void qcd_setRngStateArrayDrand48(uint8_t *v);

#endif