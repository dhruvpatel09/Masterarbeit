
/*******************************************************************************
*
* file qcd_ranlux.h
*
* Copyright (C) 2005 Martin Luescher
*
* adjusted for qcdlib: Tomasz Korzec 2016
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef H_QCD_RNGRANLUX
#define H_QCD_RNGRANLUX

/**
 * @fn int qcd_initRanlux(int seed)
 * @brief Initializes the ranlux random number generator with luxury level 2.
 * @details 
 * @param[in] seed determines the initial condition of ranlux
 *          
 */
int qcd_initRanlux(int seed);

/**
 * @fn void qcd_randRanlux(double *v, int N)
 * @brief generates double precision random numbers, flatly distributed between 0 and 1
 * @details 
 * @param[out] v pointer to allocated memory with space for N*sizeof(double)
 * @param[in] N number of random numbers to generate
 */
void qcd_randRanlux(double *v, int N);

/**
 * @fn void qcd_getRngTypeRanlux(int *rngType, int *n)
 * @brief Returns type and size of the ranlux RNG
 * @param[out] rngType #qcd_RNG_RANLUX
 * @param[out] n 
 */
void qcd_getRngTypeRanlux(int *rngType, int *n);

/**
 * @fn void qcd_getRngStateArrayRanlux(uint8_t *v)
 * @brief Returns a byte-array with the state of the RNG on the current process
 * @details The state of the ranlux random number generator is given by ...
 * @param v pointer to allocated space for at least ... bytes
 */
void qcd_getRngStateArrayRanlux(uint8_t *v);

/**
 * @fn void qcd_setRngStateArrayRanlux(uint8_t *v)
 * @brief initializes the ranlux RNG to the state encoded by a given byte-array.
 * @details The state of the ranlux random number generator is given by ...
 * @param v array with 6 bytes
 */
void qcd_setRngStateArrayRanlux(uint8_t *v);

#endif