/** @file qcd_random.h
 *  @brief header file for qcd_random.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2016
 *  @copyright GNU Public License.
 */

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#ifndef H_QCD_RANDOM
#define H_QCD_RANDOM 1

#include<stdint.h>

/**
 * @def qcd_RNG_DRAND48
 * @brief RNG type: linear congruential generator
 */
#define qcd_RNG_DRAND48 0

/**
 * @def qcd_RNG_RANLUX
 * @brief RNG type: ranlux with luxury level 2
 */
#define qcd_RNG_RANLUX 10

/**
 * @def qcd_RNG_RANLUX_LVL1
 * @brief RNG type: ranlux with luxury level 1
 */
#define qcd_RNG_RANLUX_LVL1 11

/**
 * @def qcd_RNG_MLF
 * @brief RNG type: multiplicative lagged Fibonacci
 */
#define qcd_RNG_MLF 20

/**
 * @def qcd_NOISE_UNIFORM
 * @brief Noise type: uniform distribution between 0 and 1
 */
#define qcd_NOISE_UNIFORM 3

/**
 * @def qcd_NOISE_U1
 * @brief Noise type: U(1) valued noise
 */
#define qcd_NOISE_U1 0

/**
 * @def qcd_NOISE_Z2
 * @brief Noise type: Z2 valued noise
 */
#define qcd_NOISE_Z2 1

/**
 * @def qcd_NOISE_GAUSS
 * @brief Noise type: real, normal distributed noise, mu=0, sig=1
 */
#define qcd_NOISE_GAUSS 2

/**
 * @def qcd_NOISE_Z4
 * @brief Noise type: Z4 valued noise
 */
#define qcd_NOISE_Z4 4


/** @struct qcd_rngState
 *  @brief Structure to store the state of a RNG.
 *  @details The type, size and state of a random number generator can be stored
 *           in this structure. The state is a collection of states on all
 *           processes.
 */
typedef struct
{
    int rng;             ///< Type of the RNG, e.g. qcd_RNG_DRAND48
    int size;            ///< Size of the state vector in bytes (single process)
    int nproc;           ///< number of processes
    uint8_t **dat;       ///< collection of state vectors from all processes
} qcd_rngState;


#include <qcd_rngDrand48.h>
#include <qcd_rngRanlux.h>

/**
 * @fn char* qcd_getNoiseName(int noisetype)
 * @brief Returns the name of the noise for a given integer noise-ID
 * @param[in] noisetype noise-ID, e.g. qcd_NOISE_U1
 * @return pointer to a string with the noise name on success, NULL pointer on 
 *         failure
 */
char* qcd_getNoiseName(int noisetype);

/**
 * @fn int qcd_initRng(int rngType, int seed)
 * @brief Initializes a psudo random number generator on all processes.
 * @details The pseudo random number generator is set to type rngType, e.g.
 *          qcd_RNG_RANLUX. Each process gets a different seed, that is computed
 *          from seed and the process id.
 * @param[in] rngType type of pseudo random number generator to use
 * @param[in] seed value that is used to compute the initial state on all processes
 * @return 0 on success, 1 on failure
 */
int qcd_initRng(int rngType, int seed);

/**
 * @brief function pointer to the RNG in use.
 * @details The random number generator of choice will be called to fill N entries
 *          in the array v with numbers that are flatly distributed between 0 and
 *          1.
 * @param[out] v points to an allocated memory, large enough to hold 
 *             N*sizeof(double) bytes
 * @param[in] N number of random numbers to generate
 */
extern void (*qcd_rand)(double *v, int N);

/**
 * @brief function pointer to a function that returns the state-array of a single process.
 * @details Each RNG has a function that returns its state as a string of bytes.
 *          qcd_getRngStateArray points to this function of the currently used
 *          RNG.
 * @param[out] v pointer to allocated memory with enough space to accommodate a state array 
 *             of the current RNG
 */
extern void (*qcd_getRngStateArray)(unsigned char *v);

/**
 * @brief function pointer to a function that sets the state-array on a single process.
 * @details Each RNG has a function sets its state to a string of bytes.
 *          qcd_setRngStateArray points to this function of the currently used
 *          RNG.
 * @param[out] v string of bytes that describes the state of the RNG
 */
extern void (*qcd_setRngStateArray)(unsigned char *v);

/**
 * @brief function pointer a function that returns type and size of current RNG
 * @details Random number generators have functions that return their type, e.g. 
 *          #qcd_RNG_DRAND48, and the size of their state per process in byte.
 *          These information are written into rngType and n.
 * @param[out] rngType type of the RNG
 * @param[out] n size of the state in byte (per process)
 */
extern void (*qcd_getRngType)(int *rngType, int *n);

/**
 * @fn void qcd_randGauss(double *v, int N)
 * @brief Generate random numbers with Gauss distribution.
 * @details The array v is filled with N random numbers drawn from a Gaussian
 *        distribution with mean zero and variance one. The
 *        <a href = "https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform">
 *        Box-Muller-Transform</a> is used to map pairs of flatly distributed 
 *        random numbers, created with the current RNG, onto pairs of 
 *        Gauss distributed ones.
 *        If N is odd, there is a small overhead.
 * @param v pointer to allocated memory with space for N*sizeof(double)
 * @param N number of gaussian random numbers to generate
 * @warning if qcd_rand creates an exact 0, qcd_randGauss will fail
 */
void qcd_randGauss(double *v, int N);

/**
 * @fn void qcd_randZ2(double *v, int N)
 * @brief Generate random numbers with Z2 distribution
 * @details The array v is filled with N random numbers drawn from a Z2 
 *          distribution, i.e. the values are -1 or +1 with equal probability.
 * @param v pointer to allocated memory with space for N*sizeof(double)
 * @param N number of random Z2 numbers to generate
 */
void qcd_randZ2(double *v, int N);

/**
 * @fn void qcd_randZ4(qcd_complex_16 *v, int N)
 * @brief Generate random numbers with Z4 distribution
 * @details The array v is filled with N random numbers drawn from a Z4
 *          distribution, i.e. the values are -1,+1,-i,+i with equal probability.
 * @param v pointer to allocated memory with space for N*sizeof(qcd_complex_16)
 * @param N number of random Z4 numbers to generate
 */
void qcd_randZ4(qcd_complex_16 *v, int N);


/**
 * @fn void qcd_randU1(qcd_complex_16 *v, int N)
 * @brief Generate random numbers with U(1) distribution
 * @details The complex array v is filled with N random numbers drawn from a
 *          U(1) distribution, i.e. with values exp(2*pi*i*phi), where phi
 *          is uniformly distributed between 0 and 1.
 * @param v pointer to allocated memory with space for N*sizeof(qcd_complex_16)
 * @param N number of random U(1) numbers to generate
 */
void qcd_randU1(qcd_complex_16 *v, int N);



/**
 * @fn int qcd_getRngState(qcd_rngState *state)
 * @brief Collects the RNG state on all processes and stores it in state.
 * @details The type and size of the current random number generator as well as
 *          its state on all the processes is collected and stored in state.
 *          Space for the state vectors state->dat is allocated during the 
 *          process and can be released with #qcd_freeRngState.
 *          All processes must call this function simultaneously. All processes
 *          will obtain the same result.
 * @param state pointer to a qcd_rngState structure
 * @return 0 on success, 1 on failure
 */
int qcd_getRngState(qcd_rngState *state);

/**
 * @fn int qcd_setRngState(qcd_rngState *state)
 * @brief The RNG is initialized to a given state
 * @details The values in state on process 0 are used to set the
 *          type of the random generator and 
 *          its states on all processes.
 *          All processes must call this function simultaneously.
 * @param state pointer to a properly filled qcd_rngState structure
 * @return 0 on success, 1 on failure
 */
int qcd_setRngState(qcd_rngState *state);

/**
 * @fn void qcd_freeRngState(qcd_rngState *state)
 * @brief Releases the space used for RNG state arrays.
 * @param state pointer to a properly filled qcd_rngState structure
 */
void qcd_freeRngState(qcd_rngState *state);

/**
 * @fn void qcd_randomVector(qcd_vector *v, int noisetype)
 * @brief Create a random spinor field
 * @details The (local) spinor field v is filled with random numbers from a
 *          distribution given by noisetype. Currently supported are:
 *          * qcd_NOISE_UNIFORM: real, uniformly distributed numbers between 0 and 1
 *          * qcd_NOISE_Z2: real numbers with values -1 or +1
 *          * qcd_NOISE_U1: complex phases exp(i*2*pi*phi) with uniformly 
 *                          distributed phase angles phi between 0 and 1
 *          * qcd_NOISE_GAUSS: real, normally distributed numbers with mean 0 and
 *                             variance 1
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in] noisetype type of noise to fill v with. E.g. qcd_NOISE GAUSS
 */
void qcd_randomVector(qcd_vector *v, int noisetype);

/**
 * @fn void qcd_randomVector3d(qcd_vector3d *v, int noisetype)
 * @brief Create a random 3d spinor field
 * @details The (local) 3d spinor field v is filled with random numbers from a
 *          distribution given by noisetype. Currently supported are:
 *          * qcd_NOISE_UNIFORM: real, uniformly distributed numbers between 0 and 1
 *          * qcd_NOISE_Z2: real numbers with values -1 or +1
 *          * qcd_NOISE_U1: complex phases exp(i*2*pi*phi) with uniformly 
 *                          distributed phase angles phi between 0 and 1
 *          * qcd_NOISE_GAUSS: real, normally distributed numbers with mean 0 and
 *                             variance 1
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @param[in] noisetype type of noise to fill v with. E.g. qcd_NOISE GAUSS
 */
void qcd_randomVector3d(qcd_vector3d *v, int noisetype);


/**
 * @fn void qcd_randomGaugeField(qcd_gaugeField *u)
 * @brief Create a random gauge field
 * @details on every link a 3x3 matrix with random entries is created. 
 *          Both real and imaginary part of each matrix element are 
 *          normally distributed. The matrix is then projected to 
 *          SU(3).
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 */
void qcd_randomGaugeField(qcd_gaugeField *u);

/**
 * @fn void qcd_randomGaugeTransformation(qcd_gaugeTransformation *g)
 * @brief Create a random gauge transformation field
 * @details on every site a 3x3 matrix with random entries is created. 
 *          Both real and imaginary part of each matrix element are 
 *          normally distributed. The matrix is then projected to 
 *          SU(3).
 * @param[in,out] g pointer to an initialized qcd_gaugeTransformation structure
 */
void qcd_randomGaugeTransformation(qcd_gaugeTransformation *g);
#endif 
