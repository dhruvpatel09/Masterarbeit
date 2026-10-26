/** @file qcd_lanczos.h
 *  @brief header file for qcd_lanczos.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2011-2015
 *  @copyright GNU Public License.
 */


#ifndef H_QCD_LANCZOS
#define H_QCD_LANCZOS 1

/**
 * @def QCD_LANCZOS_PRECISION_GOAL
 * @brief Stopping criterion for the Lanczos process
 */
#define QCD_LANCZOS_PRECISION_GOAL 1e-9


/* structures */

/** @struct qcd_lanczosParams
 *  @brief Parameters used in Lanczos algorithms
 *  @details qcd_lanczosReortho will use only the fields res and nLanczosFactor,
 *           qcd_lanczosCheby will use all fields.
 */
typedef struct
{
   // input parameters
   double res;           ///< precision goal: all eigenvectors ||A*vec-lambda*vec||/lambda <= res
   double redTol;        ///< precision goal for reduced eigensystems
   int deg;              ///< degree of Chebyshev polynomial
   double upper;         ///< Everything between [upperLimit, 15] gets supressed so desired eigenvalues should be smaller than upper
   double tolReortho;    ///< Tolerance for loss of orthogonality in selective reorthogonalization
   int maxRestarts;      ///< Maximal number of allowed restarts.
   int nLanczosFactor;   ///< Max. number of iterations = number of desired eigenvectors times nLanczosFactor
   int redFreqFactor;    ///< solve reduced eigensystem when iteration is multiple of redFreqFactor times number of desired eigenvectors
} qcd_lanczosParams;

/** @struct qcd_lanczosStatus
 *  @brief Diagnostic output of last Lanczos call
 *  @details qcd_lanczosReortho will use only the field res,
 *           qcd_lanczosCheby will use all fields.
 */
typedef struct
{
   // diagnostic output of last Lanczos call
   int restarts;          ///< number of restarts
   int nRed;              ///< number of solved reduced eigensystems
   int maxRed;            ///< largest small eigensystem solved
   int usedLapack;        ///< 1 if lapacke was used, 0 otherwise
   double maxRes;         ///< maximal relative residuum
   int solverType;        ///< what eigenvalue solver has been used: 1=qcd_lanczosReortho, 2=lanczosCheby, otherwise: none 
} qcd_lanczosStatus;



/* prototypes */
/**
 * @fn int qcd_fsign(double f)
 * @brief Returns 1 if f>=0, and -1 otherwise.
 * @param[in] f
 * @return sign(f)
 */
int qcd_fsign(double f);

/**
 * @fn void qcd_printLanczosStatus(qcd_lanczosStatus stat)
 * @brief summarizes the result of the last lanczos call
 * @param[in] stat lanczos status structure
 */
void qcd_printLanczosStatus(qcd_lanczosStatus stat);

/**
 * @fn void qcd_bubblesort(int *I, double *a, int n)
 * @brief Sorts an array of double precision numbers.
 * @details Given an array a of length n, the index I is computed such that \n
 *          a[I[0]] <= a[I[1]] <= a[I[2]] <= ... <= a[I[n-1]] \n
 *
 * @todo the rather primitive bubble-sort should be replaced
 *       by something efficient in the future
 *
 * @param[out] I pointer to the idnex array. Needs to have space for n integers
 * @param[in] a input array of length n
 * @param[in] n length of the input/output arrays
 */
void qcd_bubblesort(int *i, double *a, int n);

/**
 * @fn void qcd_givens(double *c, double *s, double a, double b)
 * @brief A Givens rotation.
 * @details Computes c and s such that \n
 * \verbatim
      [ c  s ]T  [ a ]   [ r ]
      [-s  c ]   [ b ] = [ 0 ]
   \endverbatim
 * Givens-rotation as in Golub & van Loan p. 216
 *
 * @param[out] c pointer to cos(phi) of the Givens rotation
 * @param[out] s pointer to sin(phi) of the Givens rotation
 * @param[in] a first component of the input vector
 * @param[in] b second component of the input vector
 */
void qcd_givens(double *c, double *s, double a, double b);

/**
 * @fn void qcd_implicitSymmetricQRStepWithShift(double *alpha, double *beta, double **Q, int N, int start, int stop)
 * @brief implicit Symmetric QR Step with Wilkinson shift
 * @details As in Golub & van Loan algorithm 8.3.2 \n
 * alpha and beta define a tridiagonal N x N matrix
 * \verbatim
     T=  [alpha[0] beta[0]  0...            ]
         [beta[0]  alpha[1] beta[1]  0 ...  ]
         [0        .        .        .      ]
         [0               beta[n-1] alpha[n]]
   \endverbatim
 * Q is a full N x N matrix \n
 * start, stop defines sub matrix \n
 * t = T(start:stop, start:stop) \n
 * and the algorithm overwrites t with Z^T t Z with Z being givens rotations with property
 * Z^T(t-mu) is upper triangular and mu is eigenvalue of t's trailing 2x2 principal submatrix
 * Q is updated: Q <- Q diag(1,Z,1)
 * @param[in] alpha array with diagonal elements of a tridiagonal symmetric matrix T
 * @param[in] beta array with offdiagonal elements of a tridiagonal symmetric matrix T
 * @param[in,out] Q N x N matrix to be updated
 * @param[in] N size of Q
 * @param start defines a submatrix T(start:stop, start:stop)
 * @param stop defines a submatrix T(start:stop, start:stop)
 */
void qcd_implicitSymmetricQRStepWithShift(double *alpha, double *beta, double *Q, int N, int start, int stop);

/**
 * @fn qcd_symmetricQR(double *al, double *be, int N, double *eval, double *Q, double tol)
 * @brief Computes Eigenvectors of symmetric tridiagonal matrices.
 * @details Symmetric QR algorithm as in Golub & van Loan alg. 8.3.3 
 * alpha, beta: diagonal & sub/super diagonal elements of a symmetric matrix. 
 * Both have length N.\n
 * N eigenvectors are stored in Q (N x N). The result is sorted in ascending order
 * @param[in] alpha array with diagonal elements of a tridiagonal symmetric matrix
 * @param[in] beta array with offdiagonal elements of a tridiagonal symmetric matrix
 * @param[out] Q the columns of Q are the eigenvectors of the tridiagonal matrix
 * @param[in] N size of the problem
 * @param[in] tol absolute tolerance (off diagonals have to drop below tol)
 */
void qcd_symmetricQR(double *al, double *be, int N, double *eval, double *Q, double tol);
//void qcd_symmetricQR(double *alpha, double *beta, double **Q, int N);


/**
 * @fn void qcd_defaultLanczosParams(qcd_lanczosParams *par)
 * @brief create a Lanczos parameter structure with default values
 * @details 
 * @author Tomasz Korzec
 * @param[in,out] pointer to parameter structure that will be initialized with default values
 */
void qcd_defaultLanczosParams(qcd_lanczosParams *par);

/**
 * @fn int qcd_lanczosCheby(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par)
 * @brief Lanczos algorithm with selective re-orthogonalization and Chebyshev accelleration.
 * @details Calculates the lowest n eigenvalues eval and eigenvectors evec of
 *          the 3d covariant Laplacian
 *          u is the gauge-field of the slice. \n
 *
 *          At most the first par.nLanczosFactor*n Lanczos iterations with complete re-orthogonalization
 *          are carried out (memory!)\n
 * 
 *          Eigenvectors are normalized and phase conventions such that the first entry
 *          is real and positive.
 *
 *          The quality of an eigenvalue/eigenvector pair is \n
 *          q = ||A*evec - eval*evec||/eval \n
 *          If the biggest q out of the first n is larger than the precision goal an error code is returned.
 *           All processes in a 3d communicator have
 *          to call this function simultaneously.
 * @author Juan Andres Urrea Nino, adapted to qcdlib by Tomasz Korzec
 * @param[out] eval array of eigenvalues
 * @param[in,out] evec array of n initialized qcd_spinorComponent3d structures
 * @param[out] stat pointer to diagnostic structure. Can be NULL.
 * @param[in] u pointer to an initialized qcd_gaugeField3d structure
 * @param[in] nvec number of eigenvalues/eigenvectors to compute
 * @param[in,out] par input parameters for lanczos, and diagnostic data for last Lanczos run
 * @return 0 on success, error-code else: 1 out-of-memory, 2 no convergence, 3 error in DSYEVR, 4 error in DSTEMR, 5 upper too small
 */
int qcd_lanczosCheby(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par);

/**
 * @fn int qcd_lanczosReortho(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par)
 * @brief Lanczos algorithm with complete re-orthogonalization.
 * @details Calculates the lowest n eigenvalues eval and eigenvectors evec of
 *          the 3d covariant Laplacian with mass-shift mm
 *          u is the gauge-field of the slice. \n
 *
 *          The first par.nLanczosFactor*n Lanczos iterations with complete re-orthogonalization
 *          are carried out (memory!)\n
 *
 *          The quality of an eigenvalue/eigenvector pair is \n
 *          q = (A*evec - eval*evec) \n
 *          If the biggest q out of the first n is larger than the precision goal, an error is returned.
 *          All processes have
 *          to call this function simultaneously.
 *
 * @param[out] eval array of eigenvalues
 * @param[in,out] evec array of n initialized qcd_spinorComponent3d structures
 * @param[out] stat pointer to diagnostic structure. Can be NULL.
 * @param[in] u pointer to an initialized qcd_gaugeField3d structure
 * @param[in] nvec number of eigenvalues/eigenvectors to compute
 * @param[in,out] par input parameters for lanczos, and diagnostic data for last Lanczos run
 * @return  0 on success, error-code else: 1 out-of-memory, 2 no convergence, 3 error in DSYEVR, 4 error in DSTEMR
 */
int qcd_lanczosReortho(double *eval, qcd_spinorComponent3d *evec, qcd_lanczosStatus *stat, qcd_gaugeField3d *u, int nvec, qcd_lanczosParams *par);
#endif
