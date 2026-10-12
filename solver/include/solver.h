/** @file solver.h
 *  @brief header file for solver.c
 *  @details Details & license
 *  @version 2.4.2
 *  @author Tomasz Korzec
 *  @date 2015-2026
 *  @copyright GNU Public License.
 */

#ifndef H_SOLVER
#define H_SOLVER 1


/** @def solver_NONE
 *  @brief label for "no solver" 
 */
#define solver_NONE 0
/** @def solver_CGNE
 *  @brief label for "conjugate gradient on normal equation solver" 
 */
#define solver_CGNE 1
/** @def solver_MSCG
 *  @brief label for "conjugate gradient solver with multiple mass shifts" 
 */
#define solver_MSCG 2
/** @def solver_SAP_GCR
 *  @brief label for "generalized conjugate residua solver with a Schwartz 
 *         alternating procedure as preconditioner" 
 */
#define solver_SAP_GCR 3
/** @def solver_DFL_SAP_GCR
 *  @brief label for "deflated generalized conjugate residua solver with a Schwartz 
 *         alternating procedure as preconditioner"
 */
#define solver_DFL_SAP_GCR 4

/** @def solver_ISPMAX 32
 *  @brief maximal number of definable solvers
 */
#define solver_ISPMAX 32

/** @struct solver_dflParameters
 *  @brief Deflation parameters.
 *  @details Gathers all paraameters that are needed for inexact deflation.
 *           These include the deflation subspace generation, deflation projection
 *           and deflation subspace parameters.
 */
typedef struct
{
   /* deflation parameters */
   int bs[4];  ///<  Block size
   int Ns;     ///< Number of basis functions in each block

   /* deflation projection */
   int nkv;        ///< Number of Krylov vectors, before restart
   int nmx;        ///< Maximal number of iterations
   double res;     ///< Relative residuum goal
   int nmx_gcr;    ///< Max Krylov vectors for precond. of little Dirac eq
   double res_gcr; ///< Desired rel. residue for precond.

   /* deflation subspace generation */
   int ninv;     ///< Number of inverse iterations
   int nmr;      ///< Number of SAP block solver iterations
   int ncy;      ///< Number of SAP cycles
   double kappa; ///< Hopping parameter
   double mu;    ///< Twisted mass parameter
} solver_dflParameters;

/** @struct solver_sapParameters
 *  @brief SAP parameters.
 */
typedef struct
{
   int bs[4]; ///< block size for Schwartz alternating procedure
} solver_sapParameters;

/** @struct solver_solverParameters
 *  @brief Solver parameters for Krylov space solvers.
 */
typedef struct
{
   int solver;             ///< Solver label, e.g. #solver_CGNE
   int nmx;                ///< Maximal number of iterations
   int nkv;                ///< Number of Krylov vectors before a restart
   int isolv;              ///< 1 if the block solver is E/O preconditioned, 0 else
   int nmr;                ///< Number of SAP block solver iterations
   int ncy;                ///< Number of SAP cycles
   int istop;              ///< stopping criterion 0:L_2 norm, 1: uniform norm based
   double res;             ///< Relative residuum goal
} solver_solverParameters;

/** @struct solver_parameters
 *  @brief Parameters of up to ISPMAX solvers and the common SAP and deflation parameters.
 */
typedef struct
{
  double kappa;                     ///< Hopping parameter
  int    isw;                       ///< type of sw-term. 0=traditional, 1=exponential
  double csw;                       ///< Improvement coefficient
  int bc;                           ///< Boundary conditions. 0=open,1=SF,2=open-SF,3=antiperiodic
  double cF;                        ///< Fermionic boundary improvement coefficient
  double cF_prime;                  ///< Fermionic boundary improvement coefficient
  double phi[2];                    ///< SF boundary angles phi1,phi2
  double phi_prime[2];              ///< SF boundary angles phi1',phi2'
  double theta[3];                  ///< fermionic spatial boundary phases  
  solver_solverParameters sol[solver_ISPMAX];  ///< Array with parameters for up to ISPMAX different solvers
  solver_sapParameters sap;         ///< Schwartz alternating procedure parameters
  solver_dflParameters dfl;         ///< Deflation projection/subspace generation parameters
} solver_parameters;









/* prototypes */
/** @fn solver_mpiInit(int argc,char **argv)
 * @brief replaces MPI_Init
 * @details Makes sure, that MPI works with openMP
 *
 * 
 * @param[in] argc number of command line parameters
 * @param[in] argv pointers to command line parameters
 */
void solver_mpiInit(int argc,char **argv);

/** @fn int solver_getNThread()
 * @brief returns number of openMP threads
 *
 * @return number of openMP threads
 */
int solver_getNThread(void);


/** @fn solver_init(qcd_geometry* geo, solver_parameters* par)
 * @brief Initializes the openQCD solvers
 * @details Given a qcdlib geometry and a collection of solver parameters,
 *          the following actions are carried out:\n
 * * initialize the openQCD geometry
 * * test for compatibility with qcdlib's geometry
 * * allocate space for openQCD gauge fields
 * * set up openQCD solvers
 * * if necessary: set up SAP and deflation
 * * allocate space for openQCD spinors\n
 * This function must be called from all processes, with 
 * identical values in the geometry and parameter structures
 *
 * 
 * @param[in] geo Pointer to an initialized qcdlib geometry structure
 * @param[in] par pointer to an initialized solver structure, that contains
 *            parameters for the ISPMAX solvers supported by openQCD
 * @return 0 upon success, 1 upon failure. Often, upon failure MPI_Abort is called 
 *           before the end of the function is reached
 */
int solver_init(qcd_geometry* geo, solver_parameters* par);


/**
 * @fn int solver_updateGaugeField(qcd_gaugeField* u)
 * @brief Updates the openQCD gauge field to the values given by a qcdlib gauge field.
 * @details The gauge field u is copied to the double precision
 *          gauge field of openQCD, and the modification flag is set.
 *          This function must be called from all processes
 *
 * @param[in] u pointer to an initialized openQCD gauge-field structure
 * @return 0 upon success, 1 upon failure
 * 
 */
int solver_updateGaugeField(qcd_gaugeField* u);

/**
 * @fn int solver_hypLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter)
 * @brief One iteration of HYP-smearing of the openQCD gauge field.
 * @details The HYP-smearing parameters alpha1, alpha2, alpha3 are used and the smearing is followed
 *          by a projection back to SU(3), that uses piter iterations. 
 *          All these steps are realized by calls to Björn Leder's
 *          functions.
 *          The smeared openQCD gauge field is finally copied to the qcdlib gauge field qcdU.
 * @param[in,out] qcdU pointer to an initialized qcdlib gauge-field structure
 * @param[in] alpha1 parameter for HYP smearing
 * @param[in] alpha2 parameter for HYP smearing
 * @param[in] alpha3 parameter for HYP smearing
 * @param[in] piter SU(3) projection uses piter iterations
 * @return 0 upon success, 1 upon failure
 * 
 */
int solver_hypLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter);

int solver_hypSpaceLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter);

int solver_hypTimeLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter);

int solver_hypLinks3d(qcd_gaugeField *qcdU, double alpha1, double alpha2, int piter);

void solver_wilsonFlow(int dnms, double eps);
double solver_computeWsl(double *wsl);
double solver_computeYsl(double *ysl);
double solver_computeQsl(double *qsl);
double **solver_allocFlowObs(int nn);
void solver_freeFlowObs(double **obs);


/**
 * @fn int solver_applyDw(qcd_vector *Dv, qcd_vector *v, double muTM)
 * @brief Applies the improved Wilson operator to a spinor field.
 * @details The Wilson (clover) operator with twisted mass muTM is 
 *          applied to the spinor v. The result is stored in Dv.
 *          This function must be called from all processes.
 * 
 * @param [in,out] Dv pointer to an initialized openQCD spinor structure
 * @param [in]      v pointer to an initialized openQCD spinor structure
 * @param [in] muTM twisted mass parameter
 * @return 0 upon success, 1 upon failure.
 *         Failure can cause an MPI_Abort before the end of the
 *         function is reached.
 *
 */
int solver_applyDw(qcd_vector *Dv, qcd_vector *v, double muTM);



/** 
 * @fn int solver_solveDw(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp)
 * @brief Solves a Dirac equation involving the improved Wilson operator.
 * @details The Dirac equation\n
 * \f$ [D_W+i\, \gamma_5\, \mathtt{muTM}]\, \mathtt{Dinvv} = \mathtt{v} \f$ \n
 * is solved using solver number isp. The result is stored
 * in Dinvv. Dw is the Wilson (clover) operator.
 * This function must be called from all processes simultaneously.
 * 
 * @param[in,out] Dinvv pointer to an initialized openQCD spinor structure
 * @param{in] v   pointer to an initialized openQCD spinor structure
 * @param{in} muTM  twisted mass parameter
 * @return 0 upon success, 1 upon failure. Failure can cause an MPI_Abort before the end of the 
 *           function is reached.
 *
 */
int solver_solveDw(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp);


/** 
 * @fn int solver_solveDwGuessed(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp)
 * @brief Solves a Dirac equation involving the improved Wilson operator.
 * @details The Dirac equation\n
 * \f$ [D_W+i\, \gamma_5\, \mathtt{muTM}]\, \mathtt{Dinvv} = \mathtt{v} \f$ \n
 * is solved using solver number isp. The result is stored
 * in Dinvv. Dw is the Wilson (clover) operator. An initial guess for the solution,
 * passed in Dinvv, is used.
 * This function must be called from all processes simultaneously.
 * 
 * @param[in,out] Dinvv pointer to an initialized openQCD spinor structure that
 *                contains the initial guess. It gets overwritten by the solution
 * @param{in] v   pointer to an initialized openQCD spinor structure
 * @param{in} muTM  twisted mass parameter
 * @return 0 upon success, 1 upon failure. Failure can cause an MPI_Abort before the end of the 
 *           function is reached.
 *
 */
int solver_solveDwGuessed(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp);

/**
 * @fn int solver_dflGeneration(void)
 * @brief Generates a basis for the deflation subspace.
 * @details A basis for the deflation subspace is generated
 *           with the current gauge field and the settings
 *           that were passed to #solver_init.
 *           This function must be called from all processes.
 * @return 0 upon success, 1 upon failure
 *
 */
int solver_dflGeneration(void);

int solver_setKappa(double kappa);
int solver_readDFL(solver_parameters *spar);
int solver_readSAP(solver_parameters *spar);
int solver_readSolver(solver_parameters *spar, int isp);
void solver_printSolvers(solver_parameters *par);
void solver_getProcessGrid(int *p);
void solver_getLatticeSize(int *l);
int solver_getOpenQCDGaugeField(qcd_gaugeField *u);

/**
 * @fn void solver_getStatus(int *stat)
 * @brief Provides a copy of the openQCD status array
 * @details stat must be an array with at least 4 elements.
 * After calls to the solver or to deflation subspace generation
 * solver_getStatus can be used to retrieve the openQCD status of 
 * that operation. After solves these are the iteration counts
 * or error codes.
 */

void solver_getStatus(int *stat);


/**
 * @fn int solver_free(void)
 * @brief Frees up space reserved by #solver_init.
 * @details Memory for the global single and double precision 
 * gauge fields and for the workspace is freed.
 * This function must be called from all processes.
 * @return 0 upon success, 1 upon failure
 */
int solver_free(void);
#endif
