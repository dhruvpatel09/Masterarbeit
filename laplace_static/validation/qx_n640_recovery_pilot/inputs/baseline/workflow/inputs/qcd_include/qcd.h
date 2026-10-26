/** @file qcd.h
 *  @brief main header for qcd-library
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2008-2019
 *  @copyright GNU Public License.
 */

/** @mainpage
 * @section intro Introduction
 * Lattice QCD simulations are computationally extremely demanding. Therefore
 * simulation programs are highly optimized. This often leads to data layouts
 * that are not very comfortable to use (even-odd ordering, cache blocks etc).
 * The focus of this library is to provide a set of simple to use functions
 * that allow quick development of "measurement software", i.e. programs that
 * read in existing gauge field configurations and evaluate complicated 
 * observables on them. For instance a baryon correlator. Some features that
 * are currently available are
 *  * Parallel I/O of gauge fields and propagators in various formats
 *  * Non-blocking communication
 *  * Simple manipulation of spinor fields
 *  * Various smearing functions (for correlation functions)
 *  * An interface to the (extremely fast) 
 *    <a href="http://luscher.web.cern.ch/luscher/openQCD/">openQCD</a> solvers.
 *  * The geometry is set up at run time. Various geometries may be used 
 *    simultaneously.
 *  * The library is written in C
 *  * A rather complete documentation exists
 *  
 * @section dep Dependencies
 * * For the I/O of ILDG files the lime library is required. It can be downloaded
 *   from <a href="http://usqcd-software.github.io/c-lime/">http://usqcd-software.github.io/c-lime/</a>
 *   follow the instructions given on the website to install it.
 *   This is OPTIONAL and the support of ILDG files can be turned on/off in qcd_config.h
 * * The inverter of openQCD is used. Currently the version 2.4.2 is required. The
 *   source code is available on
 *   <a href="http://luscher.web.cern.ch/luscher/openQCD/">http://luscher.web.cern.ch/luscher/openQCD/</a>.
 *   To compile it with qcdlib, see section @ref install.
 * 
 * 
 * @section install Compilation
 *  * Make sure you have a C compiler and MPI available.
 *  * Optionally: Download and compile <a href="http://usqcd-software.github.io/c-lime/">c-lime</a>.
 *    The pathts to liblime.a and to lime's include directory will be needed later.
 *  * Download <a href="http://www.qcdlib.com">qcdlib</a>. 
 *  * Edit the Makefile in ./build/. Necessary changes are
 *     - CC: your mpi C compiler, e.g. mpicc
 *     - CFLAGS: your prefered compiler options
 *     - QCD_DIR: path to qcdlib
 *     - QCD_INCLUDE_DIR: path to the qcdlib include files
 *     - QCD_LIB_DIR: path where to copy the compiled libqcd.a
 *     - QCD_TEST_DIR: path whith qcdlib test files
 *     - LIME_INCLUDE_DIR: path to lime's include folder
 *     - LIME_LIB_DIR: path to folder that contains liblime.a
 *  * compile with\n
 *    make clean\n
 *    make\n
 *    make test \#optional
 *  * install (i.e. copy the libqcd.a to its target folder)\n
 *    make install
 *  * download <a href="https://luscher.web.cern.ch/luscher/openQCD/openQCD-2.4.2.tar.gz">openQCD-2.4.2</a>
 *  * download the <a href="">solver interface</a>
 *  * Edit the Makefile in the solver interface's build folder
 *     -  
* 
 * @section ex Example
 * A simple program using qcdlib can look like this:
 * <a href="xxx">pion.c</a> \n
 * compilation will be similar to\n
 * mpicc -O2 -o pion pion.c -lqcd -lsolver -llime -lm
 * @section conventions Conventions
 * @subsection naming Naming conventions
 * * Every function/variable/preprocessor script/definition that is externally 
 *   visible starts with qcd_
 * * Functions, variables and definitions are (except for few exceptions) in lower case, e.g.
 *   qcd_vector
 * 
 * * If function, variable or definition names are composed of several words, capital 
 *   letters are used to differentiate them, as in qcd_gaugeField
 * 
 * * Pre-processor definitions and scripts are all capital letters, e.g.
 *   #qcd_CADD
 * 
 * @subsection functions Function calls
 * The generic case of a function in qcdlib looks like: \n
 * int qcd_someFuntionInQcdLib(*output1, *output2, *output3, ... , input1, input2, input3,...)\n
 * i.e. pointers to output arguments come first, input arguments follow. The return value
 * usually indicates whether the call was successful (see @ref errors).
 * @subsection errors Error handling
 * A library should in principle never terminate a program. In this spirit,
 * in case of failure, most functions write an error message to stderr and 
 * indicate failure
 * through the return value. It is up to the main program to decide whether to
 * exit or not. That being said, there are many exceptions to this
 * * MPI errors are fatal. qcdlib is based on MPI, so almost any function can
 *   terminate the program, if MPI encounters a problem
 * * openQCD errors are fatal. If the openQCD solver encounters an error, 
 *   the execution is terminated.
 * * Many functions don't test all the input arguments for sanity. Passing
 *   pointers to un-initialized structures will often result in segmentation
 *   faults and undefined behavior.
 *
 */

#ifndef H_QCD
#define H_QCD 1


// some constants:

/** @def qcd_MAX_STRING_LENGTH
 *  @brief Maximal allowed length for file-names etc.
 */
#define qcd_MAX_STRING_LENGTH 256

/** @struct qcd_complex_16 qcd.h
 *  @brief double precision complex numbers
 */
typedef struct
{
  double re; ///< real part
  double im; ///< imaginary part
} qcd_complex_16; 

/** @struct qcd_complex_8 qcd.h
 *  @brief single precision complex numbers
 */
typedef struct 
{
  float re; ///< real part
  float im; ///< imaginary part
} qcd_complex_8; 

/** @struct qcd_geometry qcd.h
 *  @brief Geometry structure
 *  @details Contains all information on 
 *  parallelization and geometry.
 */
typedef struct
{
   long V;                         ///< 4d volume
   long V3;                        ///< spatial volume
   int L[4];                       ///< lattice dimensions t,x,y,z = 0,1,2,3
   long lV;                        ///< local volume
   long lV3;                       ///< local spatial volume
   int lL[4];                      ///< dimensions of local lattice
   int Pos[4];                     ///< position of local lattice
   double theta[4];                ///< @brief boundary conditions of fermion fields
                                   ///< @details psi(x+L[mu]) = exp(i*theta[mu])*psi(x)
   int myid;                       ///< lexicographical index of position   
   int nproc;                      ///< number of processes
   int Pplus[8];                   ///< process IDs of neighbor blocks in + directions
   int *Pminus;                    ///< process IDs of neighbor blocks in - directions
   long (*plus)[8];                ///< @brief neighborhood relations on local lattice
                                   ///< @details x+mu = plus[x][mu] \n
                                   ///<          x-mu = plus[x][mu+4]
   long (*minus)[8];               ///< @brief neighborhood relations on local lattice
                                   ///< @details minus[x][0..3] == plus[x][4..7] \n
                                   ///< @warning minus[x][4..7] are undefined
   long (*plus3d)[8];              ///< @brief neighborhood relations on a local time-slice
                                   ///< @details x+mu = plus3d[x][mu] \n
                                   ///<          x-mu = plus3d[x][mu+4]
   long (*minus3d)[8];             ///< @brief neighborhood relations on local time-slice
                                   ///< @details minus3d[x][0..3] == plus3d[x][4..7] \n
                                   ///< @warning minus3d[x][4..7] are undefined
   long (*eo)[2];                  ///< @brief partitioning into even and odd 
                                   ///< @details eo[k][0/1] = lexicographical
                                   ///<  index of kth even/odd site
   long startBminus[4];      ///< index where - boundaries start
   long startBplus[4];       ///< index where + boundaries start
   long startBminus3d[4];    ///< index where - boundaries start (timeslice)
   long startBplus3d[4];     ///< index where + boundaries start (timeslice)
   long startEminus[4];      ///< inner field boundaries start index
   long startEplus[4];       ///< for communication to left (M) and right (P) neighbor
   long startEminus3d[4];    ///< inner 3d field boundaries start index
   long startEplus3d[4];     ///< for communication to left (minus) and right (plus) neighbor
   long boundarySize;        ///< number of points in outer boundaries
   long boundarySize3d;      ///< number of points in outer boundaries (timeslice)
   long *edge;               ///< lexicographical indices of inner boundaries
   long edgePoints;          ///< number of inner boundary points
   long *edge0;              ///< lexicographical indices of inner b. of a time-slice
   long edge0Points;         ///< number of inner time-slice boundary points
   int planeMu[6];           ///< @brief mu indices of the six planes
                             ///< @details There are 6 planes in 4d the nth (n=0..5) plane
                             ///< has mu<nu given by: planeMu[n],planeNu[n]
   int planeNu[6];           ///< @brief mu indices of the six planes
                             ///< @details There are 6 planes in 4d the nth (n=0..5) plane
                             ///< has mu<nu given by: planeMu[n],planeNu[n]
   char kbu2d;               ///< kbu2d=1 to keep outer boundaries up 2 date, 0 if not
   long lVb;                 ///< local volume possibly plus outer boundary size (if kbu2d=1)
   long lV3b;                ///< local spatial volume possibly plus outer boundary size (if kbu2d=1)

   MPI_Comm comm3d;          ///< 3d communicator
   int *root3d;              ///< rank 0 in comm3d is rank root3d[x0] in MPI_COMM_WORLD, x0=global time

   MPI_Request requests[32];       ///< needed in communication routines
   int numOfRequests;              ///< non 0 during ongoing communication
   MPI_Datatype stypeV[4];         ///< @brief Data-type that allows sends of strided data.
                                   ///< @details vector field type (for sending)
   MPI_Datatype rtypeV[4];         ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details vector field type (for receiving)
   MPI_Datatype stypeU[4];         ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Gauge-field type (for sending)
   MPI_Datatype rtypeU[4];         ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Gauge-field type (for receiving)
   MPI_Datatype stypeT[4];         ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Gauge-transformation type (for sending)
   MPI_Datatype rtypeT[4];         ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Gauge-transformation type (for receiving)
   MPI_Datatype stypeP[4];         ///< @brief Data-type that allows sends of strided data.
                                   ///< @details propagator type (for sending)
   MPI_Datatype rtypeP[4];         ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details propagator type (for receiving)
   MPI_Datatype stypeQ[4];         ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Clover field type (for sending)
   MPI_Datatype rtypeQ[4];         ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Clover field type (for receiving)
   MPI_Datatype stypeV3d[4];       ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Time-slice vector type (for sending)
   MPI_Datatype rtypeV3d[4];       ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Time-slice vector type (for receiving)
   MPI_Datatype stypeU3d[4];       ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Time-slice gaugeField type (for sending)
   MPI_Datatype rtypeU3d[4];       ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Time-slice gaugeField type (for receiving)
   MPI_Datatype stypeSC3d[4];      ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Time-slice spinor component type (for sending)
   MPI_Datatype rtypeSC3d[4];      ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Time-slice spinor component type (for receiving)
   MPI_Datatype stypeL3d[4];       ///< @brief Data-type that allows sends of strided data.
                                   ///< @details Time-slice loops type (for sending)
   MPI_Datatype rtypeL3d[4];       ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details Time-slice loops type (for receiving)
   int initializedEO;              ///< set to 1 if even/odd support enabled
   long startBminusEO[4];          ///< index where - boundaries E/O start
   long startBplusEO[4];           ///< index where + boundaries E/O start
   long boundarySizeEO;            ///< number of points in E/O outer boundaries
   long *edgeO;                    ///< lexicographical indices of inner O boundaries
   long *edgeE;                    ///< lexicographical indices of inner E boundaries
   long edgePointsEO;              ///< number of inner E/O boundary points
   long VEO;                       ///< 4d E/O volume
   int LEO[4];                     ///< E/O lattice dimensions t,x,y,z = 0,1,2,3
   long lVEO;                      ///< local E/O volume
   int lLEO[4];                    ///< dimensions of local E/O lattice
   long startEminusE[4];           ///< inner even field boundaries start index
   long startEplusE[4];            ///< for communication to left (minus) and right (plus) neighbor
   long startEminusO[4];           ///< inner odd field boundaries start index
   long startEplusO[4];            ///< for communication to left (minus) and right (plus) neighbor
   MPI_Datatype stypeUO[4];        ///< @brief Data-type that allows sends of strided data.
                                   ///< @details gaugeField (openQCD) type (for sending)
   MPI_Datatype rtypeUO[4];        ///< @brief Data-type that allows receivs of strided data.
                                   ///< @details gaugeField (openQCD) type (for receiving)
   int initialized;                ///< set to 1 during initialization 
} qcd_geometry;


/** @struct qcd_vector qcd.h
 *  @brief A spinor field.
 *  @details Dirac operators act on qcd_vector structures. They have
 *          space-time, dirac and color indices. The space-time index labels
 *          a point of the local lattice or its outer boundary.
 */
typedef struct
{
   qcd_complex_16 (*D)[4][3];         ///< data D[space-time][dirac][color]
   qcd_complex_16 (*Bplus[4])[4][3];  ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][3]; ///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                 ///< pointer to geometry
   int initialized;                   ///< set to 1 during initialization
} qcd_vector;


/** @struct qcd_vector3d qcd.h
 *  @brief A 3d spinor field.
 *  @details Like qcd_vector, but defined only on a single time-slice.
 */
typedef struct
{
   qcd_complex_16 (*D)[4][3];           ///< data D[space-time][dirac][color]
   qcd_complex_16 (*Bplus[4])[4][3];    ///< boundary data in + direction (only 1..3 are used)
   qcd_complex_16 (*Bminus[4])[4][3];   ///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date   
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_vector3d;


/** @struct qcd_spinorComponent3d qcd.h
 *  @brief A single dirac component of a time-slice spinor field.
 *  @details Has just space and color indices. Space is the local
 *           lattice and its outer boundaries.
 */
typedef struct
{
   qcd_complex_16 (*D)[3];              ///< data: D[space][color]
   qcd_complex_16 (*Bplus[4])[3];       ///< boundary data in + direction (only 1..3 are used)
   qcd_complex_16 (*Bminus[4])[3];      ///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_spinorComponent3d;


/** @struct qcd_propagator qcd.h
 *  @brief Structure to store a full propagator
 *  @details for propagator S(x,y) the whole N_dirac x N_dirac * N_color x N_color
 *          matrix is stored for a fixed y and all values of x of the local
 *          lattice and its outer boundaries.
 */
typedef struct
{
   qcd_complex_16 (*D)[4][4][3][3];         ///< data D[x][Dirac1][Dirac2][color1][color2]
   qcd_complex_16 (*Bplus[4])[4][4][3][3];  ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][4][3][3]; ///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                       ///< pointer to geometry
   int initialized;                         ///< set to 1 during initialization
} qcd_propagator;


/** @struct qcd_gaugeField qcd.h
 *  @brief Structure to store a gauge field
 *  @details For every link (x, mu) of the local lattice and its 
 *           outer boundaries, a 3x3 SU(3) matrix is stored.
 */
typedef struct
{
   qcd_complex_16 (*D)[4][3][3];        ///< data D[x][mu][color1][color2]
   qcd_complex_16 (*Bplus[4])[4][3][3]; ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][3][3];///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_gaugeField;

/** @struct qcd_gaugeFieldO qcd.h
 *  @brief Structure to store a gauge field
 *  @details This structure stores the gauge field in a fashion similar to
 *          openQCD. For every odd site x the gague field on all 8 links mu
 *          touching this site is stored
 */
typedef struct
{
   qcd_complex_16 (*D)[8][3][3];        ///< data D[x][mu][color1][color2]
   qcd_complex_16 (*Bplus[4])[8][3][3]; ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[8][3][3];///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_gaugeFieldO;

/** @struct qcd_gaugeField3d qcd.h
 *  @brief Structure to store a gauge field on a time slice
 *  @details Like qcd_gaugeField but defined only on a single time slice
 */
typedef struct
{
   qcd_complex_16 (*D)[4][3][3];        ///< data X[x3][mu][color1][color2]
   qcd_complex_16 (*Bplus[4])[4][3][3]; ///< boundary data in + direction (only 1..3 are used)
   qcd_complex_16 (*Bminus[4])[4][3][3];///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_gaugeField3d;

/** @struct qcd_loops3d qcd.h
 *  @brief Structure to store a product of gauge field along a loop on a time slice
 *  @details qcd_loops3d are 3x3 SU(3) matrices that transform like l(x) -> T'(x)l(x)T(x)
 */
typedef struct
{
   qcd_complex_16 (*D)[3][3];           ///< data X[x3][color1][color2]
   qcd_complex_16 (*Bplus[4])[3][3];    ///< boundary data in + direction (only 1..3 are used)
   qcd_complex_16 (*Bminus[4])[3][3];   ///< boundary data in - direction
   char BPup2date;                      ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                      ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_loops3d;

/** @struct qcd_gaugeTransformation qcd.h
 *  @brief Structure to store a 3x3 matrix on every site
 *  @details For site x on the local lattice or its outer boundary a 
 *           3x3 SU(3) matrix is stored.
 */
typedef struct
{
   qcd_complex_16 (*D)[3][3];        ///< data D[x][color1][color2]
   qcd_complex_16 (*Bplus[4])[3][3]; ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[3][3];///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                ///< pointer to geometry
   int initialized;                  ///< set to 1 during initialization
} qcd_gaugeTransformation;

/** @struct qcd_cloverField qcd.h
 *  @brief Structure to store the 'clover' field.
 *  @details For site x on the local lattice or its outer boundary 6
 *           3x3 matrices are stored, corresponding to 
 *           mu/nu = 0/1,0/2,0/3,1/2,1/3,2/3
 */
typedef struct
{
   qcd_complex_16 (*D)[6][3][3];        ///< data D[x][mu/nu][color1][color2]
   qcd_complex_16 (*Bplus[4])[6][3][3]; ///< boundary data in + direction
   qcd_complex_16 (*Bminus[4])[6][3][3];///< boundary data in - direction
   char BPup2date;                    ///< flag. 1 (0): positive outer boundaries are (not) up to date
   char BMup2date;                    ///< flag. 1 (0): negative outer boundaries are (not) up to date
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_cloverField;


/** @struct qcd_cloverDiagonal qcd.h
 *  @brief Structure to store the space-time diagonal part of a clover Dirac operator
 *  @details Due to the structure m +- i mu + hermitian
 *           not all entries need to be stored. 
 *           
 */
typedef struct
{
   qcd_complex_16 (*D)[2][21];
   /**< @brief Data. 
    * @details The matrix is block-diagonal in spin in D[x][pm][ind], x denotes
    *          the local lattice point, pm, which of the spin-blocks we are in
    *          and ind an index for the remaining 6x6 matrix, such that
    *\verbatim
      3*d1+c1,  3*d2+c2 --> ind
      0         0       --> 0
      0         1       --> 1
      0         2       --> 2
      ...
      0         5       --> 5
      1         1       --> 6
      1         2       --> 7
      ...
      1         5       --> 10
      ...
      5         5       --> 21
     \endverbatim 
    *  
    */
   qcd_geometry *geo;                   ///< pointer to geometry
   int initialized;                     ///< set to 1 during initialization
} qcd_cloverDiagonal;


//////////////////////////////////////////////////////////////////////////////
// some pre-processor macros. Use carefully!!
/** @def qcd_LEXIC(t,x,y,z,L)
 *  @brief Computes the lexicographical index. 
 *  @details Given a 4d point with coordinates (t,x,y,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + y*L[3] + x*L[2]*L[3] + t*L[1]*L[2]*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC(t,x,y,z,L) ( (long) ((z)+L[3]*((y)+L[2]*((x)+L[1]*(t)))) )

/** @def qcd_LEXIC0(x,y,z,L)
 *  @brief Computes the lexicographical index in outer boundaries in t direction
 *  @details Given a 3d point with coordinates (x,y,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + y*L[3] + x*L[2]*L[3] \n
 *           is computed.
 *           This is also the lexicographical index of time-slice fields
 *           like qcd_vector3d or qcd_gaugeField3d.
 */
#define qcd_LEXIC0(x,y,z,L) ( (long) ((z)+L[3]*((y)+L[2]*(x))) )

/** @def qcd_LEXIC1(t,y,z,L)
 *  @brief Computes the lexicographical index in outer boundaries in x direction
 *  @details Given a 3d point with coordinates (t,y,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + y*L[3] + t*L[2]*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC1(t,y,z,L) ( (long) ((z)+L[3]*((y)+L[2]*(t))) )

/** @def qcd_LEXIC2(t,x,z,L)
 *  @brief Computes the lexicographical index in outer boundaries in y direction
 *  @details Given a 3d point with coordinates (t,x,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + x*L[3] + t*L[1]*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC2(t,x,z,L) ( (long) ((z)+L[3]*((x)+L[1]*(t))) )

/** @def qcd_LEXIC3(t,x,y,L)
 *  @brief Computes the lexicographical index in outer boundaries in z direction
 *  @details Given a 3d point with coordinates (t,x,y) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = y + x*L[2] + t*L[1]*L[2] \n
 *           is computed.
 */
#define qcd_LEXIC3(t,x,y,L) ( (long) ((y)+L[2]*((x)+L[1]*(t))) )

/** @def qcd_LEXIC01(y,z,L)
 *  @brief Computes the lexicographical index in outer boundaries  in x direction
 *         of time-slice fields.
 *  @details Given a 2d point with coordinates (y,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + y*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC01(y,z,L)  ( (long) ((z)+L[3]*(y)) )

/** @def qcd_LEXIC02(x,z,L)
 *  @brief Computes the lexicographical index in outer boundaries in y direction
 *         of time-slice fields.
 *  @details Given a 2d point with coordinates (x,z) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = z + x*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC02(x,z,L)  ( (long) ((z)+L[3]*(x)) )

/** @def qcd_LEXIC03(x,y,L)
 *  @brief Computes the lexicographical index in outer boundaries in z direction
 *         of time-slice fields .
 *  @details Given a 2d point with coordinates (x,y) and the dimensions of
 *           the lattice L[0],..,L[3], the lexicographical index\n
 *           l = y + x*L[2] \n
 *           is computed.
 */
#define qcd_LEXIC03(x,y,L)  ( (long) ((y)+L[2]*(x)) )

/** @def qcd_LEXIC0_TO_LEXIC(v3,t,L)
 *  @brief Computes the full lexicographical index from a 3d index and
 *         a time.
 *  @details Given a 3D lexicographical index v3 (as computed by qcd_LEXIC0),
 *           a time coordinate t, and the lattice dimensions L[0]..L[4], the 
 *           full lexicographical index\n
 *           l = v3 + t*L[1]*L[2]*L[3] \n
 *           is computed.
 */
#define qcd_LEXIC0_TO_LEXIC(v3,t,L) ( (long) ((t)*L[1]*L[2]*L[3] + (v3) ))

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <qcd_config.h>
#include <qcd_gamma.h>
#include <qcd_communication.h>
#include <qcd_COps.h>
#include <qcd_3x3Ops.h>
#include <qcd_2x2Ops.h>
#include <qcd_blas.h>
#include <qcd_io.h>
#include <qcd_utils.h>
#include <qcd_smearing.h>
#include <qcd_gaugeFixing.h>
#include <qcd_wilson.h>
#include <qcd_observables.h>
#include <qcd_operators3d.h>
#include <qcd_loopsOps.h>
#include <qcd_lanczos.h>
#include <qcd_random.h>
#include <qcd_utils.h>

/* prototypes for qcd_init.c*/

/**
 * @fn void qcd_antilexic(int x[], long l, const int dim[])
 * @brief Given a position index l, the coordinates x[0]...x[3] are returned.
 * @param[out] x pointer to the coordinates
 * @param[in] l position index
 * @param[in] dim array of dimensions
 */
void qcd_antilexic(int x[], long l, const int dim[]);

/**
 * @fn qcd_antilexicOdd(int x[], long l, const int dim[])
 * @brief Given an odd-point-index l, the coordinates x[0]...x[3] are returned.
 * @param[out] x pointer to the coordinates
 * @param[in] l odd-point position index
 * @param[in] dim array of dimensions
 */
void qcd_antilexicOdd(int x[], long l, const int dim[]);

/**
 * @fn qcd_antilexicEven(int x[], long l, const int dim[])
 * @brief Given an even-point-index l, the coordinates x[0]...x[3] are returned.
 * @param[out] x pointer to the coordinates
 * @param[in] l even-point position index
 * @param[in] dim array of dimensions
 */
void qcd_antilexicEven(int x[], long l, const int dim[]);

/**
 * @fn qcd_initGeometry(qcd_geometry *geo, const int L[], const int P[], const double theta[], const int myid, const int numprocs)
 * @brief Initialize the main geometry structure.
 * @details Allocate memory, compute neighborhood relations, 
 *          calculate positions of outer 
 *          boundaries, define data types for communication.
 *          Needs to be called on all processes.
 * 
 * @param[out] geo pointer to the geometry structure that is initialized
 * @param[in] L vector of global lattice dimensions L[0]..L[3]
 * @param[in] P vector of processes to use in each direction P[0]..P[3]
 * @param[in] theta vector of angles for generalized periodic boundaries
 *        for fermion fields
 * @param[in] myid process ID
 * @param[in] numprocs number of processes
 * @return 0 on success, 1 on failure
 */
int qcd_initGeometry(qcd_geometry *geo, const int L[], const int P[], const double theta[], const int myid, const int numprocs);

/**
 * @fn void qcd_keepBndUp2Date(qcd_geometry *geo, char kbu2d)
 * @brief Control whether outer boundaries of fields are kept up to date when possible, or not.
 * @details If kbu2d=1, the loops over local volume in functions like addition,
 *          subtraction, copying etc. are extended to include the outer 
 *          boundaries.
 * 
 * @param[in,out] geo pointer to the geometry that is affected
 * @param[in] kbu2d flag. kbu2d=1 or kbu2d=0.
 */
void qcd_keepBndUp2Date(qcd_geometry *geo, char kbu2d);

/**
 * @fn qcd_destroyGeometry(const qcd_geometry *geo)
 * @brief Free all memory allocated by qcd_initGeometry.
 * @details Needs to be called on all
 *          processes.
 * @param[in] geo pointer to geometry structure to be destroyed
 */
void qcd_destroyGeometry(qcd_geometry *geo);

/**
 * @fn int qcd_initEO(qcd_geometry *geo)
 * @brief Add even/odd decomposition support to the geometry structure geo.
 * @details Needed for instance for I/O of openQCD gauge fields. Needs to be 
 *          called on all processes.
 * @param[in] geo pointer to a qcd_geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initEO(qcd_geometry *geo);

/**
 * @fn qcd_destroyEO(qcd_geometry *geo)
 * @brief Free all memory allocated by qcd_initEO.
 * @details Needs to be called on all
 *          processes.
 * @param[in] geo pointer to geometry structure from which e/o support is to 
 *            be removed.
 */
void qcd_destroyEO(qcd_geometry *geo);

/**
 * @fn int qcd_initVector(qcd_vector *vec, qcd_geometry *geo)
 * @brief Initializes a qcd_vector.
 * @details Allocates memory for a vector and its boundaries. Must be called on all 
 * processes before vec is used in any calculation or communication routines.
 * @param[out] vec pointer to a vector structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initVector(qcd_vector *vec, qcd_geometry *geo);

/**
 * @fn void qcd_destroyVector(qcd_vector *vec)
 * @brief Free all memory allocated by qcd_initVector.
 * @details Needs to be called on all processes.
 * @param[in] vec pointer to a vector structure
 */
void qcd_destroyVector(qcd_vector *vec);

/**
 * @fn int qcd_initVector3d(qcd_vector3d *vec, qcd_geometry *geo)
 * @brief Initializes a qcd_vector3d.
 * @details Allocates memory for a time-slice vector and its boundaries. 
 *          Must be called on all 
 *          processes before vec is used in any calculation 
 *          or communication routines.
 * @param[out] vec pointer to a 3d vector structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initVector3d(qcd_vector3d *vec, qcd_geometry *geo);

/**
 * @fn void qcd_destroyVector3d(qcd_vector3d *vec)
 * @brief Free all memory allocated by qcd_initVector3d.
 * @details Needs to be called on all processes.
 * @param[in] vec pointer to a 3d vector structure
 */
void qcd_destroyVector3d(qcd_vector3d *vec);

/**
 * @fn int qcd_initSpinorComponent3d(qcd_spinorComponent3d *spc, qcd_geometry *geo)
 * @brief Initializes a qcd_spiorComponent3d.
 * @details Allocates memory for a time-slice spinorComponent and its boundaries. 
 *          Must be called on all 
 *          processes before spc is used in any calculation 
 *          or communication routines.
 * @param[out] spc pointer to a 3d spinor component structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initSpinorComponent3d(qcd_spinorComponent3d *spc, qcd_geometry *geo);

/**
 * @fn void qcd_destroySpinorComponent3d(qcd_spinorComponent3d *spc)
 * @brief Free all memory allocated by qcd_initSpinorComponent3d.
 * @details Needs to be called on all processes.
 * @param[in] spc pointer to a 3d spinor component structure
 */
void qcd_destroySpinorComponent3d(qcd_spinorComponent3d *spc);

/**
 * @fn int qcd_initPropagator(qcd_propagator *prp, qcd_geometry *geo)
 * @brief Initializes a qcd_propagator.
 * @details Allocates memory for a propagator and its boundaries. 
 *          Must be called on all 
 *          processes before prp is used in any calculation 
 *          or communication routines.
 * @param[out] prp pointer to a propagator structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initPropagator(qcd_propagator *prp, qcd_geometry *geo);

/**
 * @fn void qcd_destroyPropagator(qcd_propagator *prp)
 * @brief Free all memory allocated by qcd_initPropagator.
 * @details Needs to be called on all processes.
 * @param[in] prp pointer to a propagator structure
 */
void qcd_destroyPropagator(qcd_propagator *prp);

/**
 * @fn int qcd_initGaugeField(qcd_gaugeField *u, qcd_geometry *geo)
 * @brief Initializes a qcd_gaugeField.
 * @details Allocates memory for a gauge field and its boundaries. 
 *          Must be called on all 
 *          processes before u is used in any calculation 
 *          or communication routines.
 * @param[out] u pointer to a gauge field structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initGaugeField(qcd_gaugeField *u, qcd_geometry *geo);

/**
 * @fn void qcd_destroyGaugeField(qcd_gaugeField *u)
 * @brief Free all memory allocated by qcd_initGaugeField.
 * @details Needs to be called on all processes.
 * @param[in] u pointer to a gauge field structure
 */
void qcd_destroyGaugeField(qcd_gaugeField *u);

/**
 * @fn int qcd_initGaugeFieldO(qcd_gaugeFieldO *u, qcd_geometry *geo)
 * @brief Initializes a qcd_gaugeFieldO.
 * @details Allocates memory for a gauge field (in openQCD layout)
 *          and its boundaries. 
 *          Must be called on all 
 *          processes before u is used in any calculation 
 *          or communication routines.
 * @param[out] u pointer to a gauge field structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initGaugeFieldO(qcd_gaugeFieldO *u, qcd_geometry *geo);

/**
 * @fn void qcd_destroyGaugeFieldO(qcd_gaugeFieldO *u)
 * @brief Free all memory allocated by qcd_initGaugeFieldO.
 * @details Needs to be called on all processes.
 * @param[in] u pointer to a gauge field structure
 */
void qcd_destroyGaugeFieldO(qcd_gaugeFieldO *u);

/**
 * @fn int qcd_initGaugeField3d(qcd_gaugeField3d *u, qcd_geometry *geo)
 * @brief Initializes a qcd_gaugeField3d.
 * @details Allocates memory for a time-slice gauge field 
 *          and its boundaries. 
 *          Must be called on all 
 *          processes before u is used in any calculation 
 *          or communication routines.
 * @param[out] u pointer to a time-slice gauge field structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initGaugeField3d(qcd_gaugeField3d *u, qcd_geometry *geo);

/**
 * @fn void qcd_destroyGaugeField3d(qcd_gaugeField3d *u)
 * @brief Free all memory allocated by qcd_initGaugeField3d.
 * @details Needs to be called on all processes.
 * @param[in] u pointer to a time-slice gauge field structure
 */
void qcd_destroyGaugeField3d(qcd_gaugeField3d *u);

/**
 * @fn int qcd_initLoops3d(qcd_loops3d *u, qcd_geometry *geo)
 * @brief Initializes a qcd_loops3d.
 * @details Allocates memory for a time-slice loops field 
 *          and its boundaries. 
 *          Must be called on all 
 *          processes before u is used in any calculation 
 *          or communication routines.
 * @param[out] u pointer to a time-slice loops field structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initLoops3d(qcd_loops3d *u, qcd_geometry *geo);

/**
 * @fn void qcd_destroyLoops3d(qcd_loops3d *u)
 * @brief Free all memory allocated by qcd_initLoops3d.
 * @details Needs to be called on all processes.
 * @param[in] u pointer to a time-slice loops structure
 */
void qcd_destroyLoops3d(qcd_loops3d *u);

/**
 * @fn int qcd_initCloverField(qcd_cloverField *clv, qcd_geometry *geo)
 * @brief Initializes a qcd_cloverField.
 * @details Allocates memory for a clover field
 *          and its boundaries. 
 *          Must be called on all 
 *          processes before clv is used in any calculation 
 *          or communication routines.
 * @param[out] clv pointer to a clover field structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initCloverField(qcd_cloverField *u, qcd_geometry *geo);

/**
 * @fn void qcd_destroyCloverField(qcd_cloverField *clv)
 * @brief Free all memory allocated by qcd_initCloverField.
 * @details Needs to be called on all processes.
 * @param[in] clv pointer to a clover field structure
 */
void qcd_destroyCloverField(qcd_cloverField *clv);

/**
 * @fn int qcd_initGaugeTransformation(qcd_gaugeTransformation *gt, qcd_geometry *geo)
 * @brief Initializes a qcd_gaugeTransformation.
 * @details Allocates memory for a gauge Transformation field
 *          and its boundaries. 
 *          Must be called on all 
 *          processes before clv is used in any calculation 
 *          or communication routines.
 * @param[out] gt pointer to a gauge transformation structure
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initGaugeTransformation(qcd_gaugeTransformation *gt, qcd_geometry *geo);

/**
 * @fn void qcd_destroyGaugeTransformation(qcd_gaugeTransformation *gt)
 * @brief Free all memory allocated by qcd_initGaugeTransformation.
 * @details Needs to be called on all processes.
 * @param[in] gt pointer to a gauge transformation field structure
 */
void qcd_destroyGaugeTransformation(qcd_gaugeTransformation *gt);

/**
 * @fn int qcd_initCloverDiagonal(qcd_cloverDiagonal *d, qcd_geometry *geo)
 * @brief Initializes a qcd_cloverDiagonal.
 * @details Allocates memory for a clover-diagonal field.
 *          Must be called on all 
 *          processes before d is used in any calculation routines.
 * @param[out] d pointer to a clover-diagonal field
 * @param[in] geo pointer to an initialized geometry structure
 * @return 0 on success, 1 on failure
 */
int qcd_initCloverDiagonal(qcd_cloverDiagonal *d, qcd_geometry *geo);
 
/**
 * @fn void qcd_destroyCloverDiagonal(qcd_cloverDiagonal *d)
 * @brief Free all memory allocated by qcd_initCloverDiagonal.
 * @details Needs to be called on all processes.
 * @param[in] d pointer to a clover-diagonal field structure
 */
void qcd_destroyCloverDiagonal(qcd_cloverDiagonal *d);
#endif
