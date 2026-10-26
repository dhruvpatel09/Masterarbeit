#ifndef H_QCD_CONFIG
#define H_QCD_CONFIG 1

/** @def qcd_NEEDS_LIME
 *  @brief If defined, the ildg configuration file format is supported, which requires the lime library.
 */
/* #define qcd_NEED_LIME */
#undef qcd_NEEDS_LIME


/** @def qcd_KEEP_BND_UP2DATE
 *  @brief If defined, the default behavior is to keep outer boundaries up-to-date whenever possible,
 *         i.e. operations like additions and subtractions also affect outer boundary points.
 */
#define qcd_KEEP_BND_UP2DATE 1
/*#define qcd_KEEP_BND_UP2DATE 0 */

/** @def qcd_USE_CBLAS
 *  @brief If defined, the default behavior is to use cblas for basic linear algebra
 *         operations.
 */
/*#define qcd_USE_CBLAS 1*/
#undef qcd_USE_CBLAS


/** @def qcd_USE_LAPACKE
 *  @brief If defined, the default behavior is to use cblas for basic linear algebra
 *         operations, implies qcd_USE_CBLAS.
 */
#define qcd_USE_LAPACKE 1
/* #undef qcd_USE_LAPACKE */


/** @def qcd_USE_MKL
 *  @brief If defined, cblas and lapacke from intel's math-kernel library are used.
 *         Implies qcd_USE_CBLAS and qcd_USE_LAPACKE
 */
/*#define qcd_USE_MKL 1*/
#undef qcd_USE_MKL



/** @def qcd_MEMORY_ALIGNMENT
 *  @brief Allocated memory for data structures is aligned to 2^qcd_MEMORY_ALIGNMENT 
 *         byte boundaries.
 */
#define qcd_MEMORY_ALIGNMENT 4

/** @def qcd_HAS_MALLINFO
 * @brief if defined, mallinfo is used to provide mem usage statistics via qcd_printMemInfo
 */
/*#define qcd_HAS_MALLINFO 1 */

/** @def qcd_HAS_MALLINFO2
 * @brief if defined, mallinfo is used to provide mem usage statistics via qcd_printMemInfo.
 *        use this with gclib 2.33 and newer
 */
/*#define qcd_HAS_MALLINFO2 1*/



/** @def qcd_ALPHA_SU3_PROJECTION
 *  @brief If defined, the iterative method of SU(3) projection, as used by the
 *         alpha collaboration,  is employed
 * 
 *  @warning: this does not seem to work
 */
//#define qcd_ALPHA_SU3_PROJECTION 1
//#define qcd_ALPHA_SU3_PROJECTIN_ITERATIONS 5
#undef qcd_ALPHA_SU3_PROJECTION

/** @def qcd_ETMC_SU3_PROJECTION
 *  @brief If defined, the method of SU(3) projection, as used by the
 *         ETMC collaboration,  is employed
 * 
 */
#define qcd_ETMC_SU3_PROJECTION



/************************ do not edit beyond this point *****************************/
#ifdef qcd_USE_LAPACKE
   #define qcd_USE_CBLAS 1
#endif
#ifdef qcd_USE_MKL
   #define qcd_USE_CBLAS 1
   #define qcd_USE_LAPACKE 1
#endif


#endif
