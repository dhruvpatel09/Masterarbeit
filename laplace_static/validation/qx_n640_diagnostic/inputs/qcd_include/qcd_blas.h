/** @file qcd_blas.h
 *  @brief header file for qcd_blas.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2022
 *  @copyright GNU Public License.
 *  @todo Check in all functions whether parameters can be declared "const"
 *  @todo qcd_raxpy_vector
 *  @todo axpy and raxpy for other types: propagators, gauge-fields clover etc
 *  @todo 3x3: inverse, eigenvalues, det
 */
 
#ifndef H_QCD_BLAS
#define H_QCD_BLAS 1 
 
 
 




/******************************************************************************
 * prototypes 
 ******************************************************************************/



/**
 * @fn void qcd_zeroVector(qcd_vector *vec)
 * @brief Sets a spinor field to zero.
 * @details Sets the spinor field vec to zero on the local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @todo add initialization check
 */
void qcd_zeroVector(qcd_vector *vec);


/**
 * @fn qcd_unitVector(qcd_vector *vec, int x0, int x1, int x2, int x3, int alpha, int a);
 * @brief Creates a point source
 * @details Sets the spinor field vec to zero on the local lattice.
 *          On global coordinate x0,x1,x2,x3 the real part of the spin component alpha 
 *          and color component a is set to 1.0
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] x0 global time coordinate
 * @param[in] x1 glabal x coordinate
 * @param[in] x2 global y coordinate
 * @param[in] x3 global z coordinate
 * @param[in] alpha spin component 0-3
 * @param[in] a color component 0-2
 * @todo Maybe set outer boundaries correctly if the flag is set.
 * @todo add initialization check
 */
void qcd_unitVector(qcd_vector *vec, int x0, int x1, int x2, int x3, int alpha, int a);

/**
 * @fn void qcd_zeroVector3d(qcd_vector3d *vec)
 * @brief Sets a 3d spinor field to zero.
 * @details Sets the time-slice spinor field vec to zero on the local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector3d structure
 * @warning The outer boundaries are not set to zero.
 * @todo Maybe (?) set also the outer boundaries to 0
 * @todo add initialization check
 */
void qcd_zeroVector3d(qcd_vector3d *vec);

/**
 * @fn void qcd_zeroSpinorComponent3d(qcd_spinorComponent3d *sp)
 * @brief Sets a 3d spinor-component field to zero.
 * @details Sets the time-slice spinor-component field sp to zero on the local
 *          lattice.
 * @param[in,out] sp pointer to an initialized 3d spinor-component structure
 * @warning The outer boundaries are not set to zero.
 * @todo Maybe (?) set also the outer boundaries to 0
 * @todo add initialization check
 */
void qcd_zeroSpinorComponent3d(qcd_spinorComponent3d *sp);

/**
 * @fn void qcd_zeroPropagator(qcd_propagator *prop)
 * @brief Sets a propagator to zero.
 * @details Sets the propagator prop to zero on the local lattice.
 * @param[in,out] prop pointer to an initialized qcd_propagator structure
 * @warning The outer boundaries are not set to zero.
 * @todo Maybe (?) set also the outer boundaries to 0
 * @todo add initialization check
 */
void qcd_zeroPropagator(qcd_propagator *prop);

/**
 * @fn void qcd_zeroGaugeField(qcd_gaugeField *u)
 * @brief Sets a gauge field to zero.
 * @details Sets a gauge field u to zero 3x3 matrix on the local lattice.
 * @param u pointer to an initialized gauge field structure
 * @warning The outer boundaries are not set to zero.
 * @todo Maybe (?) set also the outer boundaries to 0
 * @todo add initialization check
 */
void qcd_zeroGaugeField(qcd_gaugeField *u);

/**
 * @fn void qcd_unitGaugeField(qcd_gaugeField *u)
 * @brief Sets a gauge field to the unit matrix.
 * @details Sets a gauge field u to unit 3x3 matrix on the local lattice.
 * @param u pointer to an initialized gauge field structure
 * @warning The outer boundaries are not set to 1.
 * @todo Maybe (?) set also the outer boundaries to 1
 * @todo add initialization check
 */
void qcd_unitGaugeField(qcd_gaugeField *u);

/**
 * @fn void qcd_setVector(qcd_vector *vec, qcd_complex_16 c)
 * @brief Sets a spinor field to a (complex) constant value.
 * @details The spinor field vec is set to the value c on the local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param c a complex number
 * @warning The outer boundaries are not set to c.
 * @todo Maybe (?) set also the outer boundaries to c
 * @todo add initialization check
 */
void qcd_setVector(qcd_vector *vec, qcd_complex_16 c);

/**
 * @fn void qcd_setVector3d(qcd_vector3d *vec, qcd_complex_16 c)
 * @brief Sets a 3d spinor field to a (complex) constant value.
 * @details The time-slice spinor field vec is set to the value c on
 *          the local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector3d structure
 * @param[in] c a complex number
 * @warning The outer boundaries are not set to c.
 * @todo Maybe (?) set also the outer boundaries to c
 * @todo add initialization check
 */
void qcd_setVector3d(qcd_vector3d *vec, qcd_complex_16 c);

/**
 * @fn void qcd_setSpinorComponent3d(qcd_spinorComponent3d *sp, qcd_complex_16 c)
 * @brief Sets the 3d spinor-component field to a (complex) constant value.
 * @details The time-slice spionor-component field sp is set to the value c on
 *          the local lattice.
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] c a complex number
 * @warning The outer boundaries are not set to c.
 * @todo Maybe (?) set also the outer boundaries to c
 * @todo add initialization check
 */
void qcd_setSpinorComponent3d(qcd_spinorComponent3d *sp, qcd_complex_16 c);

/**
 * @fn void qcd_copyVector(qcd_vector *dest, qcd_vector *src)
 * @brief Copies a spinor field.
 * @details The local lattice data in the spinor field src is copied to 
 *          the spinor field dest.
 * @param[in,out] dest pointer to an initialized qcd_vector structure
 * @param[in] src pointer to an initialized qcd_vector structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVector(qcd_vector *dest, qcd_vector *src);

/**
 * @fn void qcd_copyVector3d(qcd_vector3d *dest, qcd_vector3d *src)
 * @brief Copies a 3d spinor field.
 * @details The local lattice data in the time-slice spinor field src are 
 *          copied to the time-slice spinor field dest.
 * @param[in,out] dest pointer to an initialized qcd_vector3d structure
 * @param[in] src pointer to an initialized qcd_vector3d structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVector3d(qcd_vector3d *dest, qcd_vector3d *src);

/**
 * @fn void qcd_copySpinorComponent3d(qcd_spinorComponent3d *dest, qcd_spinorComponent3d *src)
 * @brief Copies a 3d spinor-component field.
 * @details The local lattice data in the time-slice spinor-component field src
 *          are copied to the time-slice spinor-component field dest.
 * @param[in,out] dest pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] src dest pointer to an initialized qcd_spinorComponent3d structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copySpinorComponent3d(qcd_spinorComponent3d *dest, qcd_spinorComponent3d *src);

/**
 * @fn void qcd_copyGaugeField(qcd_gaugeField *dest, qcd_gaugeField *src)
 * @brief Copies a gauge field.
 * @details The local lattice data of the gauge field src is copied to the 
 *          gauge field dest.
 * @param[in,out] dest pointer to an initialized qcd_gaugeField structure
 * @param[in] src pointer to an initialized qcd_gaugeField structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyGaugeField(qcd_gaugeField *dest, qcd_gaugeField *src);

/**
 * @fn void qcd_copyVectorPropagator(qcd_vector *vec, qcd_propagator *prop, int nu, int c2)
 * @brief Copies a dirac/color component of a propagator to a spinor field.
 * @details The local lattice data of the 
 *          nu (dirac) / c2 (color) component of the propagator prop is copied 
 *          to the spinor field vec.\n
 *          vec->D[i][mu][c1] = prop->D[i][mu][nu][c1][c2] for all i,mu,c1
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] prop pointer to an initialized qcd_propagator structure
 * @param[in] nu Dirac component 0..3 of the propagator which to copy
 * @param[in] c2 color component 0..2 of the propagator which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVectorPropagator(qcd_vector *vec, qcd_propagator *prop, int nu, int c2);

/**
 * @fn void qcd_copyPropagatorVector(qcd_propagator *prop, qcd_vector *vec, int nu, int c2)
 * @brief Copies a spinor field to a dirac/color component of a propagator.
 * @details The local lattice data of the spinor field vec is copied to the 
 *          nu (dirac) / c2 (color) component of the propagator prop.\n
 *          prop->D[i][mu][nu][c1][c2] = vec->D[i][mu][c1]   for all i,mu,c1
 * @param[in,out] prop pointer to an initialized qcd_propagator structure
 * @param[in] vec pointer to an initialized qcd_vector structure
 * @param[in] nu Dirac component 0..3 of the propagator to which to copy
 * @param[in] c2 color component 0..2 of the propagator to which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyPropagatorVector(qcd_propagator *prop, qcd_vector *vec, int nu, int c2);

/**
 * @fn void qcd_copyVectorPropagator2(qcd_vector *vec, qcd_propagator *prop, int mu, int c1)
 * @brief Copies a dirac/color component of a propagator to a spinor field.
 * @details The local lattice data of the 
 *          mu (dirac) / c1 (color) component of the propagator prop is copied 
 *          to the spinor field vec.\n
 *          vec->D[i][nu][c2] = prop->D[i][mu][nu][c1][c2] for all i,nu,c2
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] prop pointer to an initialized qcd_propagator structure
 * @param[in] mu Dirac component 0..3 of the propagator which to copy
 * @param[in] c1 color component 0..2 of the propagator which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVectorPropagator2(qcd_vector *vec, qcd_propagator *prop, int mu, int c1);

/**
 * @fn void qcd_copyPropagatorVector2(qcd_propagator *prop, qcd_vector *vec, int mu, int c1)
 * @brief Copies a spinor field to a dirac/color component of a propagator.
 * @details The local lattice data of the spinor field vec is copied to the 
 *          mu (dirac) / c1 (color) component of the propagator prop.\n
 *          prop->D[i][mu][nu][c1][c2] = vec->D[i][nu][c2]   for all i,nu,c2
 * @param[in,out] prop pointer to an initialized qcd_propagator structure
 * @param[in] vec pointer to an initialized qcd_vector structure
 * @param[in] mu Dirac component 0..3 of the propagator to which to copy
 * @param[in] c1 color component 0..2 of the propagator to which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyPropagatorVector2(qcd_propagator *prop, qcd_vector *vec, int mu, int c1);

/**
 * @fn void qcd_copyVector3dVector(qcd_vector3d *vec3d, qcd_vector *vec, int t)
 * @brief Copy a timeslice of a spinor field into a 3d spinor field.
 * @details Time slice t of spinor field vec is copied into the 3d spinor field 
 *          vec3d. All processes have to call the function simultaneously.
 * @param[in,out] vec3d pointer to an initialized qcd_vector3d structure
 * @param[in] vec pointer to an initialized qcd_vector structure
 * @param t time-slice (on the global lattice) which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVector3dVector(qcd_vector3d *vec3d, qcd_vector *vec, int t);

/**
 * @fn void qcd_copyVectorVector3d(qcd_vector *vec, qcd_vector3d *vec3d, int t)
 * @brief Copy a 3d spinor field into a timeslice of a spinor field.
 * @details The 3d spinor field vec3d is copied into time slice t of spinor 
 *          field vec.
 *          All processes have to call the function simultaneously.
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] vec3d pointer to an initialized qcd_vector3d structure
 * @param t time-slice (on the global lattice) which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVectorVector3d(qcd_vector *vec, qcd_vector3d *vec3d, int t);

/**
 * @fn void qcd_copySpinorComponent3dVector(qcd_spinorComponent3d *sp3d, qcd_vector *vec, int t, int mu);
 * @brief Copies a component of a time slice of a spinor field into a 3d spinor-component field.
 * @details The mu (Dirac) component of the time-slice t of the spinor field 
 *          vec is copied into the 3d spinor-component field sp3d. All processes
 *          have to call this function simultaneously.
 * @param[in,out] sp3d pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] vec pointer to a qcd_vector structure
 * @param[in] t time-slice (on the global lattice) which to copy
 * @param[in] mu Dirac component of the spinor which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copySpinorComponent3dVector(qcd_spinorComponent3d *sp3d, qcd_vector *vec, int t, int mu);

/**
 * @fn void qcd_copyVectorSpinorComponent3d(qcd_vector *vec, qcd_spinorComponent3d *sp3d, int t, int mu);
 * @brief Copies a 3d spinor component field to  a component of a time slice of a spinor field.
 * @details The 3d spinor-component field sp3d is copied into the mu (Dirac) 
 *          component of the time-slice t of the spinor field 
 *          vec. All processes
 *          have to call this function simultaneously.
 * @param[in,out] vec pointer to a qcd_vector structure
 * @param[in] sp3d pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] t time-slice (on the global lattice) which to copy
 * @param[in] mu Dirac component of the spinor which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyVectorSpinorComponent3d(qcd_vector *vec, qcd_spinorComponent3d *sp3d, int t, int mu);

/**
 * @fn void qcd_copyGaugeField3dGaugeField(qcd_gaugeField3d *u3d, qcd_gaugeField *u, int t)
 * @brief Copies a time slice of a gauge field into a 3d gauge field.
 * @details Time slice t of the gauge field u is copied to the 3d gauge field
 *          u3d. All processes have to call this function simultaneously
 * @param[in,out] u3d pointer to an initialized qcd_gaugeField3d structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @param[in] t time-slice (on the global lattice) which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyGaugeField3dGaugeField(qcd_gaugeField3d *u3d, qcd_gaugeField *u, int t);

/**
 * @fn void qcd_copyGaugeFieldGaugeField3d(qcd_gaugeField *u, qcd_gaugeField3d *u3d, int t)
 * @brief Copies a 3d gauge Field into a time slice of a gauge field.
 * @details The 3d gauge field u3d is copied to time slice t of the gauge field u.
 *          All processes have to call this function simultaneously
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @param[in] u3d pointer to an initialized qcd_gaugeField3d structure
 * @param[in] t time-slice (on the global lattice) to which to copy
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyGaugeFieldGaugeField3d(qcd_gaugeField *u, qcd_gaugeField3d *u3d, int t);

/**
 * @fn void qcd_copyGaugeFieldGaugeFieldO(qcd_gaugeField *u, qcd_gaugeFieldO *uO)
 * @brief Copies a gauge field in openQCD layout to the default layout.
 * @details The gauge field uO that is defined on the 8 links emerging from odd
 *          sites (openQCD layout), is copied to the standard gauge field u.
 *          All processes have to call this function.
 * @param[in,out] u pointer to an initialized qcd_gaugeFiel structure
 * @param[in] uO pointer to an initialized qcd_gaugeFieldO structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyGaugeFieldGaugeFieldO(qcd_gaugeField *u, qcd_gaugeFieldO *uO);

/**
 * @fn void qcd_copyGaugeFieldOGaugeField(qcd_gaugeFieldO *uO, qcd_gaugeField *u)
 * @brief Copies a gauge field to a gauge field in openQCD layout.
 * @details The (standard) gauge field u is copied to uO which is defined on 
 *          the 8 links emerging from odd
 *          sites (openQCD layout).
 *          All processes have to call this function.
 * @param[in,out] uO pointer to an initialized qcd_gaugeFielO structure
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @warning The outer boundaries are not copied.
 * @todo Copy also the outer boundaries
 * @todo add initialization check
 */
void qcd_copyGaugeFieldOGaugeField(qcd_gaugeFieldO *uO, qcd_gaugeField *u);

/**
 * @fn void qcd_copyCloverField(qcd_cloverField *dest, qcd_cloverField *src)
 * @brief Copies a clover field.
 * @details The clover field src is copied to dest.
 * @param[in,out] dest pointer to an initialized qcd_cloverField structure
 * @param[in] src pointer to an initialized qcd_cloverField structure
 * @todo add initialization check
 */
void qcd_copyCloverField(qcd_cloverField *dest, qcd_cloverField *src);

/**
 * @fn void qcd_mulVectorC(qcd_vector *vec, qcd_complex_16 c)
 * @brief Multiplies a spinor field by a complex number
 * @details The spinor field vec gets multiplied by the complex number c. This
 *          is done on the whole local lattice.\n
 *          vec->D[x][dirac][color] = vec->D[x][dirac][color] * c
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] c complex number 
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 */
void qcd_mulVectorC(qcd_vector *vec, qcd_complex_16 c);

/**
 * @fn void qcd_mulVector3dC(qcd_vector3d *vec, qcd_complex_16 c)
 * @brief Multiplies a 3d spinor field by a complex number
 * @details The time-slice spinor field vec gets 
 *          multiplied by the complex number c. This
 *          is done on the whole local lattice.\n
 *          vec->D[x][dirac][color] = vec->D[x][dirac][color] * c
 * @param[in,out] vec pointer to an initialized qcd_vector3d structure
 * @param[in] c complex number 
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_mulVector3dC(qcd_vector3d *vec, qcd_complex_16 c);

/**
 * @fn void qcd_mulSpinorComponent3dC(qcd_spinorComponent3d *sp, qcd_complex_16 c)
 * @brief Multiplies a 3d spinor component field by a complex number
 * @details The time-slice spinor component field sp gets 
 *          multiplied by the complex number c. This
 *          is done on the whole local lattice.\n
 *          sp->D[x][color] = sp->D[x][color] * c
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] c complex number 
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_mulSpinorComponent3dC(qcd_spinorComponent3d *sp, qcd_complex_16 c);

/**
 * @fn void qcd_mulPropagatorC(qcd_propagator *prop, qcd_complex_16 c)
 * @brief Multiplies a propagator by a complex number
 * @details The propagator prop gets 
 *          multiplied by the complex number c. This
 *          is done on the whole local lattice.\n
 *          prop->D[x][d1][d2][c1][c2] = prop->D[x][d1][d2][c1][c2] * c
 * @param[in,out] prop pointer to an initialized qcd_propagator structure
 * @param[in] c complex number 
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_mulPropagatorC(qcd_propagator *prop, qcd_complex_16 c);

/**
 * @fn void qcd_mulSliceOfPropagatorC(qcd_propagator *prop, qcd_complex_16 c, int t)
 * @brief Multiplies a time-slice of a propagator by a complex number
 * @details The time-slice t (global lattice) of the propagator prop gets 
 *          multiplied by the complex number c. If t lies in the local lattice,
 *          the whole local time slice is multiplied.
 * @param[in,out] prop pointer to an initialized qcd_propagator structure
 * @param[in] c complex number 
 * @param[in] t time-slice (global lattice)
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_mulSliceOfPropagatorC(qcd_propagator *prop, qcd_complex_16 c, int t);

/**
 * @fn void qcd_mulSliceOfVectorC(qcd_vector *vec, qcd_complex_16 c, int t)
 * @brief Multiplies a time-slice of a spinor field by a complex number
 * @details The time-slice t (global lattice) of the spinor field vec gets
 *          multiplied by the complex number c. If t lies in the local lattice,
 *          the whole local time slice is multiplied.
 * @param[in,out] prop pointer to an initialized qcd_vector structure
 * @param[in] c complex number
 * @param[in] t time-slice (global lattice)
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_mulSliceOfVectorC(qcd_vector *vec, qcd_complex_16 c, int t);

/**
 * @fn void qcd_scaleVector(qcd_vector *vec, double r)
 * @brief Multiplies a spinor field by a real number.
 * @details The spinor field vec gets multiplied by the real number r on the 
 *          whole local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] r real rescaling factor
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_scaleVector(qcd_vector *vec, double r);

/**
 * @fn void qcd_scaleVector3d(qcd_vector3d *vec, double r)
 * @brief Multiplies a 3d spinor field by a real number.
 * @details The time-slice spinor field vec gets multiplied by the real number
 *          r on the whole local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector3d structure
 * @param[in] r real rescaling factor
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_scaleVector3d(qcd_vector3d *vec, double r);

/**
 * @fn void qcd_scaleSpinorComponent3d(qcd_spinorComponent3d *sp, double r)
 * @brief Multiplies a 3d spinor-component field by a real number
 * @details The time-slice spinor-component field sp gets multiplied by the real number r
 *          on the local lattice.
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @param r real rescaling factor
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_scaleSpinorComponent3d(qcd_spinorComponent3d *sp, double r);

/**
 * @fn void qcd_scaleSliceOfVector(qcd_vector *vec, double r, int t)
 * @brief Multiply a time-slice of a spinor field by a real number.
 * @details Time slice t (global lattice) of the spinor field vec gets multiplied
 *          by the real number r on the whole local lattice.
 * @param[in,out] vec pointer to an initialized qcd_vector structure
 * @param[in] r real rescaling factor
 * @param[in] t time slice (global lattice) which to rescale
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_scaleSliceOfVector(qcd_vector *vec, double r, int t);

/**
 * @fn void qcd_scaleGaugeField(qcd_gaugeField *u, double alpha)
 * @brief Multiplies a gauge field by a real number
 * @details The gaugefield u gets multiplied by the real number alpha on the 
 *          whole local lattice. After this operation u won't be a SU(3) matrix
 *          anymore.
 * @param[in,out] u pointer to an initialized gaugeField structure
 * @param[in] alpha real rescaling factor
 * @warning The outer boundary points are not multiplied
 * @todo maybe (?) multiply also outer boundary points
 * @todo add initialization check
 */
void qcd_scaleGaugeField(qcd_gaugeField *u, double alpha);

/**
 * @fn void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend)
 * @brief Subtracts two spinor fields from each other.
 * @details Spinor field subtrahend gets subtracted from spinor field minuend.
 *          The result is stored in spinor field difference. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] difference pointer to an initialized qcd_vector structure
 * @param[in] minuend pointer to an initialized qcd_vector structure
 * @param[in] subtrahend pointer to an initialized qcd_vector structure
 * @todo add initialization check
 */
void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend);

/**
 * @fn void qcd_subPropagator(qcd_propagator *difference, qcd_propagator *minuend, qcd_propagator subtrahend)
 * @brief Subtracts two propagator fields from each other.
 * @details Proppagator subtrahend gets subtracted from propagator minuend.
 *          The result is stored in propagator difference. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] difference pointer to an initialized qcd_propagator structure
 * @param[in] minuend pointer to an initialized qcd_propagator structure
 * @param[in] subtrahend pointer to an initialized qcd_propagator structure
 * @todo add initialization check
 */
void qcd_subPropagator(qcd_propagator *difference, qcd_propagator *minuend, qcd_propagator *subtrahend);

/**
 * @fn void qcd_subVector3d(qcd_vector3d *difference, qcd_vector3d *minuend, qcd_vector3d *subtrahend)
 * @brief Subtracts two 3d spinor fields from each other.
 * @details Time-slice spinor field subtrahend gets subtracted from time-slice 
 *          spinor field minuend.
 *          The result is stored in time-slice spinor field difference. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] difference pointer to an initialized qcd_vector3d structure
 * @param[in] minuend pointer to an initialized qcd_vector3d structure
 * @param[in] subtrahend pointer to an initialized qcd_vector3d structure
 * @warning The outer boundary points are not subtracted
 * @todo maybe (?) subtract also outer boundary points
 * @todo add initialization check
 */
void qcd_subVector3d(qcd_vector3d *difference, qcd_vector3d *minuend, qcd_vector3d *subtrahend);

/**
 * @fn void qcd_subSpinorComponent3d(qcd_spinorComponent3d *difference, qcd_spinorComponent3d *minuend, qcd_spinorComponent3d *subtrahend)
 * @brief Subtracts two 3d spinor-component fields from each other.
 * @details Time-slice spinor-component field subtrahend gets subtracted from time-slice 
 *          spinor-component field minuend.
 *          The result is stored in time-slice spinor-component field difference. 
 *          The operation is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] difference pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] minuend pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] subtrahend pointer to an initialized qcd_spinorComponent3d structure
 * @warning The outer boundary points are not subtracted
 * @todo maybe (?) subtract also outer boundary points
 * @todo add initialization check
 */
void qcd_subSpinorComponent3d(qcd_spinorComponent3d *difference, qcd_spinorComponent3d *minuend, qcd_spinorComponent3d *subtrahend);

/**
 * @fn void qcd_subGaugeField(qcd_gaugeField *difference, qcd_gaugeField *minuend, qcd_gaugeField *subtrahend);
 * @brief Subtracts two gauge fields from each other.
 * @details Gauge field subtrahend gets subtracted from gauge field minuend.
 *          The result is stored in gauge field difference. 
 *          The operation is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide. The gauge field difference, is not a SU(3) 
 *          matrix anymore.
 * @param[in,out] difference pointer to an initialized qcd_gaugeField structure
 * @param[in] minuend pointer to an initialized qcd_gaugeField structure
 * @param[in] subtrahend pointer to an initialized qcd_gaugeField structure
 * @warning The outer boundary points are not subtracted
 * @todo maybe (?) subtract also outer boundary points
 * @todo add initialization check
 * @todo optimization: single loop
 */
void qcd_subGaugeField(qcd_gaugeField *difference, qcd_gaugeField *minuend, qcd_gaugeField *subtrahend);

/**
 * @fn void qcd_subCloverField(qcd_cloverField *difference, qcd_cloverField *minuend, qcd_cloverField *subtrahend)
 * @brief Subtracts two clover fields from each other.
 * @details Clover field subtrahend gets subtracted from clover field minuend.
 *          The result is stored in clover field difference. 
 *          The operation is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] difference pointer to an initialized qcd_cloverField structure
 * @param[in] minuend pointer to an initialized qcd_cloverField structure
 * @param[in] subtrahend pointer to an initialized qcd_cloverField structure
 * @todo add initialization check
 * @todo optimization: single loop
 */
void qcd_subCloverField(qcd_cloverField *difference, qcd_cloverField *minuend, qcd_cloverField *subtrahend);

/**
 * @fn void qcd_addVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2);
 * @brief Adds two spinor fields to each other.
 * @details Spinor field summand1 is added to spinor field summand2.
 *          The result is stored in spinor field sum. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] sum pointer to an initialized qcd_vector structure
 * @param[in] summand1 pointer to an initialized qcd_vector structure
 * @param[in] summand2 pointer to an initialized qcd_vector structure
 * @warning The outer boundary points are not summed
 * @todo maybe (?) sum also outer boundary points
 * @todo add initialization check
 */
void qcd_addVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2);

/**
 * @fn void qcd_addVector3d(qcd_vector3d *sum, qcd_vector3d *summand1, qcd_vector3d *summand2);
 * @brief Adds two 3d spinor fields to each other.
 * @details Time-slice spinor field summand1 is added to time-slice spinor field summand2.
 *          The result is stored in time-slice spinor field sum. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] sum pointer to an initialized qcd_vector3d structure
 * @param[in] summand1 pointer to an initialized qcd_vector3d structure
 * @param[in] summand2 pointer to an initialized qcd_vector3d structure
 * @warning The outer boundary points are not summed
 * @todo maybe (?) sum also outer boundary points
 * @todo add initialization check
 */
void qcd_addVector3d(qcd_vector3d *sum, qcd_vector3d *summand1, qcd_vector3d *summand2);

/**
 * @fn void qcd_addSpinorComponent3d(qcd_spinorComponent3d *sum, qcd_spinorComponent3d *summand1, qcd_spinorComponent3d *summand2)
 * @brief Adds two 3d spinor component fields to each other.
 * @details Time-slice spinor-component field summand1 is added 
 *          to time-slice spinor-component field summand2.
 *          The result is stored in time-slice spinor-component field sum. The operation
 *          is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide.
 * @param[in,out] sum pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] summand1 pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] summand2 pointer to an initialized qcd_spinorComponent3d structure
 * @warning The outer boundary points are not summed
 * @todo maybe (?) sum also outer boundary points
 * @todo add initialization check
 */
void qcd_addSpinorComponent3d(qcd_spinorComponent3d *sum, qcd_spinorComponent3d *summand1, qcd_spinorComponent3d *summand2);

/**
 * @fn void qcd_addSliceOfVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2, int t)
 * @brief Adds two time-slices of a spinor field to each other.
 * @details Time slice t (global lattice) of the spinor fields  summand1 and
 *          summand2 are summed and the result is stored in the corresponding
 *          time-slice of spinor field sum. The operation is carried out on the
 *          whole local lattice, if t falls into it.
 * @param[in,out] sum pointer to an initialized qcd_vector structure
 * @param[in] summand1 pointer to an initialized qcd_vector structure
 * @param[in] summand2 pointer to an initialized qcd_vector structure
 * @param[in] t time-slice (global lattice)
 * @warning The outer boundary points are not summed
 * @todo maybe (?) sum also outer boundary points
 * @todo add initialization check
 */
void qcd_addSliceOfVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2, int t);

/**
 * @fn void qcd_addGaugeField(qcd_gaugeField *sum, qcd_gaugeField *summand1, qcd_gaugeField *summand2)
 * @brief Adds two gauge fields to each other.
 * @details Gauge field summand1 is added to gauge field summand2.
 *          The result is stored in gauge field sum. 
 *          The operation is performed on the whole local lattice. Two or three of the 
 *          arguments may coincide. The gauge field sum, is not a SU(3) 
 *          matrix anymore.
 * @param[in,out] sum pointer to an initialized qcd_gaugeField structure
 * @param[in] summand1 pointer to an initialized qcd_gaugeField structure
 * @param[in] summand2 pointer to an initialized qcd_gaugeField structure
 * @warning The outer boundary points are not summed
 * @todo maybe (?) sum also outer boundary points
 * @todo add initialization check
 * @todo optimization: single loop
 */
void qcd_addGaugeField(qcd_gaugeField *sum, qcd_gaugeField *summand1, qcd_gaugeField *summand2);

/**
 * @fn void qcd_axpyVector(qcd_vector *axpy, qcd_complex_16 a, qcd_vector *x, qcd_vector *y)
 * @brief Adds a rescaled spinor field to another.
 * @details The spinor field x is multiplied by the complex number a and added to
 *          the spinor field y. The result is stored in the spinor field axpy.\n
 *          axpy->D[x][d][c] =  a*x->D[x][d][c] + y->D[x][d][c]\n
 *          The operation is performed on the whole local lattice. Spinor fields
 *          x and y may be identical (but not the result).
 * @param[in,out] axpy pointer to an initialized qcd_vector structure
 * @param[in] a complex factor
 * @param[in] x pointer to an initialized qcd_vector structure
 * @param[in] y pointer to an initialized qcd_vector structure
 * @warning The outer boundary points are not axpy'd
 * @todo maybe (?) axpy also outer boundary points
 * @todo add initialization check
 */
void qcd_axpyVector(qcd_vector *axpy, qcd_complex_16 a, qcd_vector *x, qcd_vector *y);

/**
 * @fn void qcd_axpyVector3d(qcd_vector3d *axpy, qcd_complex_16 a, qcd_vector3d *x, qcd_vector3d *y)
 * @brief Adds a rescaled 3d spinor field to another.
 * @details The time-slice spinor field x is multiplied by the complex number a and added to
 *          the time-slice spinor field y. The result is stored in the time-slice 
 *          spinor field axpy.\n
 *          axpy->D[x][d][c] =  a*x->D[x][d][c] + y->D[x][d][c]\n
 *          The operation is performed on the whole local lattice. Spinor fields
 *          x and y may be identical (but not the result).
 * @param[in,out] axpy pointer to an initialized qcd_vector3d structure
 * @param[in] a complex factor
 * @param[in] x pointer to an initialized qcd_vector3d structure
 * @param[in] y pointer to an initialized qcd_vector3d structure
 * @warning The outer boundary points are not axpy'd
 * @todo maybe (?) axpy also outer boundary points
 * @todo add initialization check
 */
void qcd_axpyVector3d(qcd_vector3d *axpy, qcd_complex_16 a, qcd_vector3d *x, qcd_vector3d *y);

/**
 * @fn void qcd_axpySpinorComponent3d(qcd_spinorComponent3d *axpy, qcd_complex_16 a, qcd_spinorComponent3d *x, qcd_spinorComponent3d *y)
 * @brief Adds a rescaled 3d spinor-component field to another.
 * @details The time-slice spinor-component field x is 
 *          multiplied by the complex number a and added to
 *          the time-slice spinor-component field y. The result is stored in 
 *          the time-slice spinor-component field axpy.\n
 *          axpy->D[x][c] =  a*x->D[x][c] + y->D[x][c]\n
 *          The operation is performed on the whole local lattice. Spinor-component fields
 *          x and y may be identical (but not the result).
 * @param[in,out] axpy pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] a complex factor
 * @param[in] x pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] y pointer to an initialized qcd_spinorComponent3d structure
 * @warning The outer boundary points are not axpy'd
 * @todo maybe (?) axpy also outer boundary points
 * @todo add initialization check
 */
void qcd_axpySpinorComponent3d(qcd_spinorComponent3d *axpy, qcd_complex_16 a, qcd_spinorComponent3d *x, qcd_spinorComponent3d *y);

/**
 * @fn void qcd_raxpyVector3d(qcd_vector3d *axpy, double a, qcd_vector3d *x, qcd_vector3d *y)
 * @brief Adds a rescaled 3d spinor field to another.
 * @details The time-slice spinor field x is multiplied by the real number a and added to
 *          the time-slice spinor field y. The result is stored in the time-slice 
 *          spinor field axpy.\n
 *          axpy->D[x][d][c] =  a*x->D[x][d][c] + y->D[x][d][c]\n
 *          The operation is performed on the whole local lattice. Spinor fields
 *          x and y and axpy may be identical.
 * @param[in,out] axpy pointer to an initialized qcd_vector3d structure
 * @param[in] a real factor
 * @param[in] x pointer to an initialized qcd_vector3d structure
 * @param[in] y pointer to an initialized qcd_vector3d structure
 * @warning The outer boundary points are not axpy'd
 * @todo maybe (?) axpy also outer boundary points
 * @todo add initialization check
 */
void qcd_raxpyVector3d(qcd_vector3d *axpy, double a, qcd_vector3d *x, qcd_vector3d *y);

/**
 * @fn void qcd_raxpySpinorComponent3d(qcd_spinorComponent3d *axpy, double a, qcd_spinorComponent3d *x, qcd_spinorComponent3d *y)
 * @brief Adds a rescaled 3d spinor-component field to another.
 * @details The time-slice spinor-component field x is 
 *          multiplied by the real number a and added to
 *          the time-slice spinor-component field y. The result is stored in 
 *          the time-slice spinor-component field axpy.\n
 *          axpy->D[x][c] =  a*x->D[x][c] + y->D[x][c]\n
 *          The operation is performed on the whole local lattice. Spinor-component fields
 *          x and y and axpy may be identical.
 * @param[in,out] axpy pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] a real factor
 * @param[in] x pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] y pointer to an initialized qcd_spinorComponent3d structure
 * @warning The outer boundary points are not axpy'd
 * @todo maybe (?) axpy also outer boundary points
 * @todo add initialization check
 */
void qcd_raxpySpinorComponent3d(qcd_spinorComponent3d *axpy, double a, qcd_spinorComponent3d *x, qcd_spinorComponent3d *y);

/**
 * @fn qcd_complex_16 qcd_mulAdjointVectorVector(qcd_vector *factor1, qcd_vector *factor2);
 * @brief Complex inner product of two spinor fields.
 * @details The inner product of spinor field factor1 and spinor field factor 2 is 
 *          computed and returned (complex).\n
 *          \f$ {\tt result} = \sum\limits_{x,d,c} {\rm conj}({\tt factor1->D[x][d][c]})\, {\tt factor2->D[x][d][c]}\f$ \n
 *          All processes must call the function simultaneously. Spinor fields
 *          x and y may be identical.
 * @param[in] factor1 pointer to an initialized qcd_vector structure
 * @param[in] factor2 pointer to an initialized qcd_vector structure
 * @return the (complex) inner product
 */
qcd_complex_16 qcd_mulAdjointVectorVector(qcd_vector *factor1, qcd_vector *factor2);

/**
 * @fn qcd_complex_16 qcd_mulAdjointVector3dVector3d(qcd_vector3d *factor1, qcd_vector3d *factor2);
 * @brief Complex inner product of two 3d spinor fields.
 * @details The inner product of time-slice spinor field factor1 and time-slice 
 *          spinor field factor 2 is 
 *          computed and returned (complex).\n
 *          \f$ {\tt result} = \sum\limits_{x,d,c} {\rm conj}({\tt factor1->D[x][d][c]})\, {\tt factor2->D[x][d][c]}\f$ \n
 *          All processes must call the function simultaneously. Spinor fields
 *          x and y may be identical.
 * @param[in] factor1 pointer to an initialized qcd_vector3d structure
 * @param[in] factor2 pointer to an initialized qcd_vector3d structure
 * @return the (complex) inner product
 */
qcd_complex_16 qcd_mulAdjointVector3dVector3d(qcd_vector3d *factor1, qcd_vector3d *factor2);

/**
 * @fn qcd_complex_16 qcd_mulAdjointSpinorComponent3dSpinorComponent3d(qcd_spinorComponent3d *factor1, qcd_spinorComponent3d *factor2)
 * @brief Complex inner product of two 3d spinor-component fields.
 * @details The inner product of time-slice spinor-component field factor1 
 *          and time-slice spinor-component field factor 2 is 
 *          computed and returned (complex).\n
 *          \f$ {\tt result} = \sum\limits_{x,c} {\rm conj}({\tt factor1->D[x][c]})\, {\tt factor2->D[x][c]}\f$ \n
 *          All processes must call the function simultaneously. Spinor-component
 *          fields x and y may be identical.
 * @param[in] factor1 pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] factor2 pointer to an initialized qcd_spinorComponent3d structure
 * @return the (complex) inner product
 */
qcd_complex_16 qcd_mulAdjointSpinorComponent3dSpinorComponent3d(qcd_spinorComponent3d *factor1, qcd_spinorComponent3d *factor2);

/**
 * @fn double qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(qcd_spinorComponent3d *factor1, qcd_spinorComponent3d *factor2);
 * @brief Real part of the inner product of two 3d spinor-component fields.
 * @details The inner product of time-slice spinor-component field factor1 
 *          and time-slice spinor-component field factor 2 is 
 *          computed and its real part returned.\n
 *          \f$ {\tt result} = \sum\limits_{x,c} {\rm Re}( {\rm conj}({\tt factor1->D[x][c]})\, {\tt factor2->D[x][c]} )\f$ \n
 *          All processes must call the function simultaneously. Spinor-component
 *          fields x and y may be identical.
 * @param[in] factor1 pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] factor2 pointer to an initialized qcd_spinorComponent3d structure
 * @return the real part of the inner product
 */
double qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(qcd_spinorComponent3d *factor1, qcd_spinorComponent3d *factor2);

/**
 * @fn double qcd_normVector(qcd_vector *vec)
 * @brief Computes the 2-norm of a spinor field.
 * @details The 2-norm of spinor-field vec is computed and returned.\n
 *          \f$ {\tt result} = \sqrt{ \sum\limits_{x,d,c} {\rm conj}({\tt vec->D[x][d][c]})\, {\tt vec->D[x][d][c]} }\f$ \n
 *          All processes must call this function.
 * @param[in] vec pointer to an initialized qcd_vector structure
 * @return The (real, positive) 2-norm
 */
double qcd_normVector(qcd_vector *vec);

/**
 * @fn double qcd_normPropagator(qcd_propagator *prop)
 * @brief Computes the 2-norm of a propagator.
 * @details The 2-norm of propagator prop is computed and returned.\n
 *          \f$ {\tt result} = \sqrt{ \sum\limits_{x,d1,d2,c1,c2} {\rm conj}({\tt prop->D[x][d1][d2][c1][c2]})\, {\tt prop->D[x][d1][d2][c1][c2]} }\f$ \n
 *          All processes must call this function.
 * @param[in] prop pointer to an initialized qcd_propagator structure
 * @return The (real, positive) 2-norm
 */
double qcd_normPropagator(qcd_propagator *prop);

/**
 * @fn double qcd_normVector3d(qcd_vector3d *vec)
 * @brief Computes the 2-norm of a 3d spinor field.
 * @details The 2-norm of time-slice spinor-field vec is computed and returned.\n
 *          \f$ {\tt result} = \sqrt{ \sum\limits_{x,d,c} {\rm conj}({\tt vec->D[x][d][c]})\, {\tt vec->D[x][d][c]} }\f$ \n
 *          All processes must call this function.
 * @param[in] vec pointer to an initialized qcd_vector3d structure
 * @return The (real, positive) 2-norm
 */
double qcd_normVector3d(qcd_vector3d *vec);

/**
 * @fn double qcd_normSpinorComponent3d(qcd_spinorComponent3d *sp)
 * @brief Computes the 2-norm of a 3d spinor-component field.
 * @details The 2-norm of time-slice spinor-component field sp is 
 *          computed and returned.\n
 *          \f$ {\tt result} = \sqrt{ \sum\limits_{x,c} {\rm conj}({\tt sp->D[x][c]}) ({\tt sp->D[x][c]}) }\f$ \n
 *          All processes must call this function.
 * @param[in] sp pointer to an initialized qcd_vector3d structure
 * @return The (real, positive) 2-norm
 */
double qcd_normSpinorComponent3d(qcd_spinorComponent3d *sp);

/**
 * @fn double qcd_normGaugeField(qcd_gaugeField *u)
 * @brief Compute the norm of a gauge field.
 * @details The norms of gauge fields on all links is 
 *          summed and returned. As norm on a single link,
 *          the Frobenius norm is used. 
 *          All processes must call this function simultaneously.
 * @param[in] u pointer to an initialized qcd_gaugeField structure
 * @return a (real, positive) norm
 */
double qcd_normGaugeField(qcd_gaugeField *u);

/**
 * @fn double qcd_normGaugeFieldO(qcd_gaugeFieldO *u)
 * @brief Compute the norm of a gauge field.
 * @details The norms of gauge fields on all links is 
 *          summed and returned. As norm on a single link,
 *          the Frobenius norm is used. 
 *          All processes must call this function simultaneously.
 *          This variant works on gauge fields in openQCD layout.
 * @param[in] u pointer to an initialized qcd_gaugeFieldO structure
 * @return a (real, positive) norm
 */
double qcd_normGaugeFieldO(qcd_gaugeFieldO *u);

/**
 * @fn double qcd_normCloverField(qcd_cloverField *q);
 * @brief Compute the norm of a clover field.
 * @details The norms of the clover fields on all plaquettes are
 *          summed and returned. As norm on a plaquette,
 *          the Frobenius norm of the clover field on this
 *          plaquette is taken. 
 *          All processes must call this function simultaneously.
 * @param[in] q pointer to an initialized qcd_cloverField structure
 * @return a (real, positive) norm
 */
double qcd_normCloverField(qcd_cloverField *q);

/**
 * @fn void qcd_conjPropagator(qcd_propagator *prop)
 * @brief Complex conjugation of a propagator
 * @details All local lattice components of the propagator prop are complex
 *          conjugated.
 * @param[in,out] prop pointer to an initialized propagator structure
 */
void qcd_conjPropagator(qcd_propagator *prop);

/**
 * @fn void qcd_projectSU3(qcd_gaugeField *gf)
 * @brief Projects a gauge field to SU(3).
 * @details The gaugeField gf, that contains arbitrary 3x3 matrices, is projected
 *          to the group SU(3).
 *          By default qcd_ETMC_SU3_PROJECTION is set, in which case the projection
 *          proceeds in the following way:
 *          If M is an arbitrary 3x3 matrix it can be factorized as the product
 *          of a hermitian and a unitary matrix, M = U H. We will calculate U 
 *          as the projection of M to SU(3). Now M^dag M = H^2. By
 *          going to diagonal basis for H^2 we calculate \n
 *          H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^dag \n
 *          where v is the 3x3 matrix whos columns are the eigenvectors of H^2 
 *          and L is the diagonal representation of H^2. Thus: U = M H^-1 \n
 *          [1] Phys. Lett. B307, 375-382, (1993)
 * @author Giannis Koutsou
 * @param[in,out] gf  pointer to an initialized qcd_gaugeField structure
 */
void qcd_projectSU3(qcd_gaugeField *gf);
 
/**
 * @fn void qcd_projectSU33d(qcd_gaugeField *gf)
 * @brief Projects spatial links of a gauge field to SU(3).
 * @details The gaugeField gf, that contains arbitrary 3x3 matrices, 
 *          is projected to the group SU(3) on all spatial links.
 *          By default qcd_ETMC_SU3_PROJECTION is set, in which case the projection
 *          proceeds in the following way:
 *          If M is an arbitrary 3x3 matrix it can be factorized as the product
 *          of a hermitian and a unitary matrix, M = U H. We will calculate U 
 *          as the projection of M to SU(3). Now M^dag M = H^2. By
 *          going to diagonal basis for H^2 we calculate \n
 *          H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^dag \n
 *          where v is the 3x3 matrix whos columns are the eigenvectors of H^2 
 *          and L is the diagonal representation of H^2. Thus: U = M H^-1 \n
 *          [1] Phys. Lett. B307, 375-382, (1993)
 * @author Giannis Koutsou
 * @param[in,out] gf  pointer to an initialized qcd_gaugeField3d structure
 * @warning The name suggests that this works on time-slice fields. This is 
 *          not the case, which is confusing. The name may change in the future.
 */
void qcd_projectSU33d(qcd_gaugeField *gf);


/**
 * @fn void qcd_projectSU3GaugeTransformation(qcd_gaugeTransformation *g);
 * @brief Projects a gauge transformation field to SU(3).
 * @details The gaugeTransformation g, that contains arbitrary 3x3 matrices, is projected
 *          to the group SU(3).
 * 
 *          By default qcd_ETMC_SU3_PROJECTION is set, in which case the projection
 *          proceeds in the following way:
 *          If M is an arbitrary 3x3 matrix it can be factorized as the product
 *          of a hermitian and a unitary matrix, M = U H. We will calculate U 
 *          as the projection of M to SU(3). Now M^dag M = H^2. By
 *          going to diagonal basis for H^2 we calculate \n
 *          H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^dag \n
 *          where v is the 3x3 matrix whos columns are the eigenvectors of H^2 
 *          and L is the diagonal representation of H^2. Thus: U = M H^-1 \n
 *          [1] Phys. Lett. B307, 375-382, (1993)
 * @author Giannis Koutsou, Tomasz Korzec
 * @param[in,out] g  pointer to an initialized qcd_gaugeField structure
 */
void qcd_projectSU3GaugeTransformation(qcd_gaugeTransformation *g);


#endif
