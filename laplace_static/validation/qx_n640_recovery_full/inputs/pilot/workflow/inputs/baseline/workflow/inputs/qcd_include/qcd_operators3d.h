/** @file qcd_operators3d.h
 *  @brief header file for qcd_operators3d.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2010-2015
 *  @copyright GNU Public License.
 */

#ifndef H_QCD_OPERATORS3D
#define H_QCD_OPERATORS3D 1

/* preprocessor macros */
/**
 * @def qcd_SUM_UP_HOPP
 * @brief computes t[i] = a[i] + b[i] for i=0...23
 */
#define qcd_SUM_UP_HOPP(t,a,b)  \
    { t[0]  += a[0] + b[0];  \
      t[1]  += a[1] + b[1];  \
      t[2]  += a[2] + b[2];  \
      t[3]  += a[3] + b[3];  \
      t[4]  += a[4] + b[4];  \
      t[5]  += a[5] + b[5];  \
      t[6]  += a[6] + b[6];  \
      t[7]  += a[7] + b[7];  \
      t[8]  += a[8] + b[8];  \
      t[9]  += a[9] + b[9];  \
      t[10] += a[10]+ b[10]; \
      t[11] += a[11]+ b[11]; \
      t[12] += a[12]+ b[12]; \
      t[13] += a[13]+ b[13]; \
      t[14] += a[14]+ b[14]; \
      t[15] += a[15]+ b[15]; \
      t[16] += a[16]+ b[16]; \
      t[17] += a[17]+ b[17]; \
      t[18] += a[18]+ b[18]; \
      t[19] += a[19]+ b[19]; \
      t[20] += a[20]+ b[20]; \
      t[21] += a[21]+ b[21]; \
      t[22] += a[22]+ b[22]; \
      t[23] += a[23]+ b[23]; \
   }

/**
 * @def qcd_SUM_UP_HOPP_SPINORCOMPONENT
 * @brief computes t[i] = a[i] + b[i] for i=0...5
 */
#define qcd_SUM_UP_HOPP_SPINORCOMPONENT(t,a,b)  \
    { t[0]  += a[0] + b[0];  \
      t[1]  += a[1] + b[1];  \
      t[2]  += a[2] + b[2];  \
      t[3]  += a[3] + b[3];  \
      t[4]  += a[4] + b[4];  \
      t[5]  += a[5] + b[5];  \
    }

/* prototypes */

 /**
  * @fn int qcd_DVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, int mu)
  * @brief Applies the covariant forward derivative to a 3d spinor field
  * @details Dv(x) = U_mu(x) v(x+mu) - v(x)
  * @param[in,out] Dv pointer to an initialized qcd_vector3d structure
  * @param[in] v pointer to an initialized qcd_vector3d structure
  * @param[in] u pointer to an initialized qcd_gaugeField3d structure
  * @param[in] mu direction of the derivative (1...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, int mu);

/**
  * @fn int qcd_DBackVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, int mu)
  * @brief Applies the backward covariant forward derivative to a 3d spinor field
  * @details Dv(x) = v(x) - U^dag_mu(x-mu) v(x-mu)
  * @param[in,out] Dv pointer to an initialized qcd_vector3d structure
  * @param[in] v pointer to an initialized qcd_vector3d structure
  * @param[in] u pointer to an initialized qcd_gaugeField3d structure
  * @param[in] mu direction of the derivative (1...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DBackVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, int mu);
 
  /**
  * @fn int qcd_DSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu)
  * @brief Applies the covariant forward derivative to a 3d spinor-component field
  * @details Dv(x) = U_mu(x) v(x+mu) - v(x)
  * @param[in,out] Dv pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] v pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] u pointer to an initialized qcd_gaugeField3d structure
  * @param[in] mu direction of the derivative (1...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu);

  /**
  * @fn int qcd_DBackSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu)
  * @brief Applies the backward covariant derivative to a 3d spinor-component field
  * @details Dv(x) = v(x) - U^dag_mu(x-mu) v(x-mu)
  * @param[in,out] Dv pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] v pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] u pointer to an initialized qcd_gaugeField3d structure
  * @param[in] mu direction of the derivative (1...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DBackSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu);

 
  /**
  * @fn int qcd_DSymSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu)
  * @brief Applies the symmetric covariant derivative to a 3d spinor-component field
  * @details Dv(x) = ( U_mu(x) v(x+mu) - U^dag_mu(x-mu) v(x-mu) )/2
  * @param[in,out] Dv pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] v pointer to an initialized qcd_spinorComponent3d structure
  * @param[in] u pointer to an initialized qcd_gaugeField3d structure
  * @param[in] mu direction of the derivative (1...3)
  * @return 0 on success, 1 on failure
  */
 int qcd_DSymSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, int mu);

/**
 * @fn int qcd_applyCovariantLaplaceVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, double mm)
 * @brief Applies the shifted 3d covariant Laplace operator to a spinor field.
 * @details Given the time-slice spinor field v the time-slice gauge field u and 
 *          a real shift mm,\n
 *  \f$ {\tt Dv} = [-(\nabla^*_k \nabla_k) + {\tt mm}] {\tt v} \f$\n
 *    is computed, 
 *    where \f$ \nabla_k / \nabla^*_k\f$ are covariant forward/backward
 *    derivatives in k'th direction. k is summed 1..3. All processes must 
 *    call this function simultaneously.
 * @param[in,out] Dv pointer to an initialized qcd_vector3d structure
 * @param[in] v pointer to an initialized qcd_vector3d structure
 * @param[in] u pointer to an initialized qcd_gaugeField3d structure
 * @param[in] mm real shift, can be negative
 * @return 0 on success, 1 on failure
 * @warning It is assumed, that the gauge-field communication has been performed
 *          (using #qcd_communicateGauge3dP)
 *          and the negative outer boundaries are filled with the correct entries.
 */
int qcd_applyCovariantLaplaceVector3d(qcd_vector3d *Dv, qcd_vector3d *v, qcd_gaugeField3d *u, double mm);

/**
 * @fn int qcd_applyCovariantLaplaceSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, double mm)
 * @brief Applies the shifted 3d covariant Laplace operator to a 3d spinor-component field.
 * @details Given the time-slice spinor-component field v the time-slice gauge field u and 
 *          a real shift mm,\n
 *  \f$ {\tt Dv} = [-(\nabla^*_k \nabla_k) + {\tt mm}] {\tt v} \f$\n
 *    is computed, 
 *    where \f$ \nabla_k / \nabla^*_k\f$ are covariant forward/backward
 *    derivatives in k'th direction. k is summed 1..3. All processes must 
 *    call this function simultaneously.
 * @param[in,out] Dv pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] v pointer to an initialized qcd_spinorComponent3d structure
 * @param[in] u pointer to an initialized qcd_gaugeField3d structure
 * @param[in] mm real shift, can be negative
 * @return 0 on success, 1 on failure
 * @warning It is assumed, that the gauge-field communication has been performed
 *          (using #qcd_communicateGauge3dP)
 *          and the negative outer boundaries are filled with the correct entries.
 */
int qcd_applyCovariantLaplaceSpinorComponent3d(qcd_spinorComponent3d *Dv, qcd_spinorComponent3d *v, qcd_gaugeField3d *u, double mm);
#endif
