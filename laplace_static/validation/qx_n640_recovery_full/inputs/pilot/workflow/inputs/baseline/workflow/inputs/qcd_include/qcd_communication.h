/** @file qcd_communication.h
 *  @brief header file for qcd_communication.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2009-2015
 *  @copyright GNU Public License.
 *  @warning If two or more geometries are defined, it is in principle 
 *           possible to overlap communications of fields that have
 *           different geometries. This is however at the moment not supported
 *           and leads to undefined behavior.
 */

#ifndef H_QCD_COMMUNICATION
#define H_QCD_COMMUNICATION 1

/* prototypes */

/**
 * @fn int qcd_communicateVectorPMGaugeP(qcd_vector *v, qcd_gaugeField *u)
 * @brief Sends spinor fields to +/- and gauge fields to + neighbors.
 * @details The inner boundaries of spinor field v is sent to + and - neighbors.
 *          The outer boundaries are received from - and + neighbors.
 *          The inner boundaries of gauge field u is sent to + neighbors.
 *          The outer boundaries of gauge field u are received from - neighbors.
 *          This is what is needed e.g. for Dirac operators.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVectorPMGaugeP(qcd_vector *v, qcd_gaugeField *u);

/**
 * @fn int qcd_communicateVectorPGaugeP(qcd_vector *v, qcd_gaugeField *u)
 * @brief Sends spinor fields to + and gauge fields to + neighbors.
 * @details The inner boundaries of spinor field v is sent to + neighbors.
 *          The outer boundaries are received from - neighbors.
 *          The inner boundaries of gauge field u is sent to + neighbors.
 *          The outer boundaries of gauge field u are received from - neighbors.
 *          This is what is needed e.g. for the backward covariant derivative.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVectorPGaugeP(qcd_vector *v, qcd_gaugeField *u);

/**
 * @fn int qcd_communicateVectorPM(qcd_vector *v);
 * @brief Sends spinor fields to +/- neighbors.
 * @details The inner boundaries of spinor field v is sent to + and - neighbors.
 *          The outer boundaries are received from - and + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVectorPM(qcd_vector *v);

/**
 * @fn int qcd_communicateVectorM(qcd_vector *v);
 * @brief Sends spinor fields to - neighbors.
 * @details The inner boundaries of spinor field v is sent to - neighbors.
 *          The outer boundaries are received from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVectorM(qcd_vector *v);

/**
 * @fn int qcd_communicateVectorP(qcd_vector *v);
 * @brief Sends spinor fields to + neighbors.
 * @details The inner boundaries of spinor field v is sent to + neighbors.
 *          The outer boundaries are received from - neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVectorP(qcd_vector *v);

/**
 * @fn int qcd_communicateGaugeM(qcd_gaugeField *u)
 * @brief Sends gauge fields to - neighbors.
 * @details The inner boundaries of gauge field u is sent to - neighbors.
 *          The outer boundaries of gauge field u are received from + neighbors.
 *          This is what is needed e.g. for the computation of the plaquette.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGaugeM(qcd_gaugeField *u);

/**
 * @fn int qcd_communicateGaugeP(qcd_gaugeField *u)
 * @brief Sends gauge fields to + neighbors.
 * @details The inner boundaries of gauge field u is sent to + neighbors.
 *          The outer boundaries of gauge field u are received from - neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGaugeP(qcd_gaugeField *u);

/**
 * @fn int qcd_communicateGaugePInverse(qcd_gaugeField *u)
 * @brief Sends outer boundaries of gauge fields to + neighbors.
 * @details The outer boundaries of gauge field u is sent to + neighbors.
 *          The inner boundaries of gauge field u are received from - neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGaugePInverse(qcd_gaugeField *u);

/**
 * @fn int qcd_communicateGaugePM(qcd_gaugeField *u)
 * @brief Sends gauge fields to + and - neighbors.
 * @details The inner boundaries of gauge field u is sent to +  and - neighbors.
 *          The outer boundaries of gauge field u are received from - and + 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGaugePM(qcd_gaugeField *u);

/**
 * @fn int qcd_communicateCloverPM(qcd_cloverField *clv);
 * @brief Sends colver fields to + and - neighbors.
 * @details The inner boundaries of clover field clv is sent to + and - neighbors.
 *          The outer boundaries of clover field clv are received from - and + 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] clv pointer to an initialized qcd_cloverField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateCloverPM(qcd_cloverField *clv);

/**
 * @fn int qcd_communicateCloverP(qcd_cloverField *clv);
 * @brief Sends colver fields to + neighbors.
 * @details The inner boundaries of clover field clv is sent to + neighbors.
 *          The outer boundaries of clover field clv are received from - 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] clv pointer to an initialized qcd_cloverField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateCloverP(qcd_cloverField *clv);

/**
 * @fn int qcd_communicateCloverM(qcd_cloverField *clv);
 * @brief Sends colver fields to - neighbors.
 * @details The inner boundaries of clover field clv is sent to - neighbors.
 *          The outer boundaries of clover field clv are received from + 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] clv pointer to an initialized qcd_cloverField structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateCloverM(qcd_cloverField *clv);

/**
 * @fn int qcd_communicateTransformationM(qcd_gaugeTransformation *gt)
 * @brief Sends gauge transformations to - neighbors.
 * @details The inner boundaries of gauge transformation field gt is sent to
 *          - neighbors.
 *          The outer boundaries of gauge transformation field gt are 
 *          received from  + neighbors.
 *          Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] gt pointer to an initialized qcd_gaugeTransformation structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateTransformationM(qcd_gaugeTransformation *u);

/**
 * @fn int qcd_communicatePropagatorP(qcd_propagator *p)
 * @brief Sends propagator to + neighbors.
 * @details The inner boundaries of propagator p is sent to + neighbors.
 *          The outer boundaries of propagator p are received from - 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicatePropagatorP(qcd_propagator *p);

/**
 * @fn int qcd_communicatePropagatorPM(qcd_propagator *p)
 * @brief Sends propagator to + and - neighbors.
 * @details The inner boundaries of propagator p is sent to + and - neighbors.
 *          The outer boundaries of propagator p are received from - and + 
 *          neighbors. Non blocking communication is started and needs 
 *          to be finished with a call to qcd_waitall. All processes need to 
 *          call this function simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] p pointer to an initialized qcd_propagator structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicatePropagatorPM(qcd_propagator *p);

/**
 * @fn int qcd_communicateVector3dPM(qcd_vector3d *v);
 * @brief Sends 3d spinor fields to +/- neighbors.
 * @details The inner boundaries of time-slice spinor field v is sent 
 *          to + and - neighbors.
 *          The outer boundaries are received from - and + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVector3dPM(qcd_vector3d *v);

/**
 * @fn int qcd_communicateVector3dM(qcd_vector3d *v);
 * @brief Sends 3d spinor fields to - neighbors.
 * @details The inner boundaries of time-slice spinor field v is sent 
 *          to - neighbors.
 *          The outer boundaries are received from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVector3dM(qcd_vector3d *v);

/**
 * @fn int qcd_communicateVector3dP(qcd_vector3d *v);
 * @brief Sends 3d spinor fields to + neighbors.
 * @details The inner boundaries of time-slice spinor field v is sent 
 *          to +  neighbors.
 *          The outer boundaries are received from -  neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVector3dP(qcd_vector3d *v);

/**
 * @fn int qcd_communicateSpinorComponent3dPM(qcd_spinorComponent3d *sp);
 * @brief Sends 3d spinor-component fields to +/- neighbors.
 * @details The inner boundaries of time-slice spinor-component field sp is sent 
 *          to + and - neighbors.
 *          The outer boundaries are received from - and + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateSpinorComponent3dPM(qcd_spinorComponent3d *sp);

/**
 * @fn int qcd_communicateSpinorComponent3dM(qcd_spinorComponent3d *sp);
 * @brief Sends 3d spinor-component fields to - neighbors.
 * @details The inner boundaries of time-slice spinor-component field sp is sent 
 *          to - neighbors.
 *          The outer boundaries are received from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateSpinorComponent3dM(qcd_spinorComponent3d *sp);

/**
 * @fn int qcd_communicateSpinorComponent3dP(qcd_spinorComponent3d *sp);
 * @brief Sends 3d spinor-component fields to + neighbors.
 * @details The inner boundaries of time-slice spinor-component field sp is sent 
 *          to +  neighbors.
 *          The outer boundaries are received from -  neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] sp pointer to an initialized qcd_spinorComponent3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateSpinorComponent3dP(qcd_spinorComponent3d *sp);

/**
 * @fn int qcd_communicateLoops3dPM(qcd_loops3d *l);
 * @brief Sends 3d loops fields to +/- neighbors.
 * @details The inner boundaries of loops field l is sent 
 *          to + and - neighbors.
 *          The outer boundaries are received from - and + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets the appropriate boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] l pointer to an initialized qcd_loops3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateLoops3dPM(qcd_loops3d *l);

/**
 * @fn int qcd_communicateLoops3dM(qcd_loops3d *l);
 * @brief Sends 3d loops fields to - neighbors.
 * @details The inner boundaries of loops field l is sent 
 *          to - neighbors.
 *          The outer boundaries are received from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets the appropriate boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] l pointer to an initialized qcd_loops3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateLoops3dM(qcd_loops3d *l);

/**
 * @fn int qcd_communicateLoops3dP(qcd_loops3d *l);
 * @brief Sends 3d loops fields to + neighbors.
 * @details The inner boundaries of time-slice loops field l is sent 
 *          to +  neighbors.
 *          The outer boundaries are received from -  neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets the appropriate boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] l pointer to an initialized qcd_loops3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateLoops3dP(qcd_loops3d *l);


/**
 * @fn int qcd_communicateLoops3dPMGauge3dP(qcd_loops3d *l, qcd_gaugeField3d *u)
 * @brief Sends 3d loops fields to + and - neighbors and 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          The inner boundaries of time-slice loops fields are sent to + and - neighbors
 *          The outer boundaries of time-slice loops fields are received from + and - neighbors
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] l pointer to an initialized qcd_loops3d structure
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateLoops3dPMGauge3dP(qcd_loops3d *l, qcd_gaugeField3d *u);

/**
 * @fn int qcd_communicateGauge3dM(qcd_gaugeField3d *u)
 * @brief Sends 3d gauge fields to - neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to - neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGauge3dM(qcd_gaugeField3d *u);

/**
 * @fn int qcd_communicateGauge3dP(qcd_gaugeField3d *u)
 * @brief Sends 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGauge3dP(qcd_gaugeField3d *u);



/**
 * @fn int qcd_communicateVector3dPGauge3dP(qcd_vector3d *v, qcd_gaugeField3d *u)
 * @brief Sends 3d spinor fields to +  neighbors and 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          The inner boundaries of time-slice spinor fields are sent to +  neighbors
 *          The outer boundaries of time-slice spinor fields are received from - neighbors
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVector3dPGauge3dP(qcd_vector3d *v, qcd_gaugeField3d *u);


/**
 * @fn int qcd_communicateVector3dPMGauge3dP(qcd_vector3d *v, qcd_gaugeField3d *u)
 * @brief Sends 3d spinor fields to + and - neighbors and 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          The inner boundaries of time-slice spinor fields are sent to + and - neighbors
 *          The outer boundaries of time-slice spinor fields are received from + and - neighbors
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] v pointer to an initialized qcd_vector3d structure
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateVector3dPMGauge3dP(qcd_vector3d *v, qcd_gaugeField3d *u);

/**
 * @fn int qcd_communicateSpinorComponent3dPGauge3dP(qcd_spinorComponent3d *v, qcd_gaugeField3d *u)
 * @brief Sends 3d spinor component fields to +  neighbors and 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          The inner boundaries of time-slice spinor component fields are sent to +  neighbors
 *          The outer boundaries of time-slice spinor component fields are received from - neighbors
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateSpinorComponent3dPGauge3dP(qcd_spinorComponent3d *v, qcd_gaugeField3d *u);

/**
 * @fn int qcd_communicateSpinorComponent3dPMGauge3dP(qcd_spinorComponent3d *v, qcd_gaugeField3d *u)
 * @brief Sends 3d spinor component fields to + and - neighbors and 3d gauge fields to + neighbors.
 * @details The inner boundaries of time-slice gauge field u is sent to + neighbors.
 *          The outer boundaries of time-slice gauge field u are received 
 *          from - neighbors.
 *          The inner boundaries of time-slice spinor component fields are sent to + and - neighbors
 *          The outer boundaries of time-slice spinor component fields are received from -  and + neighbors
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeField3d structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateSpinorComponent3dPMGauge3dP(qcd_spinorComponent3d *v, qcd_gaugeField3d *u);

/**
 * @fn int qcd_communicateGaugeOM(qcd_gaugeFieldO *u)
 * @brief Sends gauge fields to - neighbors.
 * @details The inner boundaries of gauge field u (in openQCD layout) is sent 
 *          to - neighbors.
 *          The outer boundaries of gauge field u are received 
 *          from + neighbors.
 *          Non blocking communication is started and needs to be finished with
 *          a call to qcd_waitall. All processes need to call this function
 *          simultaneously.
 *          This function sets appropriate the boundary-up-to-date-flags to true,
 *          which will become true only after the successful call to qcd_waitall!
 * @param[in,out] u pointer to an initialized qcd_gaugeFieldO structure
 * @return 0 on success, 1 on failure
 */
int qcd_communicateGaugeOM(qcd_gaugeFieldO *u);

/**
 * @fn void qcd_waitall(qcd_geometry *geo)
 * @brief Waits untill all non-blocking communication is finished.
 * @details The MPI communication of fields, that use the geometry geo, is
 *          finished with a call to this function.
 *          All processes have to call this function simultaneously.
 * @param geo pointer to an initialized qcd_geometry structure
 */
void qcd_waitall(qcd_geometry *geo);

#endif
