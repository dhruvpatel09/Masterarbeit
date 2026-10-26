/** @file qcd_loopsOps.h
 *  @brief header file for qcd_loopsOps.c
 *  @details Details & license
 *  @version 1.0
 *  @author Tomasz Korzec
 *  @date 2024
 *  @copyright GNU Public License.
 */
 
#ifndef H_QCD_loops
#define H_QCD_loops 1

/******************************************************************************
 * prototypes 
 ******************************************************************************/


void qcd_zeroLoops3d(qcd_loops3d *l);

void qcd_unitLoops3d(qcd_loops3d *l);

void qcd_scaleLoops3d(qcd_loops3d *l, double r);

void qcd_mulLoops3dC(qcd_loops3d *l, qcd_complex_16 c);

void qcd_addLoops3d(qcd_loops3d *sum, qcd_loops3d *summand1, qcd_loops3d *summand2);

qcd_complex_16 qcd_dotLoops(qcd_loops3d *fac1, qcd_loops3d *fac2);

double qcd_normLoops(qcd_loops3d *l);

int qcd_LaplaceLoops(qcd_loops3d *Dl, qcd_loops3d *l, qcd_gaugeField3d *u);

#endif
