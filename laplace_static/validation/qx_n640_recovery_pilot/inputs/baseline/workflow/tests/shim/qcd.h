/* Dense diagonal-operator fixture, not the production qcd API implementation. */
#ifndef QX_TEST_QCD_H
#define QX_TEST_QCD_H
#include <math.h>
#include <mpi.h>
#define qcd_USE_LAPACKE 1
typedef struct {double re,im;} qcd_complex_16;
typedef struct {int myid,Pos[4],lL[4],root3d[1],comm3d;long lV3;} qcd_geometry;
typedef struct {qcd_complex_16 (*D)[3];qcd_geometry *geo;int BMup2date,BPup2date;} qcd_spinorComponent3d;
typedef struct {qcd_geometry *geo;} qcd_gaugeField3d;
#include <qcd_lanczos.h>
#define qcd_CARG(x) atan2((x).im,(x).re)
int qcd_initSpinorComponent3d(qcd_spinorComponent3d *,qcd_geometry *);
void qcd_destroySpinorComponent3d(qcd_spinorComponent3d *);
void qcd_zeroSpinorComponent3d(qcd_spinorComponent3d *);
void qcd_setSpinorComponent3d(qcd_spinorComponent3d *,qcd_complex_16);
void qcd_copySpinorComponent3d(qcd_spinorComponent3d *,qcd_spinorComponent3d *);
void qcd_scaleSpinorComponent3d(qcd_spinorComponent3d *,double);
void qcd_raxpySpinorComponent3d(qcd_spinorComponent3d *,double,qcd_spinorComponent3d *,qcd_spinorComponent3d *);
void qcd_axpySpinorComponent3d(qcd_spinorComponent3d *,qcd_complex_16,qcd_spinorComponent3d *,qcd_spinorComponent3d *);
void qcd_subSpinorComponent3d(qcd_spinorComponent3d *,qcd_spinorComponent3d *,qcd_spinorComponent3d *);
void qcd_addSpinorComponent3d(qcd_spinorComponent3d *,qcd_spinorComponent3d *,qcd_spinorComponent3d *);
void qcd_mulSpinorComponent3dC(qcd_spinorComponent3d *,qcd_complex_16);
qcd_complex_16 qcd_mulAdjointSpinorComponent3dSpinorComponent3d(qcd_spinorComponent3d *,qcd_spinorComponent3d *);
double qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(qcd_spinorComponent3d *,qcd_spinorComponent3d *);
double qcd_normSpinorComponent3d(qcd_spinorComponent3d *);
int qcd_applyCovariantLaplaceSpinorComponent3d(qcd_spinorComponent3d *,qcd_spinorComponent3d *,qcd_gaugeField3d *,double);
void qcd_rand(double *,long);
#endif
