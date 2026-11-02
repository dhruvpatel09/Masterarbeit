/* Independent rectangular sums, source/sink orientation and symmetry controls. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef NV
#define NV 4
#endif
typedef struct {double re,im;} qcd_complex_16;
#define qcd_CADD(a,b) ((qcd_complex_16){(a).re+(b).re,(a).im+(b).im})
#define qcd_CMUL(a,b) ((qcd_complex_16){(a).re*(b).re-(a).im*(b).im,(a).re*(b).im+(a).im*(b).re})
#define qcd_CONJ(a) ((qcd_complex_16){(a).re,-(a).im})
#include "kernels.h"
#include "mixed_contractions.h"

static unsigned int seed=88721;
static double random_value(void){seed=1664525U*seed+1013904223U;return ((seed>>8)/16777216.)-.5;}
static void close_value(qcd_complex_16 x,qcd_complex_16 y) {
    double error=hypot(x.re-y.re,x.im-y.im),scale=fmax(hypot(x.re,x.im),hypot(y.re,y.im));
    assert(error<=1e-12+1e-11*scale);
}
static void fill(endpoint_packet *x) {
    memset(x,0,sizeof(*x));
    for(int t=0;t<2;t++)for(int n=0;n<NV;n++)for(int c=0;c<3;c++)x->v[t][n][c]=(qcd_complex_16){random_value(),random_value()};
    for(int a=0;a<3;a++)for(int b=0;b<3;b++)x->p[a][b]=(qcd_complex_16){random_value(),random_value()};
}
static qcd_complex_16 brute(const endpoint_packet *x,const endpoint_packet *y,int src,int sink) {
    qcd_complex_16 value={0};
    for(int i=0;i<src;i++)for(int j=0;j<sink;j++) {
        qcd_complex_16 a={0},b={0};
        for(int c=0;c<3;c++)for(int d=0;d<3;d++) {
            a=qcd_CADD(a,qcd_CMUL(qcd_CONJ(x->v[0][i][c]),qcd_CMUL(x->p[c][d],x->v[1][j][d])));
            b=qcd_CADD(b,qcd_CMUL(qcd_CONJ(y->v[0][i][c]),qcd_CMUL(y->p[c][d],y->v[1][j][d])));
        }
        value=qcd_CADD(value,qcd_CMUL(b,qcd_CONJ(a)));
    }
    return value;
}
static qcd_complex_16 phase(double a){return (qcd_complex_16){cos(a),sin(a)};}
static void unitary(double offset,qcd_complex_16 u[3][3]) {
    double theta=.3+offset,phi=.7-offset;
    double real[3][3]={{cos(theta),-sin(theta)*cos(phi),sin(theta)*sin(phi)},
        {sin(theta),cos(theta)*cos(phi),-cos(theta)*sin(phi)},{0,sin(phi),cos(phi)}};
    double angles[3]={offset,.7,-offset-.7};
    for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
        qcd_complex_16 z=phase(angles[i]);u[i][j]=(qcd_complex_16){real[i][j]*z.re,real[i][j]*z.im};
    }
}
static void gauge_change(endpoint_packet *x,double offset) {
    qcd_complex_16 u[2][3][3];unitary(offset,u[0]);unitary(offset+.29,u[1]);
    for(int t=0;t<2;t++)for(int n=0;n<NV;n++) {
        qcd_complex_16 v[3]={0};
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)v[i]=qcd_CADD(v[i],qcd_CMUL(u[t][i][j],x->v[t][n][j]));
        memcpy(x->v[t][n],v,sizeof(v));
    }
    qcd_complex_16 result[3][3]={0};
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)
        result[i][j]=qcd_CADD(result[i][j],qcd_CMUL(u[0][i][a],qcd_CMUL(x->p[a][b],qcd_CONJ(u[1][j][b]))));
    memcpy(x->p,result,sizeof(result));
}
static void contractions(const int cuts[4]) {
    for(int trial=0;trial<(NV==100?5:50);trial++) {
        endpoint_packet x,y;fill(&x);fill(&y);qcd_complex_16 a[4][4],b[4][4],changed[4][4];
        assert(!packet_mixed_all(&x,&y,cuts,4,a));assert(!packet_mixed_all(&y,&x,cuts,4,b));
        for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
            close_value(a[i][j],brute(&x,&y,cuts[i],cuts[j]));
            close_value(a[i][j],packet_mixed_direct(&x,&y,cuts[i],cuts[j]));
            close_value(a[i][j],qcd_CONJ(b[i][j]));
            if(i==j){active_modes=cuts[i];close_value(a[i][j],packet_pair(&x,&y));}
        }
        gauge_change(&x,.13);gauge_change(&y,.39);assert(!packet_mixed_all(&x,&y,cuts,4,changed));
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)close_value(a[i][j],changed[i][j]);
    }
}
static void orientation(void) {
    endpoint_packet x={0};qcd_complex_16 a[4][4];int cuts[4]={1,2,3,4};
    x.v[0][0][0].re=1;x.v[1][2][0].re=2;
    for(int c=0;c<3;c++)x.p[c][c].re=1;
    assert(!packet_mixed_all(&x,&x,cuts,4,a));
    assert(a[0][2].re==4&&a[2][0].re==0);
    assert(a[0][2].im==0&&a[2][0].im==0);
}
static void equal_time(const int cuts[4]) {
    endpoint_packet x,y;fill(&x);fill(&y);qcd_complex_16 a[4][4];
    for(int n=0;n<NV;n++)for(int c=0;c<3;c++){x.v[1][n][c]=x.v[0][n][c];y.v[1][n][c]=y.v[0][n][c];}
    mat3_zero(x.p);mat3_zero(y.p);for(int c=0;c<3;c++){x.p[c][c].re=1;y.p[c][c].re=1;}
    assert(!packet_mixed_all(&x,&y,cuts,4,a));
    for(int i=0;i<4;i++){assert(a[i][i].re>=-1e-12);assert(fabs(a[i][i].im)<1e-10);for(int j=0;j<4;j++)close_value(a[i][j],qcd_CONJ(a[j][i]));}
    for(int trial=0;trial<30;trial++) {
        qcd_complex_16 z[4],v={0};for(int i=0;i<4;i++)z[i]=(qcd_complex_16){random_value(),random_value()};
        for(int i=0;i<4;i++)for(int j=0;j<4;j++)v=qcd_CADD(v,qcd_CMUL(qcd_CONJ(z[i]),qcd_CMUL(a[i][j],z[j])));
        assert(v.re>=-1e-10);assert(fabs(v.im)<1e-9);
    }
}
int main(void) {
    int cuts[4];for(int i=0;i<4;i++)cuts[i]=NV==100?(int[]){40,60,80,100}[i]:i+1;
    contractions(cuts);orientation();equal_time(cuts);
    endpoint_packet x={0};qcd_complex_16 a[4][4];int invalid[4]={1,2,2,4};
    assert(packet_mixed_all(&x,&x,invalid,4,a));assert(packet_mixed_all(&x,&x,cuts,0,a));
    assert(!isfinite(packet_mixed_direct(&x,&x,0,1).re));
    printf("PASS mixed contractions NV=%d: all rectangles, independent sums, orientation, gauge covariance, legacy diagonals and equal-time Gram checks\n",NV);
    return 0;
}
