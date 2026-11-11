/* Independent time-dependent SU(3) paths on a 192x1x1x1 synthetic field. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {double re,im;} qcd_complex_16;
typedef struct {int lL[4];} qcd_geometry;
typedef struct {qcd_geometry *geo;qcd_complex_16 (*D)[4][3][3];} qcd_gaugeField;
#define qcd_LEXIC(t,x,y,z,L) ((long)(t))
#define qcd_CADD(a,b) ((qcd_complex_16){(a).re+(b).re,(a).im+(b).im})
#define qcd_CMUL(a,b) ((qcd_complex_16){(a).re*(b).re-(a).im*(b).im,(a).re*(b).im+(a).im*(b).re})
#include "time_functions.h"
#include "source_time_plan.h"

int main(void) {
    qcd_geometry g={{192,1,1,1}};
    qcd_gaugeField u={&g,calloc(192,sizeof(*u.D))};assert(u.D);
    for(int t=0;t<192;t++) {
        double theta=.002*(t+1),alpha=.001*(t+1)*(t+1);
        qcd_complex_16 positive={cos(alpha),sin(alpha)},negative={cos(alpha),-sin(alpha)};
        u.D[t][0][0][0]=(qcd_complex_16){cos(theta)*positive.re,cos(theta)*positive.im};
        u.D[t][0][1][0]=(qcd_complex_16){sin(theta)*positive.re,sin(theta)*positive.im};
        u.D[t][0][0][1]=(qcd_complex_16){-sin(theta)*negative.re,-sin(theta)*negative.im};
        u.D[t][0][1][1]=(qcd_complex_16){cos(theta)*negative.re,cos(theta)*negative.im};u.D[t][0][2][2].re=1;
    }
    for(int sample=0;sample<3;sample++)for(int tau=0;tau<=32;tau++) {
        int first=(int[]){58,90,122}[sample],end=first+tau;qcd_complex_16 p[3][3];
        temporal_product(&u,0,0,0,first,end,p);
        for(int j=0;j<3;j++) {
            qcd_complex_16 v[3]={{0}};v[j].re=1;
            for(int t=end-1;t>=first;t--) {
                qcd_complex_16 next[3]={{0}};
                for(int i=0;i<3;i++)for(int k=0;k<3;k++)next[i]=qcd_CADD(next[i],qcd_CMUL(u.D[t][0][i][k],v[k]));
                memcpy(v,next,sizeof(v));
            }
            for(int i=0;i<3;i++)assert(hypot(p[i][j].re-v[i].re,p[i][j].im-v[i].im)<3e-12);
        }
    }
    for(int i=0;i<3;i++) {
        int cfg=(int[]){170,640,650}[i];assert(source_time_scan_allowed(cfg,10,58));assert(source_time_scan_allowed(cfg,10,122));
        assert(source_time_scan_allowed(cfg,10,90)==(cfg==650));assert(source_time_scan_allowed(cfg,40,122)==(cfg==650));
        assert(!source_time_scan_allowed(cfg,40,58));assert(!source_time_scan_allowed(cfg,4,90));
        assert(!source_time_scan_allowed(cfg,10,-1));assert(!source_time_scan_allowed(cfg,10,160));
    }
    assert(!source_time_scan_allowed(100,10,58));free(u.D);
    puts("PASS independent ordered temporal paths at 58/90/122, all 33 separations and pilot admission guards");return 0;
}
