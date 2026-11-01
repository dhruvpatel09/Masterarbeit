/* Execute the exact extracted kernels and face shifts on a simulated 4^3 grid. */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#define NV 4
#define NS 64
typedef struct {double re,im;} qcd_complex_16;
typedef struct {int Pos[4],lL[4];long lV3;} qcd_geometry;
#define qcd_CADD(a,b) ((qcd_complex_16){(a).re+(b).re,(a).im+(b).im})
#define qcd_CMUL(a,b) ((qcd_complex_16){(a).re*(b).re-(a).im*(b).im,(a).re*(b).im+(a).im*(b).re})
#define qcd_CONJ(a) ((qcd_complex_16){(a).re,-(a).im})
#define qcd_LEXIC0(x,y,z,L) ((long)((z)+(L)[3]*((y)+(L)[2]*(x))))
#include "kernels.h"
#include "mpi_helpers.h"

static qcd_geometry geometry;
static int current_axis,current_displacement,minus_rank,plus_rank;

static void fill(endpoint_packet *p,int x,int y,int z) {
    memset(p,0,sizeof(*p));p->xyz[0]=x;p->xyz[1]=y;p->xyz[2]=z;
    double base=1.+x+64.*y+4096.*z;
    for(int t=0;t<2;t++)for(int n=0;n<NV;n++)for(int c=0;c<3;c++) {
        p->v[t][n][c].re=base+.125*(t*12+n*3+c);
        p->v[t][n][c].im=-base+.0625*(t*12+n*3+c);
    }
    for(int a=0;a<3;a++)for(int b=0;b<3;b++) {
        p->p[a][b].re=base+9*a+b;p->p[a][b].im=-base+a-3*b;
    }
}

int MPI_Abort(MPI_Comm comm,int code){(void)comm;fprintf(stderr,"unexpected MPI_Abort\n");exit(code);return 0;}
int MPI_Comm_size(MPI_Comm comm,int *size){(void)comm;*size=64;return 0;}
int MPI_Allreduce(const void *a,void *b,int n,MPI_Datatype type,int op,MPI_Comm comm) {
    (void)op;(void)comm;assert(type==MPI_INT);memcpy(b,a,n*sizeof(int));return 0;
}
int MPI_Allgather(const void *a,int n,MPI_Datatype type,void *b,int m,MPI_Datatype other,MPI_Comm comm) {
    (void)a;(void)comm;assert(n==4&&m==4&&type==MPI_INT&&other==MPI_INT);
    int (*positions)[4]=b;
    for(int rank=0;rank<64;rank++){positions[rank][0]=0;positions[rank][1]=rank/16;positions[rank][2]=(rank/4)%4;positions[rank][3]=rank%4;}
    return 0;
}
int MPI_Sendrecv(const void *a,int n,MPI_Datatype type,int to,int tag,void *b,int m,MPI_Datatype other,int from,int tag2,MPI_Comm comm,MPI_Status *status) {
    (void)comm;(void)status;
    assert(n==256*(int)sizeof(endpoint_packet)&&m==n&&type==MPI_BYTE&&other==MPI_BYTE);
    assert(to==minus_rank&&from==plus_rank&&tag==170+current_axis&&tag2==tag);
    const endpoint_packet *send=a;endpoint_packet *recv=b;
    int axis=current_axis,jaxis=(axis+1)%3,kaxis=(axis+2)%3;
    for(int j=0;j<16;j++)for(int k=0;k<16;k++) {
        int xyz[3]={16*geometry.Pos[1],16*geometry.Pos[2],16*geometry.Pos[3]};
        xyz[jaxis]+=j;xyz[kaxis]+=k;xyz[axis]=(xyz[axis]+current_displacement)%64;
        endpoint_packet expected;fill(&expected,xyz[0],xyz[1],xyz[2]);
        assert(memcmp(send+j*16+k,&expected,sizeof(expected))==0);
        xyz[axis]=(xyz[axis]+16)%64;fill(recv+j*16+k,xyz[0],xyz[1],xyz[2]);
    }
    return 0;
}

static void shifts(void) {
    endpoint_packet *base=calloc(4096,sizeof(*base)),*a=malloc(4096*sizeof(*a)),*b=malloc(4096*sizeof(*b));
    endpoint_packet *send=malloc(256*sizeof(*send)),*recv=malloc(256*sizeof(*recv));
    assert(base&&a&&b&&send&&recv);
    for(int corner=0;corner<2;corner++) {
        geometry.Pos[0]=0;geometry.lL[0]=192;geometry.lV3=4096;
        for(int k=1;k<4;k++){geometry.Pos[k]=corner?3:0;geometry.lL[k]=16;}
        for(int x=0;x<16;x++)for(int y=0;y<16;y++)for(int z=0;z<16;z++)
            fill(base+qcd_LEXIC0(x,y,z,geometry.lL),x+16*geometry.Pos[1],y+16*geometry.Pos[2],z+16*geometry.Pos[3]);
        int minus[3],plus[3];rank_neighbours(&geometry,minus,plus);
        for(int axis=0;axis<3;axis++) {
            current_axis=axis;minus_rank=minus[axis];plus_rank=plus[axis];memcpy(a,base,4096*sizeof(*a));
            for(int r=1;r<=64;r++) {
                current_displacement=r-1;shift_endpoint(a,b,send,recv,&geometry,axis,minus_rank,plus_rank);
                endpoint_packet *tmp=a;a=b;b=tmp;verify_shift(base,a,4096,axis,r);
                if(r==1||r==16||r==17||r==31||r==32||r==64)for(int site=0;site<4096;site++) {
                    int xyz[3];for(int k=0;k<3;k++)xyz[k]=(base[site].xyz[k]+(k==axis?r:0))%64;
                    endpoint_packet expected;fill(&expected,xyz[0],xyz[1],xyz[2]);
                    assert(memcmp(a+site,&expected,sizeof(expected))==0);
                }
            }
            assert(memcmp(base,a,4096*sizeof(*a))==0);
        }
    }
    free(base);free(a);free(b);free(send);free(recv);
}

static unsigned int seed=1234567;
static double random_value(void){seed=1664525U*seed+1013904223U;return ((seed>>8)/16777216.)-.5;}
static void contractions(void) {
    for(active_modes=1;active_modes<=NV;active_modes++)for(int trial=0;trial<100;trial++) {
        endpoint_packet x={0},y={0};
        for(int t=0;t<2;t++)for(int n=0;n<NV;n++)for(int c=0;c<3;c++) {
            x.v[t][n][c]=(qcd_complex_16){random_value(),random_value()};
            y.v[t][n][c]=(qcd_complex_16){random_value(),random_value()};
        }
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
            x.p[i][j]=(qcd_complex_16){random_value(),random_value()};y.p[i][j]=(qcd_complex_16){random_value(),random_value()};
        }
        qcd_complex_16 fast=packet_pair(&x,&y),direct=packet_direct(&x,&y),reverse=packet_pair(&y,&x);
        assert(hypot(fast.re-direct.re,fast.im-direct.im)<1e-12);
        assert(hypot(fast.re-reverse.re,fast.im+reverse.im)<1e-12);
        qcd_complex_16 self=packet_pair(&x,&x);assert(self.re>=-1e-12&&fabs(self.im)<1e-12);
    }
}

static void analytic_cutoff(void) {
    endpoint_packet x={0};
    for(int n=0;n<NV;n++)for(int t=0;t<2;t++)x.v[t][n][n%3].re=n+1;
    for(int c=0;c<3;c++)x.p[c][c].re=1;
    for(active_modes=1;active_modes<=NV;active_modes++) {
        double weights[3]={0},expected=0;
        for(int n=0;n<active_modes;n++)weights[n%3]+=(n+1.)*(n+1.);
        for(int c=0;c<3;c++)expected+=weights[c]*weights[c];
        qcd_complex_16 value=packet_pair(&x,&x);
        assert(value.re==expected&&value.im==0);
    }
}

int main(void){analytic_cutoff();contractions();shifts();puts("PASS every native mode cutoff against independent quadratic contraction, conjugation and complete periodic payload shifts");return 0;}
