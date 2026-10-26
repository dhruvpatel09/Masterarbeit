/* Local observation-equivalence fixture: diagonal operator and real LAPACK DSYEV. */
#include <qcd.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>

int qcd_initSpinorComponent3d(qcd_spinorComponent3d *v,qcd_geometry *g)
{v->geo=g;v->D=calloc(g->lV3,sizeof(*v->D));v->BMup2date=v->BPup2date=0;return v->D==NULL;}
void qcd_destroySpinorComponent3d(qcd_spinorComponent3d *v){free(v->D);}
void qcd_zeroSpinorComponent3d(qcd_spinorComponent3d *v){memset(v->D,0,v->geo->lV3*sizeof(*v->D));}
void qcd_setSpinorComponent3d(qcd_spinorComponent3d *v,qcd_complex_16 a)
{for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++)v->D[i][j]=a;}
void qcd_copySpinorComponent3d(qcd_spinorComponent3d *v,qcd_spinorComponent3d *w)
{memmove(v->D,w->D,v->geo->lV3*sizeof(*v->D));}
void qcd_scaleSpinorComponent3d(qcd_spinorComponent3d *v,double a)
{for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++){v->D[i][j].re*=a;v->D[i][j].im*=a;}}
void qcd_raxpySpinorComponent3d(qcd_spinorComponent3d *v,double a,qcd_spinorComponent3d *x,qcd_spinorComponent3d *y)
{for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++){double re=a*x->D[i][j].re+y->D[i][j].re,im=a*x->D[i][j].im+y->D[i][j].im;v->D[i][j]=(qcd_complex_16){re,im};}}
void qcd_axpySpinorComponent3d(qcd_spinorComponent3d *v,qcd_complex_16 a,qcd_spinorComponent3d *x,qcd_spinorComponent3d *y)
{for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++){qcd_complex_16 b=x->D[i][j],c=y->D[i][j];v->D[i][j]=(qcd_complex_16){a.re*b.re-a.im*b.im+c.re,a.re*b.im+a.im*b.re+c.im};}}
void qcd_subSpinorComponent3d(qcd_spinorComponent3d *v,qcd_spinorComponent3d *x,qcd_spinorComponent3d *y)
{qcd_raxpySpinorComponent3d(v,-1.,y,x);}
void qcd_addSpinorComponent3d(qcd_spinorComponent3d *v,qcd_spinorComponent3d *x,qcd_spinorComponent3d *y)
{qcd_raxpySpinorComponent3d(v,1.,y,x);}
void qcd_mulSpinorComponent3dC(qcd_spinorComponent3d *v,qcd_complex_16 a)
{for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++){qcd_complex_16 b=v->D[i][j];v->D[i][j]=(qcd_complex_16){a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};}}
qcd_complex_16 qcd_mulAdjointSpinorComponent3dSpinorComponent3d(qcd_spinorComponent3d *v,qcd_spinorComponent3d *w)
{qcd_complex_16 out={0.,0.};for(long i=0;i<v->geo->lV3;i++)for(int j=0;j<3;j++){qcd_complex_16 a=v->D[i][j],b=w->D[i][j];out.re+=a.re*b.re+a.im*b.im;out.im+=a.re*b.im-a.im*b.re;}return out;}
double qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(qcd_spinorComponent3d *v,qcd_spinorComponent3d *w)
{return qcd_mulAdjointSpinorComponent3dSpinorComponent3d(v,w).re;}
double qcd_normSpinorComponent3d(qcd_spinorComponent3d *v)
{return sqrt(qcd_mulAdjointSpinorComponent3dSpinorComponent3dR(v,v));}
int qcd_applyCovariantLaplaceSpinorComponent3d(qcd_spinorComponent3d *out,qcd_spinorComponent3d *in,qcd_gaugeField3d *g,double mm)
{(void)g;for(long i=0;i<in->geo->lV3;i++)for(int j=0;j<3;j++){long k=3*i+j;double d=k<3?.01*(k+1):.2+10.*(k-3)/(3*in->geo->lV3-3);out->D[i][j].re=(d+mm)*in->D[i][j].re;out->D[i][j].im=(d+mm)*in->D[i][j].im;}return 0;}
void qcd_rand(double *x,long n){for(long i=0;i<n;i++)x[i]=sin((double)i+1.);}

extern void dsyev_(char *,char *,int *,double *,int *,double *,double *,int *,int *);
static int dense(double *a,int n,double *eval,double *vectors)
{
    char job='V',triangle='U';int info,lwork=-1;double query;
    dsyev_(&job,&triangle,&n,a,&n,eval,&query,&lwork,&info);
    if(info)return info;
    lwork=(int)query;double *work=malloc(lwork*sizeof(*work));
    dsyev_(&job,&triangle,&n,a,&n,eval,work,&lwork,&info);free(work);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)vectors[i*n+j]=a[j*n+i];
    return info;
}
int LAPACKE_dsyevr(int layout,char job,char range,char triangle,int n,double *a,int lda,double vl,double vu,int il,int iu,double tol,int *m,double *eval,double *vectors,int ldz,int *supp)
{
    (void)layout;(void)job;(void)range;(void)triangle;(void)lda;(void)vl;(void)vu;(void)il;(void)iu;(void)tol;(void)ldz;(void)supp;
    const char *failure=getenv("QX_TEST_REDUCED_FAILURE");if(failure && strcmp(failure,"DSYEVR")==0)return -73;
    *m=n;return dense(a,n,eval,vectors);
}
int LAPACKE_dstemr(int layout,char job,char range,int n,double *d,double *e,double vl,double vu,int il,int iu,int *m,double *eval,double *vectors,int ldz,int nzc,int *supp,int *tryrac)
{
    (void)layout;(void)job;(void)range;(void)vl;(void)vu;(void)il;(void)iu;(void)ldz;(void)nzc;(void)supp;(void)tryrac;
    const char *failure=getenv("QX_TEST_REDUCED_FAILURE");if(failure && strcmp(failure,"DSTEMR")==0)return -71;
    double *a=calloc(n*n,sizeof(*a));for(int i=0;i<n;i++){a[i*n+i]=d[i];if(i<n-1)a[i*n+i+1]=a[(i+1)*n+i]=e[i];}
    *m=n;int result=dense(a,n,eval,vectors);free(a);return result;
}
