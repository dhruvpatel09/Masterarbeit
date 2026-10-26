/* Test facade backed by Fortran DSYEV; production uses its pinned MKL facade. */
#ifndef QX_TEST_LAPACKE_H
#define QX_TEST_LAPACKE_H
typedef int lapack_int;
typedef int lapack_logical;
#define LAPACK_ROW_MAJOR 101
int LAPACKE_dsyevr(int,char,char,char,int,double *,int,double,double,int,int,double,int *,double *,double *,int,int *);
int LAPACKE_dstemr(int,char,char,int,double *,double *,double,double,int,int,int *,double *,double *,int,int,int *,int *);
#endif
