#ifndef SAVED_MODE_METRICS_H
#define SAVED_MODE_METRICS_H
#include <math.h>
#include <stddef.h>
/* Arrays are [re,im] pairs. Accumulate locally in long double, then MPI_SUM. */
static void mode_sums(const double *v, const double *av, size_t n,
                      double lambda, double out[4])
{
    long double norm=0, residual=0, ray_re=0, ray_im=0;
    for (size_t k=0;k<n;k++) {
        const long double vr=v[2*k], vi=v[2*k+1];
        const long double ar=av[2*k], ai=av[2*k+1];
        const long double rr=ar-lambda*vr, ri=ai-lambda*vi;
        norm+=vr*vr+vi*vi;
        residual+=rr*rr+ri*ri;
        ray_re+=vr*ar+vi*ai;
        ray_im+=vr*ai-vi*ar;
    }
    out[0]=(double)norm; out[1]=(double)residual;
    out[2]=(double)ray_re; out[3]=(double)ray_im;
}
static void overlap_sums(const double *a, const double *b, size_t n,
                         double out[2])
{
    long double re=0, im=0;
    for (size_t k=0;k<n;k++) {
        re+=(long double)a[2*k]*b[2*k]+(long double)a[2*k+1]*b[2*k+1];
        im+=(long double)a[2*k]*b[2*k+1]-(long double)a[2*k+1]*b[2*k];
    }
    out[0]=(double)re; out[1]=(double)im;
}
#endif

