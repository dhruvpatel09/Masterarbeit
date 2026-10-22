#ifndef QX_CLOVER_TRACE_H
#define QX_CLOVER_TRACE_H

/* F=(Q-Q^dagger)/8. All densities below sum the six unordered planes;
   there is no additional factor 1/2. T=|Tr F|^2/3 is the trace component. */
typedef struct {
    double raw, trace_component, traceless_direct;
} qx_clover_components;

static qx_clover_components qx_clover_plane(
    const double qr[3][3], const double qi[3][3])
{
    double fr[3][3], fi[3][3], tr=0.0, ti=0.0;
    qx_clover_components result={0.0,0.0,0.0};
    for(int a=0;a<3;++a) for(int b=0;b<3;++b) {
        fr[a][b]=0.125*(qr[a][b]-qr[b][a]);
        fi[a][b]=0.125*(qi[a][b]+qi[b][a]);
        result.raw+=fr[a][b]*fr[a][b]+fi[a][b]*fi[a][b];
    }
    for(int a=0;a<3;++a) { tr+=fr[a][a]; ti+=fi[a][a]; }
    result.trace_component=(tr*tr+ti*ti)/3.0;
    for(int a=0;a<3;++a) for(int b=0;b<3;++b) {
        double r=fr[a][b]-(a==b?tr/3.0:0.0);
        double i=fi[a][b]-(a==b?ti/3.0:0.0);
        result.traceless_direct+=r*r+i*i;
    }
    return result;
}
#endif
