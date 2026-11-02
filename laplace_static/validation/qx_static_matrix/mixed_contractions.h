#ifndef QX_MIXED_CONTRACTIONS_H
#define QX_MIXED_CONTRACTIONS_H
#include <math.h>
#include <string.h>
#define MIXED_MAX_CUTOFFS 4

/* out[source_slot][sink_slot]; source modes live at t0, sink modes at t1. */
static int packet_mixed_all(const endpoint_packet *x,const endpoint_packet *y,
    const int *cutoffs,int count,qcd_complex_16 out[MIXED_MAX_CUTOFFS][MIXED_MAX_CUTOFFS])
{
    if(count<1||count>MIXED_MAX_CUTOFFS)return 1;
    for(int k=0;k<count;k++)if(cutoffs[k]<1||cutoffs[k]>NV||(k&&cutoffs[k]<=cutoffs[k-1]))return 1;
    qcd_complex_16 k0[3][3],k1[3][3],sources[4][3][3],sinks[4][3][3];
    mat3_zero(k0);mat3_zero(k1);int slot=0;
    for(int n=0;n<cutoffs[count-1];n++) {
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
            k1[i][j]=qcd_CADD(k1[i][j],qcd_CMUL(y->v[1][n][i],qcd_CONJ(x->v[1][n][j])));
            k0[i][j]=qcd_CADD(k0[i][j],qcd_CMUL(x->v[0][n][i],qcd_CONJ(y->v[0][n][j])));
        }
        if(slot<count&&n+1==cutoffs[slot]) {
            memcpy(sources[slot],k0,sizeof(k0));memcpy(sinks[slot],k1,sizeof(k1));slot++;
        }
    }
    qcd_complex_16 pd[3][3],a[3][3],b[3][3],value[3][3];mat3_adjoint(pd,x->p);
    for(int sink=0;sink<count;sink++) {
        mat3_mul(a,y->p,sinks[sink]);mat3_mul(b,a,pd);
        for(int source=0;source<count;source++) {
            mat3_mul(value,b,sources[source]);out[source][sink]=mat3_trace(value);
        }
    }
    return 0;
}

/* Independent quadratic sum over the rectangular source/sink mode ranges. */
static qcd_complex_16 packet_mixed_direct(const endpoint_packet *x,const endpoint_packet *y,int source_modes,int sink_modes)
{
    qcd_complex_16 value={0};
    if(source_modes<1||sink_modes<1||source_modes>NV||sink_modes>NV)return (qcd_complex_16){NAN,NAN};
    for(int i=0;i<source_modes;i++)for(int j=0;j<sink_modes;j++)
        value=qcd_CADD(value,qcd_CMUL(packet_tau(y,i,j),qcd_CONJ(packet_tau(x,i,j))));
    return value;
}
#endif
