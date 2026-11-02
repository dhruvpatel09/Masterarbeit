/* First active_modes saved eigenvectors; I/O size remains NV. */
static int active_modes=NV;
/* Exact kernels extracted from the pinned qX profile source. */
static void mat3_zero(
    qcd_complex_16 a[3][3]
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j].re = 0.0;
            a[i][j].im = 0.0;
        }
    }
}

static void mat3_mul(
    qcd_complex_16 out[3][3],
    const qcd_complex_16 a[3][3],
    const qcd_complex_16 b[3][3]
) {
    qcd_complex_16 tmp[3][3];
    mat3_zero(tmp);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tmp[i][j] = qcd_CADD(
                    tmp[i][j],
                    qcd_CMUL(
                        a[i][k],
                        b[k][j]
                    )
                );
            }
        }
    }

    memcpy(out, tmp, sizeof(tmp));
}

static void mat3_adjoint(
    qcd_complex_16 out[3][3],
    const qcd_complex_16 in[3][3]
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = qcd_CONJ(
                in[j][i]
            );
        }
    }
}

static qcd_complex_16 mat3_trace(
    const qcd_complex_16 a[3][3]
) {
    qcd_complex_16 tr = {
        0.0,
        0.0
    };

    for (int i = 0; i < 3; i++) {
        tr = qcd_CADD(
            tr,
            a[i][i]
        );
    }

    return tr;
}
typedef struct {
    qcd_complex_16 v[2][NV][3];
    qcd_complex_16 p[3][3];
    int xyz[3];
} endpoint_packet;


static qcd_complex_16 packet_pair(const endpoint_packet *x,const endpoint_packet *y)
{
    qcd_complex_16 k1[3][3],k0[3][3],pd[3][3],a[3][3],b[3][3];
    mat3_zero(k1); mat3_zero(k0);
    for(int n=0;n<active_modes;++n) for(int i=0;i<3;++i) for(int j=0;j<3;++j) {
        k1[i][j]=qcd_CADD(k1[i][j],qcd_CMUL(y->v[1][n][i],qcd_CONJ(x->v[1][n][j])));
        k0[i][j]=qcd_CADD(k0[i][j],qcd_CMUL(x->v[0][n][i],qcd_CONJ(y->v[0][n][j])));
    }
    mat3_adjoint(pd,x->p);
    mat3_mul(a,y->p,k1); mat3_mul(b,a,pd); mat3_mul(a,b,k0);
    return mat3_trace(a);
}

static qcd_complex_16 packet_tau(const endpoint_packet *p,int i,int j)
{
    qcd_complex_16 tau={0};
    for(int a=0;a<3;++a) for(int b=0;b<3;++b)
        tau=qcd_CADD(tau,qcd_CMUL(qcd_CONJ(p->v[0][i][a]),qcd_CMUL(p->p[a][b],p->v[1][j][b])));
    return tau;
}

static qcd_complex_16 packet_direct(const endpoint_packet *x,const endpoint_packet *y)
{
    qcd_complex_16 sum={0};
    for(int i=0;i<active_modes;++i) for(int j=0;j<active_modes;++j)
        sum=qcd_CADD(sum,qcd_CMUL(packet_tau(y,i,j),qcd_CONJ(packet_tau(x,i,j))));
    return sum;
}
