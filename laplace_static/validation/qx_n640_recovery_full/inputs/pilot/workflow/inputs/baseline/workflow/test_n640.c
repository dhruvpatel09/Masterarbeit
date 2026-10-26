/* Selected-slice qX Lanczos diagnostic. No eigenvector files are written. */
#include <mpi.h>
#include <qcd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "diagnostic_config.h"
#ifndef QX_VARIANT
#define QX_VARIANT "traced"
#endif

static int rank_id;

static void fail_if(int local_error, const char *message)
{
    int any;
    MPI_Allreduce(&local_error, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (any) {
        if (rank_id == 0) { fprintf(stderr, "QX_FATAL %s\n", message); fflush(stderr); }
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }
}

static int integer(const char *s, int *out)
{
    char *end;
    long x = strtol(s, &end, 10);
    if (!*s || *end || x < 0 || x > INT_MAX) return 1;
    *out = (int)x;
    return 0;
}

static qcd_complex_16 add(qcd_complex_16 a, qcd_complex_16 b)
{ return (qcd_complex_16){a.re+b.re, a.im+b.im}; }
static qcd_complex_16 mul(qcd_complex_16 a, qcd_complex_16 b)
{ return (qcd_complex_16){a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re}; }
static qcd_complex_16 conjugate(qcd_complex_16 a)
{ return (qcd_complex_16){a.re, -a.im}; }
static double square(qcd_complex_16 a) { return a.re*a.re+a.im*a.im; }

static void link_audit(qcd_gaugeField3d *u, const char *stage)
{
    unsigned long long count=0, bad=0, total, nonfinite;
    double unit=0., det=0., global_unit, global_det;
    for (long site=0; site<u->geo->lV3; ++site) for (int mu=1; mu<4; ++mu) {
        qcd_complex_16 (*v)[3] = u->D[site][mu];
        int invalid=0;
        for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            invalid |= !isfinite(v[i][j].re) || !isfinite(v[i][j].im);
        ++count;
        if (invalid) { ++bad; continue; }
        double norm=0.;
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
            qcd_complex_16 z={0.,0.};
            for (int k=0;k<3;++k) z=add(z,mul(conjugate(v[k][i]),v[k][j]));
            if (i==j) z.re-=1.;
            norm+=square(z);
        }
        unit=fmax(unit,sqrt(norm));
        qcd_complex_16 determinant={0.,0.};
        const int permutations[6][3]={{0,1,2},{1,2,0},{2,0,1},{2,1,0},{1,0,2},{0,2,1}};
        for (int p=0;p<6;++p) {
            qcd_complex_16 z=mul(mul(v[0][permutations[p][0]],v[1][permutations[p][1]]),v[2][permutations[p][2]]);
            if(p>=3) {z.re=-z.re;z.im=-z.im;}
            determinant=add(determinant,z);
        }
        determinant.re-=1.;
        det=fmax(det,sqrt(square(determinant)));
    }
    MPI_Allreduce(&count,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&bad,&nonfinite,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&unit,&global_unit,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&det,&global_det,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(rank_id==0) {
        printf("QX_LINKS stage=%s count=%llu nonfinite_links=%llu unitarity_max=%.17e determinant_max=%.17e\n",
               stage,total,nonfinite,global_unit,global_det); fflush(stdout);
    }
    fail_if(nonfinite!=0,"nonfinite spatial gauge links");
}

int main(int argc, char **argv)
{
    int provided, provided_min, ranks, cfg=0, time=0;
    if (MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)!=MPI_SUCCESS) return 1;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    MPI_Allreduce(&provided,&provided_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    fail_if(provided_min<MPI_THREAD_FUNNELED,"MPI thread support below FUNNELED");
    fail_if(argc!=4,"usage: test_n640 GAUGE_PATH CONFIGURATION GLOBAL_TIME");
    int argument_error=integer(argv[2],&cfg) || integer(argv[3],&time);
    fail_if(argument_error || time>=192 || ranks!=64,"invalid configuration/time or rank count");
    int L[4]={192,64,64,64}, P[4]={1,4,4,4};
    double theta[4]={0.,0.,0.,0.};
    qcd_geometry geo;
    fail_if(qcd_initGeometry(&geo,L,P,theta,rank_id,ranks),"geometry initialization");
    fail_if(qcd_initEO(&geo),"even/odd geometry initialization");
    int spatial_ranks;
    MPI_Comm_size(geo.comm3d,&spatial_ranks);
    fail_if(spatial_ranks!=64 || geo.Pos[0]!=0 || geo.root3d[time]!=0,"unexpected spatial communicator");
    qcd_gaugeField original, smeared, *input=&original, *ape=&smeared, *swap;
    qcd_gaugeField3d slice;
    int allocation_error=qcd_initGaugeField(&original,&geo);
    allocation_error+=qcd_initGaugeField(&smeared,&geo);
    allocation_error+=qcd_initGaugeField3d(&slice,&geo);
    fail_if(allocation_error,"gauge field allocation");
    qcd_spinorComponent3d *vectors=calloc(QX_NV,sizeof(*vectors)), applied, residual;
    double *eval=malloc(QX_NV*sizeof(*eval));
    fail_if(vectors==NULL || eval==NULL,"eigenpair allocation");
    allocation_error=0;
    for(int i=0;i<QX_NV;++i) {allocation_error+=qcd_initSpinorComponent3d(&vectors[i],&geo);eval[i]=NAN;}
    allocation_error+=qcd_initSpinorComponent3d(&applied,&geo);
    allocation_error+=qcd_initSpinorComponent3d(&residual,&geo);
    fail_if(allocation_error,"spinor allocation");
    qcd_lanczosParams par;
    qcd_defaultLanczosParams(&par);
    par.res=QX_RES;par.redTol=QX_RED_TOL;par.deg=QX_DEG;par.upper=QX_UPPER;
    par.tolReortho=QX_REORTHO;par.maxRestarts=QX_MAX_RESTARTS;
    par.nLanczosFactor=QX_LANCZOS_FACTOR;par.redFreqFactor=QX_RED_FREQ_FACTOR;
    if(rank_id==0) {
        printf("QX_CASE configuration=%d time=%d variant=%s Nv=%d ranks=%d spatial_ranks=%d P=1,4,4,4 res=%.17e redTol=%.17e deg=%d upper=%.17e tolReortho=%.17e maxRestarts=%d nLanczosFactor=%d redFreqFactor=%d nAPE=%d alphaAPE=%.17e mpi_thread_min=%d\n",
               cfg,time,QX_VARIANT,QX_NV,ranks,spatial_ranks,par.res,par.redTol,par.deg,par.upper,
               par.tolReortho,par.maxRestarts,par.nLanczosFactor,par.redFreqFactor,QX_NAPE,QX_ALPHA_APE,provided_min);
        fflush(stdout);
    }
    double start=MPI_Wtime();
    fail_if(qcd_getGaugeField(argv[1],qcd_GF_OPENQCD,&original),"gauge field read");
    double plaquette=qcd_calculatePlaquette(&original);
    if(rank_id==0) {printf("QX_GAUGE plaquette=%.17e load_seconds=%.9f\n",plaquette,MPI_Wtime()-start);fflush(stdout);}
    fail_if(!isfinite(plaquette),"nonfinite plaquette");
    qcd_copyGaugeField3dGaugeField(&slice,&original,time);
    link_audit(&slice,"raw");
    start=MPI_Wtime();
    for(int step=0;step<QX_NAPE;++step) {
        qcd_apeSmear3d(ape,input,QX_ALPHA_APE);
        if(step!=QX_NAPE-1) {swap=input;input=ape;ape=swap;}
    }
    if(rank_id==0) {printf("QX_APE steps=%d alpha=%.17e seconds=%.9f\n",QX_NAPE,QX_ALPHA_APE,MPI_Wtime()-start);fflush(stdout);}
    qcd_copyGaugeField3dGaugeField(&slice,ape,time);
    link_audit(&slice,"ape");
    qcd_lanczosStatus status;
    memset(&status,0,sizeof(status));
    start=MPI_Wtime();
    int code=qcd_lanczosCheby(eval,vectors,&status,&slice,QX_NV,&par);
    double elapsed=MPI_Wtime()-start;
    int code_min,code_max;
    MPI_Allreduce(&code,&code_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&code,&code_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(rank_id==0) {
        printf("QX_RETURN code=%d code_min=%d code_max=%d restarts=%d nRed=%d maxRed=%d maxRes=%.17e usedLapack=%d seconds=%.9f\n",
               code,code_min,code_max,status.restarts,status.nRed,status.maxRed,status.maxRes,status.usedLapack,elapsed);fflush(stdout);
    }
    fail_if(code_min!=code_max,"rank disagreement in solver return code");
    if(code==0 || code==2 || code==5) {
        double max_res=0.,max_norm=0.,max_gram=0.;
        int invalid=0;
        for(int i=0;i<QX_NV;++i) {
            double norm=qcd_normSpinorComponent3d(&vectors[i]);
            qcd_applyCovariantLaplaceSpinorComponent3d(&applied,&vectors[i],&slice,0.);
            qcd_raxpySpinorComponent3d(&residual,-eval[i],&vectors[i],&applied);
            double res_norm=qcd_normSpinorComponent3d(&residual);
            double relative=res_norm/(fabs(eval[i])*norm);
            int finite=isfinite(norm) && isfinite(eval[i]) && isfinite(relative) && norm>0. && eval[i]>0.;
            invalid+=!finite;
            max_res=fmax(max_res,relative);max_norm=fmax(max_norm,fabs(norm-1.));
            if(rank_id==0) printf("QX_AUDIT_MODE index=%d lambda=%.17e norm=%.17e relative_residual=%.17e finite_positive=%d\n",i,eval[i],norm,relative,finite);
            for(int j=0;j<i;++j) {
                qcd_complex_16 overlap=qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&vectors[j],&vectors[i]);
                double magnitude=hypot(overlap.re,overlap.im);
                invalid+=!isfinite(magnitude);max_gram=fmax(max_gram,magnitude);
            }
        }
        int valid=(code==0 && invalid==0 && max_res<=QX_RES && max_norm<=1e-10 && max_gram<=1e-8);
        if(rank_id==0) {
            printf("QX_AUDIT modes=%d invalid=%d max_relative_residual=%.17e max_norm_defect=%.17e max_gram_offdiag=%.17e valid_selected_slice=%d\n",
                   QX_NV,invalid,max_res,max_norm,max_gram,valid);fflush(stdout);
        }
    } else if(rank_id==0) {printf("QX_AUDIT skipped=1 reason=solver_return_code\n");fflush(stdout);}
    for(int i=0;i<QX_NV;++i) qcd_destroySpinorComponent3d(&vectors[i]);
    free(vectors);free(eval);
    qcd_destroySpinorComponent3d(&applied);qcd_destroySpinorComponent3d(&residual);
    qcd_destroyGaugeField3d(&slice);qcd_destroyGaugeField(&original);qcd_destroyGaugeField(&smeared);
    qcd_destroyGeometry(&geo);
    if(rank_id==0) {printf("QX_COMPLETE configuration=%d time=%d solver_code=%d\n",cfg,time,code);fflush(stdout);}
    MPI_Finalize();
    /* A reported solver failure is a completed diagnostic, not accepted production data. */
    return 0;
}
