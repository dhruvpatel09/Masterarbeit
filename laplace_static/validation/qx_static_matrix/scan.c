/* Full source/sink cutoff matrix using the validated 100-mode endpoint payload. */
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <qcd.h>
#include "solver.h"
#include <DistEigvecsHdf5Reader.h>
#define NV 100
#define NT 192
#define NS 64
#define PATH_SIZE 4096
#include "kernels.h"
#include "matrix_plan.h"
#include "mixed_contractions.h"
#include "io_helpers.h"
#include "mpi_helpers.h"
#include "su3_window.h"

static qcd_spinorComponent3d source_modes[NV],sink_modes[NV];
static int cached_source=-1,records=0,cfg,piter,myid,numprocs;
static const char *evroot;

static void fail(int bad,const char *what) {
    int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(any){if(!myid)fprintf(stderr,"FAIL %s\n",what);MPI_Abort(MPI_COMM_WORLD,1);}
}

static void modes(int t,qcd_geometry *g,qcd_spinorComponent3d *v) {
    char path[PATH_SIZE];int n=snprintf(path,sizeof(path),"%s/qX_p2n%d_evec_t%d.h5",evroot,cfg,t);
    fail(n<0||n>=(int)sizeof(path),"eigenvector path overflow");
    read_evec_time(path,cfg,t,"qX_p2",g,v);int bad=0;
    for(int j=0;j<NV;j++)for(long x=0;x<g->lV3;x++)for(int c=0;c<3;c++)
        if(!isfinite(v[j].D[x][c].re)||!isfinite(v[j].D[x][c].im))bad=1;
    fail(bad,"nonfinite mode");
}

static void emit(int tau,int axis,int r,int source,int sink,const double value[2],unsigned long long count) {
    if(!myid){printf("MATRIX_RT cfg=%d piter=%d scheme=fixed t0=90 tau=%d R=%d axis=%d Nv_src=%d Nv_sink=%d Nsrc=%llu L_RE=%.17e L_IM=%.17e\n",
        cfg,piter,tau,r,axis,source,sink,count,value[0]/count,value[1]/count);records++;}
}

static void direct_check(const char *prefix,int tau,int r,int axis,long site,int source,int sink,qcd_complex_16 fast,qcd_complex_16 direct) {
    double error=hypot(fast.re-direct.re,fast.im-direct.im);
    double tolerance=1e-20+5e-12*fmax(hypot(fast.re,fast.im),hypot(direct.re,direct.im));
    int good=isfinite(error)&&isfinite(tolerance)&&error<=tolerance;
    printf("%s cfg=%d piter=%d tau=%d R=%d axis=%d rank=0 site=%ld Nv_src=%d Nv_sink=%d fast_re=%.17e fast_im=%.17e direct_re=%.17e direct_im=%.17e abs_error=%.17e tolerance=%.17e status=%s\n",
        prefix,cfg,piter,tau,r,axis,site,source,sink,fast.re,fast.im,direct.re,direct.im,error,tolerance,good?"PASS":"FAIL");
    if(!good)MPI_Abort(MPI_COMM_WORLD,1);
}

static void measure(const qcd_gaugeField *u,int tau) {
    qcd_geometry *g=u->geo;int t0=MATRIX_T0;
    if(cached_source!=t0){modes(t0,g,source_modes);cached_source=t0;}
    if(tau)modes(t0+tau,g,sink_modes);
    endpoint_packet *base=calloc((size_t)g->lV3,sizeof(*base));
    endpoint_packet *a=malloc((size_t)g->lV3*sizeof(*a)),*b=malloc((size_t)g->lV3*sizeof(*b));
    endpoint_packet *send=malloc(256*sizeof(*send)),*recv=malloc(256*sizeof(*recv));
    fail(!base||!a||!b||!send||!recv,"endpoint allocation");
    for(int x=0;x<16;x++)for(int y=0;y<16;y++)for(int z=0;z<16;z++) {
        long site=qcd_LEXIC0(x,y,z,g->lL);
        for(int j=0;j<NV;j++)for(int c=0;c<3;c++) {
            base[site].v[0][j][c]=source_modes[j].D[site][c];
            base[site].v[1][j][c]=(tau?sink_modes:source_modes)[j].D[site][c];
        }
        temporal_product(u,x,y,z,t0,t0+tau,base[site].p);
        base[site].xyz[0]=x+16*g->Pos[1];base[site].xyz[1]=y+16*g->Pos[2];base[site].xyz[2]=z+16*g->Pos[3];
    }
    int minus[3],plus[3];rank_neighbours(g,minus,plus);
    double combined[33][4][4][2]={{{{0}}}};
    for(int axis=0;axis<3;axis++) {
        memcpy(a,base,(size_t)g->lV3*sizeof(*a));
        for(int r=1;r<=32;r++) {
            shift_endpoint(a,b,send,recv,g,axis,minus[axis],plus[axis]);endpoint_packet *tmp=a;a=b;b=tmp;
            verify_shift(base,a,g->lV3,axis,r);
            if(!matrix_contains(matrix_radii,4,r))continue;
            double sum[4][4][2]={{{0}}},correction[4][4][2]={{{0}}},global[4][4][2];int bad=0;
            for(long site=0;site<g->lV3;site++) {
                qcd_complex_16 values[4][4];
                if(packet_mixed_all(&base[site],&a[site],matrix_cutoffs,4,values)){bad=1;continue;}
                for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
                    if(!isfinite(values[i][j].re)||!isfinite(values[i][j].im))bad=1;
                    kahan_add(values[i][j].re,&sum[i][j][0],&correction[i][j][0]);
                    kahan_add(values[i][j].im,&sum[i][j][1],&correction[i][j][1]);
                }
            }
            fail(bad,"nonfinite mixed correlator");
            if(!myid&&matrix_contains(matrix_direct_taus,4,tau)) {
                for(int sample=0;sample<2;sample++) {
                    long site=sample?g->lV3-1:0;qcd_complex_16 values[4][4];
                    if(packet_mixed_all(&base[site],&a[site],matrix_cutoffs,4,values))MPI_Abort(MPI_COMM_WORLD,1);
                    for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
                        qcd_complex_16 direct=packet_mixed_direct(&base[site],&a[site],matrix_cutoffs[i],matrix_cutoffs[j]);
                        direct_check("MATRIX_DIRECT",tau,r,axis,site,matrix_cutoffs[i],matrix_cutoffs[j],values[i][j],direct);
                        if(i==j){active_modes=matrix_cutoffs[i];direct=packet_pair(&base[site],&a[site]);
                            direct_check("MATRIX_LEGACY",tau,r,axis,site,matrix_cutoffs[i],matrix_cutoffs[j],values[i][j],direct);}
                    }
                    active_modes=NV;
                }
            }
            MPI_Allreduce(sum,global,32,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            unsigned long long local_count=(unsigned long long)g->lV3,total_count;
            MPI_Allreduce(&local_count,&total_count,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
            fail(total_count!=64ULL*64*64,"source count");bad=0;
            for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
                if(r==32&&fabs(global[i][j][1]/total_count)>1e-18+1e-10*fabs(global[i][j][0]/total_count))bad=1;
                emit(tau,axis,r,matrix_cutoffs[i],matrix_cutoffs[j],global[i][j],total_count);
                for(int k=0;k<2;k++)combined[r][i][j][k]+=global[i][j][k];
            }
            fail(bad,"half-box imaginary cancellation");
            if(matrix_contains(matrix_block_taus,2,tau)&&matrix_contains(matrix_block_radii,2,r)) {
                double sums[64][4][4][2];int positions[64][3];
                MPI_Gather(sum,32,MPI_DOUBLE,sums,32,MPI_DOUBLE,0,MPI_COMM_WORLD);
                MPI_Gather(&g->Pos[1],3,MPI_INT,positions,3,MPI_INT,0,MPI_COMM_WORLD);
                if(!myid)for(int rank=0;rank<64;rank++)for(int i=0;i<4;i++)for(int j=0;j<4;j++)
                    printf("MATRIX_BLOCK cfg=%d piter=%d t0=90 tau=%d R=%d axis=%d Nv_src=%d Nv_sink=%d rank=%d px=%d py=%d pz=%d Nsrc=4096 L_RE=%.17e L_IM=%.17e\n",
                        cfg,piter,tau,r,axis,matrix_cutoffs[i],matrix_cutoffs[j],rank,positions[rank][0],positions[rank][1],positions[rank][2],sums[rank][i][j][0]/4096.,sums[rank][i][j][1]/4096.);
            }
        }
        if(tau==0) {
            for(int r=33;r<=64;r++){shift_endpoint(a,b,send,recv,g,axis,minus[axis],plus[axis]);endpoint_packet *tmp=a;a=b;b=tmp;}
            fail(memcmp(base,a,(size_t)g->lV3*sizeof(*base))!=0,"full periodic payload roundtrip");
            if(!myid)printf("PERIODIC_ROUNDTRIP axis=%d cells=64 status=PASS\n",axis);
        }
    }
    for(int r=0;r<4;r++)for(int i=0;i<4;i++)for(int j=0;j<4;j++)
        emit(tau,-1,matrix_radii[r],matrix_cutoffs[i],matrix_cutoffs[j],combined[matrix_radii[r]][i][j],3ULL*64*64*64);
    free(base);free(a);free(b);free(send);free(recv);
}

int main(int argc,char **argv) {
    solver_mpiInit(argc,argv);setvbuf(stdout,NULL,_IONBF,0);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    cfg=read_env_int("QX_CFG",650,100,800);piter=read_env_int("HYP_PITER",10,1,1000);
    fail(numprocs!=64||cfg%10||(piter!=10&&piter!=40)||solver_getNThread()!=1,"64 ranks, one thread, valid configuration required");
    int provided;MPI_Query_thread(&provided);fail(provided<MPI_THREAD_FUNNELED,"MPI FUNNELED support");
    evroot=env_or_default("EVROOT","");fail(!*evroot,"EVROOT is required");
    int L[4]={192,64,64,64},P[4]={1,4,4,4};double theta[4]={0};qcd_geometry geo;
    fail(qcd_initGeometry(&geo,L,P,theta,myid,numprocs),"geometry");fail(qcd_initEO(&geo),"even/odd geometry");
    qcd_gaugeField original,transport;
    fail(qcd_initGaugeField(&original,&geo)||qcd_initGaugeField(&transport,&geo),"gauge allocation");
    char gauge[PATH_SIZE];int n=snprintf(gauge,sizeof(gauge),"%s/qX_p2n%d",env_or_default("QX_ROOT","/hppfs/scratch/00/di98qes3/qX"),cfg);
    fail(n<0||n>=(int)sizeof(gauge),"gauge path");fail(qcd_getGaugeField(gauge,qcd_GF_OPENQCD,&original),"gauge read");
    solver_parameters par={0};par.kappa=.125;par.bc=0;par.cF=par.cF_prime=1;
    fail(solver_init(&geo,&par)||solver_updateGaugeField(&original),"solver conversion");
    fail(solver_hypTimeLinks(&transport,1.,1.,.5,piter),"temporal HYP");
    double local[2]={0},global[2];int bad=0;
    for(int t=0;t<192;t++)for(int x=0;x<16;x++)for(int y=0;y<16;y++)for(int z=0;z<16;z++) {
        long v=qcd_LEXIC(t,x,y,z,geo.lL);
        for(int mu=0;mu<4;mu++)for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
            double dr=fabs(original.D[v][mu][i][j].re-transport.D[v][mu][i][j].re);
            double di=fabs(original.D[v][mu][i][j].im-transport.D[v][mu][i][j].im);
            if(!isfinite(dr)||!isfinite(di))bad=1;
            int k=(mu==0&&t<191)?0:1;local[k]=fmax(local[k],fmax(dr,di));
        }
    }
    fail(bad,"nonfinite HYP field");MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    fail(global[0]<=1e-13||global[1]>1e-13,"temporal HYP isolation");
    if(!myid)printf("MATRIX_BEGIN cfg=%d piter=%d saved_Nv=100 source_cutoffs=40,60,80,100 sink_cutoffs=40,60,80,100 ranks=64 grid=1x4x4x4 t0=90 tau_min=0 tau_max=32 radii=4,12,24,32 temporal_bc=open hyp_steps=1 temporal_change=%.17e preserved_change=%.17e\n",cfg,piter,global[0],global[1]);
    profile_window_su3(&transport,cfg,piter,90,122,myid);
    for(int i=0;i<NV;i++)fail(qcd_initSpinorComponent3d(&source_modes[i],&geo)||qcd_initSpinorComponent3d(&sink_modes[i],&geo),"mode allocation");
    for(int tau=0;tau<=MATRIX_TMAX;tau++)measure(&transport,tau);
    if(!myid)printf("MATRIX_END cfg=%d piter=%d records=%d status=%s\n",cfg,piter,records,records==MATRIX_RECORDS?"PASS":"FAIL");
    fail(!myid&&records!=MATRIX_RECORDS,"matrix record coverage");
    for(int i=0;i<NV;i++){qcd_destroySpinorComponent3d(&source_modes[i]);qcd_destroySpinorComponent3d(&sink_modes[i]);}
    qcd_destroyGaugeField(&transport);qcd_destroyGaugeField(&original);qcd_destroyEO(&geo);qcd_destroyGeometry(&geo);
    MPI_Finalize();return 0;
}
