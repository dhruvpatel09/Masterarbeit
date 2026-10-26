#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "qcd.h"
#include "selection.h"
#include "DistEigvecsHdf5Reader.h"
#include "metrics.h"

#define NV 100
#define NT 192
static int rank_id;
static void require(int bad, const char *message, MPI_Comm comm)
{
    int any, group_rank;
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);
    if (any) {
        MPI_Comm_rank(comm,&group_rank);
        if (!group_rank) fprintf(stderr,"FAIL rank=%d: %s\n",rank_id,message);
        MPI_Abort(MPI_COMM_WORLD,1);
        exit(1);
    }
}
/* Shape-checked scalar/array read; a malformed HDF5 shape never overfills a buffer. */
static int read_numeric(hid_t f, const char *name, hid_t type, int nd,
                         const hsize_t *dims, void *out)
{
    hid_t d=H5Dopen2(f,name,H5P_DEFAULT);
    if (d<0) return 1;
    hid_t s=H5Dget_space(d);
    if (s<0) {H5Dclose(d);return 1;}
    int bad=H5Sget_simple_extent_ndims(s)!=nd;
    hsize_t got[2]={0,0};
    if (!bad && nd) {
        if (H5Sget_simple_extent_dims(s,got,NULL)<0) bad=1;
        for (int k=0;k<nd;k++) if (got[k]!=dims[k]) bad=1;
    }
    if (!bad && H5Dread(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,out)<0) bad=1;
    H5Sclose(s); H5Dclose(d);return bad;
}
static int read_string(hid_t f,const char *name,char out[128])
{
    hid_t d=H5Dopen2(f,name,H5P_DEFAULT);
    if(d<0)return 1;
    hid_t s=H5Dget_space(d),ty=H5Dget_type(d);
    int bad=s<0 || ty<0;
    if(!bad) bad=H5Sget_simple_extent_ndims(s)!=0 ||
        H5Tget_class(ty)!=H5T_STRING || H5Tis_variable_str(ty)!=0 ||
        H5Tget_size(ty)>=128;
    memset(out,0,128);
    if(!bad && H5Dread(d,ty,H5S_ALL,H5S_ALL,H5P_DEFAULT,out)<0)bad=1;
    if(s>=0)H5Sclose(s);
    if(ty>=0)H5Tclose(ty);
    H5Dclose(d);return bad;
}
/* The eigenvalue table is small: root reads all rows, checks metadata and
   coverage, then broadcasts it. Respect stored time ordering. */
static int load_values(const char *file,int cfg,double values[NT][NV])
{
    hid_t f=H5Fopen(file,H5F_ACC_RDONLY,H5P_DEFAULT);
    if(f<0)return 1;
    int nc=0,bc=0,times[NT],seen[NT]={0};
    unsigned nt=0,nv=0,ape=0;
    double alpha=0,theta[3]={0},raw[NT][NV];
    char name[128],smear[128];
    hsize_t one[1]={NT},three[1]={3},two[2]={NT,NV};
    int bad=read_numeric(f,"cnfg_id",H5T_NATIVE_INT,0,NULL,&nc) ||
        read_numeric(f,"times_num",H5T_NATIVE_UINT,0,NULL,&nt) ||
        read_numeric(f,"disteigpairs_num",H5T_NATIVE_UINT,0,NULL,&nv) ||
        read_numeric(f,"times",H5T_NATIVE_INT,1,one,times) ||
        read_numeric(f,"disteigvals",H5T_NATIVE_DOUBLE,2,two,raw) ||
        read_numeric(f,"metadata/bc_params/type",H5T_NATIVE_INT,0,NULL,&bc) ||
        read_numeric(f,"metadata/bc_params/theta",H5T_NATIVE_DOUBLE,1,three,theta) ||
        read_numeric(f,"metadata/laplace_smear_params/itercnt",H5T_NATIVE_UINT,0,NULL,&ape) ||
        read_numeric(f,"metadata/laplace_smear_params/alpha",H5T_NATIVE_DOUBLE,0,NULL,&alpha) ||
        read_string(f,"run_name",name) ||
        read_string(f,"metadata/laplace_smear_params/type",smear);
    if(!bad) bad=nc!=cfg || nt!=NT || nv!=NV || bc!=0 || ape!=20 ||
        !isfinite(alpha) || fabs(alpha-0.5)>1e-14 || strcmp(name,"qX_p2") || strcmp(smear,"APE");
    for(int k=0;!bad && k<3;k++) if(!isfinite(theta[k]) || fabs(theta[k])>1e-14)bad=1;
    for(int row=0;!bad && row<NT;row++) {
        int t=times[row];
        if(t<0 || t>=NT || seen[t]) {bad=1;break;}
        seen[t]=1;
        if(!qx_selected(t))continue; /* Unwritten rows are not valid data. */
        for(int i=0;i<NV;i++) {
            double x=raw[row][i];
            if(!isfinite(x) || x<=0 || x>12.0+1e-8) {bad=1;break;}
            if(i && x+1e-12<raw[row][i-1]) {bad=1;break;}
            values[t][i]=x;
        }
    }
    H5Fclose(f);return bad;
}
static int selected(int t)
{
    for(size_t k=0;k<sizeof(selected_times)/sizeof(selected_times[0]);k++)
        if(t==selected_times[k])return 1;
    return 0;
}
int main(int argc,char **argv)
{
    int requested;
    if(MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&requested)!=MPI_SUCCESS)return 1;
    int ranks,provided,minprovided;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    MPI_Query_thread(&provided);
    MPI_Allreduce(&provided,&minprovided,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    require(argc!=5 || ranks!=64 || minprovided<MPI_THREAD_FUNNELED,
            "need 64 ranks, FUNNELED; CFG GAUGE EV_DIR EVAL",MPI_COMM_WORLD);
    char *end;long nc=strtol(argv[1],&end,10);
    require(*end || nc!=640,"invalid configuration",MPI_COMM_WORLD);
    int cfg=(int)nc,L[4]={NT,64,64,64},P[4]={1,4,4,4};
    double theta[4]={0,0,0,0};
    qcd_geometry geo;
    require(qcd_initGeometry(&geo,L,P,theta,rank_id,ranks),"geometry",MPI_COMM_WORLD);
    require(qcd_initEO(&geo),"even/odd geometry",MPI_COMM_WORLD);
    qcd_gaugeField u[2];qcd_gaugeField3d u3;
    qcd_spinorComponent3d v[NV],av;
    require(qcd_initGaugeField(&u[0],&geo),"gauge allocation 0",MPI_COMM_WORLD);
    require(qcd_initGaugeField(&u[1],&geo),"gauge allocation 1",MPI_COMM_WORLD);
    require(qcd_initGaugeField3d(&u3,&geo),"spatial gauge allocation",MPI_COMM_WORLD);
    for(int i=0;i<NV;i++) require(qcd_initSpinorComponent3d(&v[i],&geo),"mode allocation",MPI_COMM_WORLD);
    require(qcd_initSpinorComponent3d(&av,&geo),"operator-output allocation",MPI_COMM_WORLD);
    size_t components=(size_t)geo.lV3*3;
    DistEigvecsHdf5_ComplexDouble *buffer=malloc(NV*components*sizeof(*buffer));
    double *gram_local=calloc(2*NV*NV,sizeof(double)),*gram_global=calloc(2*NV*NV,sizeof(double));
    require(!buffer || !gram_local || !gram_global,"read/Gram buffer allocation",MPI_COMM_WORLD);
    static double values[NT][NV];int value_bad=0;
    if(!rank_id)value_bad=load_values(argv[4],cfg,values);
    MPI_Bcast(&value_bad,1,MPI_INT,0,MPI_COMM_WORLD);
    require(value_bad,"eigenvalue metadata/shape/coverage/finite/order check",MPI_COMM_WORLD);
    MPI_Bcast(values,NT*NV,MPI_DOUBLE,0,MPI_COMM_WORLD);
    require(qcd_getGaugeField(argv[2],qcd_GF_OPENQCD,&u[0]),"read gauge",MPI_COMM_WORLD);
    qcd_gaugeField *input=&u[0],*output=&u[1];
    double began=MPI_Wtime();
    for(int step=0;step<20;step++) {
        require(qcd_apeSmear3d(output,input,0.5),"spatial APE20",MPI_COMM_WORLD);
        if(step!=19) {qcd_gaugeField *tmp=input;input=output;output=tmp;}
    }
    if(!rank_id) {
        printf("QUALITY_START cfg=%d grid=1x4x4x4 Nv=100 APE=20,0.5 FUNNELED=%d smear_seconds=%.17e\n",cfg,minprovided,MPI_Wtime()-began);
        fflush(stdout);
    }
    DistEigvecsHdf5Reader reader={0};
    require(DistEigvecsHdf5Reader_Init(&reader),"reader init",MPI_COMM_WORLD);
    int spatialP[3]={4,4,4},group_rank,local_fail=0,local_slices=0;
    MPI_Comm_rank(geo.comm3d,&group_rank);
    for(int lt=0;lt<geo.lL[0];lt++) {
        int t=geo.Pos[0]*geo.lL[0]+lt;if(!selected(t))continue;
        char path[4096],*run_name=NULL;int ec=0,et=-1;unsigned count=0;
        DistEigvecsHdf5_Metadata m={0};DistEigvecsHdf5_Programdata prog={0};
        int wrote=snprintf(path,sizeof(path),"%s/qX_p2n%d_evec_t%d.h5",argv[3],cfg,t);
        require(wrote<0 || wrote>=(int)sizeof(path),"eigenvector path",geo.comm3d);
        require(DistEigvecsHdf5Reader_OpenFile(&reader,geo.comm3d,path,&run_name,&ec,
            &geo.L[1],spatialP,&et,&count,&m,&prog),"open eigenvector file",geo.comm3d);
        int meta_bad=!run_name || ec!=cfg || et!=t || count!=NV || m.bc.type!=0 ||
            !m.laplacesmear.type || strcmp(m.laplacesmear.type,"APE") || m.laplacesmear.itercnt!=20 ||
            !isfinite(m.laplacesmear.alpha) || fabs(m.laplacesmear.alpha-0.5)>1e-14;
        if(run_name && strcmp(run_name,"qX_p2"))meta_bad=1;
        for(int k=0;k<3;k++)if(!isfinite(m.bc.theta[k]) || fabs(m.bc.theta[k])>1e-14)meta_bad=1;
        require(meta_bad,"eigenvector metadata",geo.comm3d);
        require(DistEigvecsHdf5Reader_ReadDistEigvecsRange(&reader,&geo.Pos[1],0,NV,buffer),"read all 100 modes",geo.comm3d);
        require(DistEigvecsHdf5Reader_CloseFile(&reader),"close eigenvector file",geo.comm3d);
        free(run_name);free(m.laplacesmear.type);free(prog.name);free(prog.infile);
        int nonfinite=0;
        for(int i=0;i<NV;i++) {
            for(size_t x=0;x<(size_t)geo.lV3;x++)for(int c=0;c<3;c++) {
                size_t k=(size_t)i*components+x*3+(size_t)c;
                v[i].D[x][c].re=buffer[k].re;v[i].D[x][c].im=buffer[k].im;
                if(!isfinite(buffer[k].re) || !isfinite(buffer[k].im))nonfinite=1;
            }
            v[i].BPup2date=0;v[i].BMup2date=0;
        }
        require(nonfinite,"non-finite saved mode",geo.comm3d);
        qcd_copyGaugeField3dGaugeField(&u3,output,t);
        require(qcd_communicateGauge3dP(&u3),"spatial gauge halo",geo.comm3d);
        qcd_waitall(&geo);
        double maxrel=0,maxnorm=0,maxdiag=0,maxoff=0,maxray_im=0;
        for(int i=0;i<NV;i++) {
            require(qcd_applyCovariantLaplaceSpinorComponent3d(&av,&v[i],&u3,0.0),"Laplacian application",geo.comm3d);
            double sums[4],total[4];
            mode_sums((double*)&v[i].D[0][0],(double*)&av.D[0][0],components,values[t][i],sums);
            MPI_Allreduce(sums,total,4,MPI_DOUBLE,MPI_SUM,geo.comm3d);
            require(!isfinite(total[0]) || !isfinite(total[1]) || !isfinite(total[2]) ||
                    !isfinite(total[3]) || total[0]<=0 || total[1]<0,"invalid mode reductions",geo.comm3d);
            double norm=sqrt(total[0]),res=sqrt(total[1]);
            double rel=res/(fabs(values[t][i])*norm);
            double rayre=total[2]/total[0],rayim=total[3]/total[0];
            maxrel=fmax(maxrel,rel);maxnorm=fmax(maxnorm,fabs(norm-1.0));
            maxray_im=fmax(maxray_im,fabs(rayim));
            if(rel>1e-9 || fabs(norm-1.0)>1e-9)local_fail=1;
            if(!group_rank)printf("MODE cfg=%d t=%d i=%d lambda=%.17e norm=%.17e residual_l2=%.17e residual_relative=%.17e rayleigh_re=%.17e rayleigh_im=%.17e\n",cfg,t,i,values[t][i],norm,res,rel,rayre,rayim);
        }
        memset(gram_local,0,2*NV*NV*sizeof(double));
        for(int i=0;i<NV;i++)for(int j=i;j<NV;j++)
            overlap_sums((double*)&v[i].D[0][0],(double*)&v[j].D[0][0],components,&gram_local[2*(i*NV+j)]);
        MPI_Allreduce(gram_local,gram_global,2*NV*NV,MPI_DOUBLE,MPI_SUM,geo.comm3d);
        int gram_bad=0;
        for(int i=0;i<NV;i++)for(int j=i;j<NV;j++) {
            double re=gram_global[2*(i*NV+j)],im=gram_global[2*(i*NV+j)+1];
            if(!isfinite(re) || !isfinite(im))gram_bad=1;
            if(i==j)maxdiag=fmax(maxdiag,hypot(re-1.0,im));else maxoff=fmax(maxoff,hypot(re,im));
        }
        require(gram_bad,"non-finite Gram entry",geo.comm3d);
        if(maxdiag>1e-9 || maxoff>1e-9 || maxray_im>1e-10)local_fail=1;
        if(!group_rank) {
            printf("SLICE cfg=%d t=%d Nv=100 max_relative_residual=%.17e max_norm_defect=%.17e gram_diag_defect=%.17e gram_offdiag=%.17e max_rayleigh_im=%.17e lambda_min=%.17e lambda_max=%.17e status=%s\n",cfg,t,maxrel,maxnorm,maxdiag,maxoff,maxray_im,values[t][0],values[t][NV-1],
                (maxrel<=1e-9 && maxnorm<=1e-9 && maxdiag<=1e-9 && maxoff<=1e-9 && maxray_im<=1e-10)?"PASS":"FAIL");
            fflush(stdout);local_slices++;
        }
    }
    int anyfail=0,slices=0;
    MPI_Allreduce(&local_fail,&anyfail,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&local_slices,&slices,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    require(slices!=(int)(sizeof(selected_times)/sizeof(selected_times[0])),"selected slice coverage",MPI_COMM_WORLD);
    if(!rank_id) {printf("QUALITY_END cfg=%d slices=%d status=%s\n",cfg,slices,anyfail?"FAIL":"PASS");fflush(stdout);}
    DistEigvecsHdf5Reader_Free(&reader);
    for(int i=0;i<NV;i++)qcd_destroySpinorComponent3d(&v[i]);
    qcd_destroySpinorComponent3d(&av);qcd_destroyGaugeField3d(&u3);
    qcd_destroyGaugeField(&u[1]);qcd_destroyGaugeField(&u[0]);
    free(buffer);free(gram_local);free(gram_global);
    qcd_destroyEO(&geo);qcd_destroyGeometry(&geo);MPI_Finalize();return anyfail?1:0;
}

