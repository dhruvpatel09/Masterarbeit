/* Fill the standard eigenvalue container using independently audited values.
   The controller copies the pilot container first; no vector data are rewritten. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "audited_values.h"

static int read_unsigned_scalar(hid_t f,const char *name,unsigned *out)
{
    hid_t d=H5Dopen2(f,name,H5P_DEFAULT),s=d<0?-1:H5Dget_space(d);
    int bad=s<0 || H5Sget_simple_extent_ndims(s)!=0;
    if(!bad && H5Dread(d,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,H5P_DEFAULT,out)<0)bad=1;
    if(s>=0)H5Sclose(s);
    if(d>=0)H5Dclose(d);
    return bad;
}

int main(int argc,char **argv)
{
    int provided,ranks;
    if(MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)!=MPI_SUCCESS)return 1;
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if(argc!=4 || ranks!=1 || provided<MPI_THREAD_FUNNELED) {
        fprintf(stderr,"usage: assemble_values EVAL_H5 VALUES_CSV PROVENANCE_SHA256; one MPI rank\n");
        MPI_Abort(MPI_COMM_WORLD,1);return 1;
    }
    static double values[192][100];
    int bad=read_audited_values(argv[2],values);
    if(strlen(argv[3])!=64)bad=1;
    if(bad) {fprintf(stderr,"Invalid audited eigenvalue table\n");MPI_Abort(MPI_COMM_WORLD,1);return 1;}
    hid_t f=H5Fopen(argv[1],H5F_ACC_RDWR,H5P_DEFAULT);
    if(f<0) {MPI_Abort(MPI_COMM_WORLD,1);return 1;}
    hid_t d=H5Dopen2(f,"disteigvals",H5P_DEFAULT),s=d<0?-1:H5Dget_space(d);
    hsize_t dims[2]={0,0};int times[192];unsigned nt=0,nv=0;
    if(s<0 || H5Sget_simple_extent_ndims(s)!=2 || H5Sget_simple_extent_dims(s,dims,NULL)<0 || dims[0]!=192 || dims[1]!=100)bad=1;
    hid_t td=H5Dopen2(f,"times",H5P_DEFAULT),ts=td<0?-1:H5Dget_space(td);hsize_t tdims[1]={0};
    if(ts<0 || H5Sget_simple_extent_ndims(ts)!=1 || H5Sget_simple_extent_dims(ts,tdims,NULL)<0 || tdims[0]!=192)bad=1;
    if(!bad && (H5Dread(td,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,times)<0 ||
        read_unsigned_scalar(f,"times_num",&nt) ||
        read_unsigned_scalar(f,"disteigpairs_num",&nv) || nt!=192 || nv!=100))bad=1;
    for(int t=0;t<192 && !bad;++t)if(times[t]!=t)bad=1;
    if(!bad && H5Dwrite(d,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values)<0)bad=1;
    if(!bad && (H5LTset_attribute_string(f,"/","qx_recovery_values_sha256",argv[3])<0 ||
        H5LTset_attribute_string(f,"/","qx_recovery_scope","192 times assembled from audited partial files; see recovery provenance manifest")<0))bad=1;
    if(ts>=0)H5Sclose(ts);if(td>=0)H5Dclose(td);if(s>=0)H5Sclose(s);if(d>=0)H5Dclose(d);
    if(H5Fflush(f,H5F_SCOPE_GLOBAL)<0)bad=1;
    if(H5Fclose(f)<0)bad=1;
    if(bad) {fprintf(stderr,"Eigenvalue assembly failed\n");MPI_Abort(MPI_COMM_WORLD,1);return 1;}
    printf("ASSEMBLED_VALUES times=192 Nv=100 rows=19200\n");fflush(stdout);
    MPI_Finalize();return 0;
}
