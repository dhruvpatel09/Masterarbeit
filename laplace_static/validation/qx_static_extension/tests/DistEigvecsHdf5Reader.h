/* Type/interface shim for syntax checks only. Production uses the HDF5 reader. */
typedef struct {int unused;} DistEigvecsHdf5Reader;
typedef struct {double re,im;} DistEigvecsHdf5_ComplexDouble;
typedef struct {
    struct {int type;double theta[3];} bc;
    struct {char *type;unsigned int itercnt;double alpha;} laplacesmear;
} DistEigvecsHdf5_Metadata;
typedef struct {char *name,*infile;} DistEigvecsHdf5_Programdata;
int DistEigvecsHdf5Reader_Init(DistEigvecsHdf5Reader *);
int DistEigvecsHdf5Reader_OpenFile(DistEigvecsHdf5Reader *,MPI_Comm,const char *,char **,int *,int *,int *,int *,unsigned int *,DistEigvecsHdf5_Metadata *,DistEigvecsHdf5_Programdata *);
int DistEigvecsHdf5Reader_ReadDistEigvecsRange(DistEigvecsHdf5Reader *,int *,unsigned int,unsigned int,DistEigvecsHdf5_ComplexDouble *);
int DistEigvecsHdf5Reader_CloseFile(DistEigvecsHdf5Reader *);
void DistEigvecsHdf5Reader_Free(DistEigvecsHdf5Reader *);
