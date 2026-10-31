/* Declarations for local test builds only; never used by the production build. */
#ifndef QX_STATIC_TEST_MPI_H
#define QX_STATIC_TEST_MPI_H
typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Request;
typedef int MPI_Status;
#define MPI_COMM_WORLD 0
#define MPI_INT 1
#define MPI_DOUBLE 2
#define MPI_UNSIGNED_LONG_LONG 3
#define MPI_BYTE 4
#define MPI_MAX 1
#define MPI_MIN 2
#define MPI_SUM 3
#define MPI_THREAD_FUNNELED 1
#define MPI_STATUS_IGNORE ((MPI_Status *)0)
int MPI_Comm_size(MPI_Comm,int *);
int MPI_Comm_rank(MPI_Comm,int *);
int MPI_Allgather(const void *,int,MPI_Datatype,void *,int,MPI_Datatype,MPI_Comm);
int MPI_Allreduce(const void *,void *,int,MPI_Datatype,int,MPI_Comm);
int MPI_Gather(const void *,int,MPI_Datatype,void *,int,MPI_Datatype,int,MPI_Comm);
int MPI_Sendrecv(const void *,int,MPI_Datatype,int,int,void *,int,MPI_Datatype,int,int,MPI_Comm,MPI_Status *);
int MPI_Abort(MPI_Comm,int);
int MPI_Query_thread(int *);
int MPI_Finalize(void);
#endif
