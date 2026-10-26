/* Single-process test shim only. Production builds use the system MPI header. */
#ifndef QX_TEST_MPI_H
#define QX_TEST_MPI_H
#include <string.h>
#include <stdlib.h>
#include <time.h>
typedef int MPI_Comm;
typedef int MPI_Datatype;
typedef int MPI_Request;
#define MPI_COMM_WORLD 0
#define MPI_SUCCESS 0
#define MPI_THREAD_FUNNELED 1
#define MPI_INT 1
#define MPI_DOUBLE 2
#define MPI_UNSIGNED_LONG_LONG 3
#define MPI_SUM 1
#define MPI_MIN 2
#define MPI_MAX 3
static inline int MPI_Bcast(void *p,int n,int type,int root,int comm)
{(void)p;(void)n;(void)type;(void)root;(void)comm;return 0;}
static inline double MPI_Wtime(void) {return (double)clock()/CLOCKS_PER_SEC;}
static inline int MPI_Init_thread(int *argc,char ***argv,int required,int *provided)
{(void)argc;(void)argv;*provided=required;return 0;}
static inline int MPI_Comm_rank(int c,int *r){(void)c;*r=0;return 0;}
static inline int MPI_Comm_size(int c,int *r){(void)c;*r=1;return 0;}
static inline int MPI_Allreduce(const void *a,void *b,int n,int type,int op,int comm)
{(void)op;(void)comm;memcpy(b,a,n*(type==MPI_DOUBLE?sizeof(double):type==MPI_INT?sizeof(int):sizeof(unsigned long long)));return 0;}
static inline int MPI_Abort(int c,int code){(void)c;exit(code);return 0;}
static inline int MPI_Finalize(void){return 0;}
#endif
