#include <stdio.h>
#include <mpi.h>
#include "qcd.h"
#include "solver.h"

int main(void)
{
    int threads = solver_getNThread();
    printf("Compiled solver threads per MPI rank: %d\n", threads);
    return threads == 1 ? 0 : 1;
}
