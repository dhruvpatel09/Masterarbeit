changes with respect to 1.6 version:

* openMP is possible (set up via include/global.h)

* every solver has additional parameter istop

* deflation projection has additional parameters nmx_gcr, res_gcr

* when initializing MPI with openMP, please use solver_mpiInit instead of MPI_Init
