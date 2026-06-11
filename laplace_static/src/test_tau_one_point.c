#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include <qcd.h>
#include <DistEigvecsHdf5Reader.h>

#define NV 10

static void read_evec_time(
    const char *fname,
    const int nc,
    const int time,
    const char *nbase,
    qcd_geometry *geo,
    qcd_spinorComponent3d *v
) {
    DistEigvecsHdf5Reader reader;
    DistEigvecsHdf5_ComplexDouble *disteigvecs = NULL;

    int ret;
    char *tmp_run_name = NULL;
    int tmp_nc = -1;
    int tmp_time = -1;
    unsigned int tmp_nv = 0;
    int spatial_lat_size[3] = { geo->L[1], geo->L[2], geo->L[3] };
    int spatial_proclat_size[3] = {
        geo->L[1] / geo->lL[1],
        geo->L[2] / geo->lL[2],
        geo->L[3] / geo->lL[3]
    };
    int spatial_proclat_pos[3] = { geo->Pos[1], geo->Pos[2], geo->Pos[3] };

    DistEigvecsHdf5_Metadata metadata;
    DistEigvecsHdf5_Programdata programdata;

    ret = DistEigvecsHdf5Reader_Init(&reader);
    if (ret != 0) {
        fprintf(stderr, "ERROR: DistEigvecsHdf5Reader_Init failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    ret = DistEigvecsHdf5Reader_OpenFile(
        &reader,
        geo->comm3d,
        fname,
        &tmp_run_name,
        &tmp_nc,
        spatial_lat_size,
        spatial_proclat_size,
        &tmp_time,
        &tmp_nv,
        &metadata,
        &programdata
    );

    if (ret != 0) {
        fprintf(stderr, "ERROR: OpenFile failed for %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (tmp_nc != nc || tmp_time != time || tmp_nv < NV || strcmp(tmp_run_name, nbase) != 0) {
        fprintf(stderr,
                "ERROR: eigenvector metadata mismatch for %s: tmp_nc=%d tmp_time=%d tmp_nv=%u tmp_run_name=%s\n",
                fname, tmp_nc, tmp_time, tmp_nv, tmp_run_name ? tmp_run_name : "(null)");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    disteigvecs = malloc((size_t)geo->lV3 * 3 * NV * sizeof(DistEigvecsHdf5_ComplexDouble));
    if (!disteigvecs) {
        fprintf(stderr, "ERROR: malloc disteigvecs failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    ret = DistEigvecsHdf5Reader_ReadDistEigvecsRange(
        &reader,
        spatial_proclat_pos,
        0,
        NV,
        disteigvecs
    );

    if (ret != 0) {
        fprintf(stderr, "ERROR: ReadDistEigvecsRange failed for %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    size_t ncp = (size_t)geo->lV3 * 3 * sizeof(DistEigvecsHdf5_ComplexDouble);
    for (int j = 0; j < NV; j++) {
        memcpy(&(v[j].D[0][0].re), &(disteigvecs[(size_t)j * geo->lV3 * 3]), ncp);
    }

    free(disteigvecs);

    ret = DistEigvecsHdf5Reader_CloseFile(&reader);
    if (ret != 0) {
        fprintf(stderr, "ERROR: CloseFile failed for %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    DistEigvecsHdf5Reader_Free(&reader);
}

static qcd_complex_16 tau_one_link(
    const qcd_spinorComponent3d *v0,
    const qcd_spinorComponent3d *v1,
    const qcd_gaugeField *u,
    const int i,
    const int j,
    const int site4,
    const int site3
) {
    qcd_complex_16 tau = {0.0, 0.0};

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            qcd_complex_16 left = qcd_CONJ(v0[i].D[site3][a]);
            qcd_complex_16 mid  = u->D[site4][0][a][b];  // temporal link mu=0
            qcd_complex_16 right = v1[j].D[site3][b];

            qcd_complex_16 term = qcd_CMUL(left, qcd_CMUL(mid, right));
            tau = qcd_CADD(tau, term);
        }
    }

    return tau;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int myid = 0;
    int numprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (numprocs != 1) {
        if (myid == 0) {
            fprintf(stderr, "ERROR: this first test is intended for exactly 1 MPI rank.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const int nc = 1;
    const char *nbase = "Em1p4";

    int L[4] = {48, 24, 24, 24};
    int P[4] = {1, 1, 1, 1};
    double theta[3] = {0.0, 0.0, 0.0};

    qcd_geometry geo;
    qcd_gaugeField u;
    qcd_spinorComponent3d v0[NV];
    qcd_spinorComponent3d v1[NV];

    if (qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) {
        fprintf(stderr, "ERROR: qcd_initGeometry failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (qcd_initEO(&geo)) {
        fprintf(stderr, "ERROR: qcd_initEO failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (qcd_initGaugeField(&u, &geo)) {
        fprintf(stderr, "ERROR: qcd_initGaugeField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int i = 0; i < NV; i++) {
        if (qcd_initSpinorComponent3d(&v0[i], &geo) ||
            qcd_initSpinorComponent3d(&v1[i], &geo)) {
            fprintf(stderr, "ERROR: qcd_initSpinorComponent3d failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    char gauge_file[1024];
    char evec_t0_file[1024];
    char evec_t1_file[1024];

    snprintf(gauge_file, sizeof(gauge_file),
             "/home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n1");

    snprintf(evec_t0_file, sizeof(evec_t0_file),
             "/home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n1/eigenvectors/Em1p4n1_evec_t0.h5");

    snprintf(evec_t1_file, sizeof(evec_t1_file),
             "/home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n1/eigenvectors/Em1p4n1_evec_t1.h5");

    if (myid == 0) {
        printf("Reading gauge field: %s\n", gauge_file);
        fflush(stdout);
    }

    if (qcd_getGaugeField(gauge_file, qcd_GF_OPENQCD, &u)) {
        fprintf(stderr, "ERROR: qcd_getGaugeField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (myid == 0) {
        printf("Reading eigenvectors:\n");
        printf("  %s\n", evec_t0_file);
        printf("  %s\n", evec_t1_file);
        fflush(stdout);
    }

    read_evec_time(evec_t0_file, nc, 0, nbase, &geo, v0);
    read_evec_time(evec_t1_file, nc, 1, nbase, &geo, v1);

    // First test point: y=(0,0,0), t0=0.
    // With one MPI rank and lexicographic ordering, the all-zero site is local index 0.
    const int site3 = 0;
    const int site4 = 0;

    if (myid == 0) {
        printf("\nTau_ij(y=(0,0,0), t0=0, t1=1), Nv=%d\n", NV);
        printf("# i j Re(tau) Im(tau)\n");

        for (int i = 0; i < NV; i++) {
            for (int j = 0; j < NV; j++) {
                qcd_complex_16 tau = tau_one_link(v0, v1, &u, i, j, site4, site3);
                printf("%2d %2d %+.16e %+.16e\n", i, j, tau.re, tau.im);
            }
        }
    }

    for (int i = 0; i < NV; i++) {
        qcd_destroySpinorComponent3d(&v0[i]);
        qcd_destroySpinorComponent3d(&v1[i]);
    }

    qcd_destroyGaugeField(&u);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();
    return 0;
}
