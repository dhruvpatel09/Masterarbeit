#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <mpi.h>

#include <qcd.h>
#include <DistEigvecsHdf5Reader.h>

#define NV 10
#define DEFAULT_RMAX 3
#define DEFAULT_TMAX 4

static int site3_index(const int x, const int y, const int z, const int L[4]) {
    return (x * L[2] + y) * L[3] + z;
}

static int site4_index(const int t, const int x, const int y, const int z, const int L[4]) {
    return ((t * L[1] + x) * L[2] + y) * L[3] + z;
}

static int read_env_int(
    const char *name,
    const int default_value,
    const int min_value,
    const int max_value
) {
    const char *text = getenv(name);

    if (text == NULL || *text == '\0') {
        return default_value;
    }

    errno = 0;
    char *end = NULL;
    const long value = strtol(text, &end, 10);

    if (errno == ERANGE ||
        end == text ||
        *end != '\0' ||
        value < min_value ||
        value > max_value) {
        fprintf(
            stderr,
            "ERROR: %s must be an integer in [%d,%d], got '%s'\n",
            name,
            min_value,
            max_value,
            text
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return (int)value;
}

static void read_evec_time(
    const char *fname,
    const int nc,
    const int time,
    const char *nbase,
    qcd_geometry *geo,
    qcd_spinorComponent3d *v
) {
    DistEigvecsHdf5Reader reader;
    memset(&reader, 0, sizeof(reader));
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

static qcd_complex_16 tau_time_path(
    const qcd_spinorComponent3d *v0,
    const qcd_spinorComponent3d *vT,
    const qcd_gaugeField *u,
    const int i,
    const int j,
    const int T,
    const int site3,
    const int x,
    const int y,
    const int z,
    const int L[4]
) {
    qcd_complex_16 vec[3];
    qcd_complex_16 tmp[3];

    for (int a = 0; a < 3; a++) {
        vec[a] = vT[j].D[site3][a];
    }

    /*
      Compute:
        U_0(x,t=0) U_0(x,t=1) ... U_0(x,t=T-1) v_j(x,T)

      We apply links backwards to the vector:
        vec <- U(T-1) vec
        ...
        vec <- U(0) vec
    */
    for (int tt = T - 1; tt >= 0; tt--) {
        const int site4 = site4_index(tt, x, y, z, L);

        for (int a = 0; a < 3; a++) {
            tmp[a].re = 0.0;
            tmp[a].im = 0.0;

            for (int b = 0; b < 3; b++) {
                tmp[a] = qcd_CADD(tmp[a], qcd_CMUL(u->D[site4][0][a][b], vec[b]));
            }
        }

        for (int a = 0; a < 3; a++) {
            vec[a] = tmp[a];
        }
    }

    qcd_complex_16 tau = {0.0, 0.0};
    for (int a = 0; a < 3; a++) {
        tau = qcd_CADD(tau, qcd_CMUL(qcd_CONJ(v0[i].D[site3][a]), vec[a]));
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
            fprintf(stderr, "ERROR: this test is intended for exactly 1 MPI rank.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const char *nbase = "Em1p4";

    int L[4] = {48, 24, 24, 24};
    int P[4] = {1, 1, 1, 1};

    const int nc = read_env_int("NCFG", 1, 1, INT_MAX);
    const int rmax = read_env_int(
        "RMAX",
        DEFAULT_RMAX,
        0,
        L[1] / 2
    );
    const int tmax = read_env_int(
        "TMAX",
        DEFAULT_TMAX,
        1,
        L[0] - 1
    );
    double theta[3] = {0.0, 0.0, 0.0};

    qcd_geometry geo;
    qcd_gaugeField u;
    qcd_spinorComponent3d v0[NV];
    qcd_spinorComponent3d vT[NV];

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
            qcd_initSpinorComponent3d(&vT[i], &geo)) {
            fprintf(stderr, "ERROR: qcd_initSpinorComponent3d failed\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    char gauge_file[1024];
    char evec_t0_file[1024];
    char evec_t_file[1024];

    snprintf(gauge_file, sizeof(gauge_file),
             "/home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n%d", nc);

    snprintf(evec_t0_file, sizeof(evec_t0_file),
             "/home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n%d/eigenvectors/Em1p4n%d_evec_t0.h5",
             nc, nc);

    if (myid == 0) {
        printf("Reading gauge field: %s\n", gauge_file);
        fflush(stdout);
    }

    if (qcd_getGaugeField(gauge_file, qcd_GF_OPENQCD, &u)) {
        fprintf(stderr, "ERROR: qcd_getGaugeField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (myid == 0) {
        printf("Reading source-time eigenvectors: %s\n", evec_t0_file);
        fflush(stdout);
    }

    read_evec_time(evec_t0_file, nc, 0, nbase, &geo, v0);

    if (myid == 0) {
        printf("\nRT spatial average scan: cfg=%d Nv=%d R/a=0..%d T/a=1..%d rho_i=1\n",
               nc, NV, rmax, tmax);
        printf(
            "META cfg=%d Nv=%d Rmin=0 Rmax=%d Tmin=1 Tmax=%d t0=0 rho=constant\n",
            nc,
            NV,
            rmax,
            tmax
        );
        printf("# DATA cfg T R Nsrc Re[avg L(R,T)] Im[avg L(R,T)]\n");
        fflush(stdout);
    }

    for (int T = 1; T <= tmax; T++) {
        snprintf(evec_t_file, sizeof(evec_t_file),
                 "/home/m2130292/Masterarbeit/mental/runs_Em1p4_Nv10_qcdnew_full/n%d/eigenvectors/Em1p4n%d_evec_t%d.h5",
                 nc, nc, T);

        if (myid == 0) {
            printf("Reading sink-time eigenvectors T=%d: %s\n", T, evec_t_file);
            fflush(stdout);
        }

        read_evec_time(evec_t_file, nc, T, nbase, &geo, vT);

        for (int R = 0; R <= rmax; R++) {
            qcd_complex_16 Lsum = {0.0, 0.0};
            long nsrc = 0;

            for (int xx = 0; xx < L[1]; xx++) {
                for (int yy = 0; yy < L[2]; yy++) {
                    for (int zz = 0; zz < L[3]; zz++) {
                        const int yx = (xx + R) % L[1];

                        const int x_site3 = site3_index(xx, yy, zz, L);
                        const int y_site3 = site3_index(yx, yy, zz, L);

                        qcd_complex_16 Lxy = {0.0, 0.0};

                        for (int i = 0; i < NV; i++) {
                            for (int j = 0; j < NV; j++) {
                                qcd_complex_16 tau_y =
                                    tau_time_path(v0, vT, &u, i, j, T, y_site3, yx, yy, zz, L);

                                qcd_complex_16 tau_x =
                                    tau_time_path(v0, vT, &u, i, j, T, x_site3, xx, yy, zz, L);

                                Lxy = qcd_CADD(Lxy, qcd_CMUL(tau_y, qcd_CONJ(tau_x)));
                            }
                        }

                        Lsum = qcd_CADD(Lsum, Lxy);
                        nsrc++;
                    }
                }
            }

            qcd_complex_16 Lavg = {Lsum.re / (double)nsrc, Lsum.im / (double)nsrc};

            if (myid == 0) {
                printf("DATA %d %d %d %ld %+.16e %+.16e\n",
                       nc, T, R, nsrc, Lavg.re, Lavg.im);
                fflush(stdout);
            }
        }
    }

    for (int i = 0; i < NV; i++) {
        qcd_destroySpinorComponent3d(&v0[i]);
        qcd_destroySpinorComponent3d(&vT[i]);
    }

    qcd_destroyGaugeField(&u);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();
    return 0;
}
