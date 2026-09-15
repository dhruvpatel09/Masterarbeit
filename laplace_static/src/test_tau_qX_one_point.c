#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <qcd.h>
#include <qcd_communication.h>
#include <DistEigvecsHdf5Reader.h>

#define NV 100
#define NT 192
#define NS 64
#define CFG_ID 100

#define DEFAULT_T0 96
#define DEFAULT_GX 0
#define DEFAULT_GY 0
#define DEFAULT_GZ 0

#define PATH_SIZE 2048

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

static const char *env_or_default(
    const char *name,
    const char *default_value
) {
    const char *value = getenv(name);

    if (value == NULL || *value == '\0') {
        return default_value;
    }

    return value;
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

    int spatial_lat_size[3] = {
        geo->L[1],
        geo->L[2],
        geo->L[3]
    };

    int spatial_proclat_size[3] = {
        geo->L[1] / geo->lL[1],
        geo->L[2] / geo->lL[2],
        geo->L[3] / geo->lL[3]
    };

    int spatial_proclat_pos[3] = {
        geo->Pos[1],
        geo->Pos[2],
        geo->Pos[3]
    };

    DistEigvecsHdf5_Metadata metadata;
    DistEigvecsHdf5_Programdata programdata;

    ret = DistEigvecsHdf5Reader_Init(&reader);

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: DistEigvecsHdf5Reader_Init failed\n"
        );
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
        fprintf(
            stderr,
            "ERROR: OpenFile failed for %s\n",
            fname
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (tmp_nc != nc ||
        tmp_time != time ||
        tmp_nv < NV ||
        tmp_run_name == NULL ||
        strcmp(tmp_run_name, nbase) != 0) {
        fprintf(
            stderr,
            "ERROR: eigenvector metadata mismatch for %s: "
            "cfg=%d time=%d Nv=%u run_name=%s\n",
            fname,
            tmp_nc,
            tmp_time,
            tmp_nv,
            tmp_run_name ? tmp_run_name : "(null)"
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    disteigvecs = malloc(
        (size_t)geo->lV3 *
        3U *
        NV *
        sizeof(*disteigvecs)
    );

    if (disteigvecs == NULL) {
        fprintf(
            stderr,
            "ERROR: allocation of distributed eigenvectors failed\n"
        );
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
        fprintf(
            stderr,
            "ERROR: ReadDistEigvecsRange failed for %s\n",
            fname
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const size_t mode_bytes =
        (size_t)geo->lV3 *
        3U *
        sizeof(*disteigvecs);

    for (int j = 0; j < NV; j++) {
        memcpy(
            &(v[j].D[0][0].re),
            &(disteigvecs[(size_t)j * geo->lV3 * 3U]),
            mode_bytes
        );
    }

    free(disteigvecs);

    ret = DistEigvecsHdf5Reader_CloseFile(&reader);

    if (ret != 0) {
        fprintf(
            stderr,
            "ERROR: CloseFile failed for %s\n",
            fname
        );
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
    const long site4,
    const long site3
) {
    qcd_complex_16 tau = {0.0, 0.0};

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            const qcd_complex_16 left =
                qcd_CONJ(v0[i].D[site3][a]);

            const qcd_complex_16 temporal_link =
                u->D[site4][0][a][b];

            const qcd_complex_16 right =
                v1[j].D[site3][b];

            tau = qcd_CADD(
                tau,
                qcd_CMUL(
                    left,
                    qcd_CMUL(
                        temporal_link,
                        right
                    )
                )
            );
        }
    }

    return tau;
}


static void mat3_zero(
    qcd_complex_16 a[3][3]
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j].re = 0.0;
            a[i][j].im = 0.0;
        }
    }
}

static void mat3_mul(
    qcd_complex_16 out[3][3],
    const qcd_complex_16 a[3][3],
    const qcd_complex_16 b[3][3]
) {
    qcd_complex_16 tmp[3][3];
    mat3_zero(tmp);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tmp[i][j] = qcd_CADD(
                    tmp[i][j],
                    qcd_CMUL(
                        a[i][k],
                        b[k][j]
                    )
                );
            }
        }
    }

    memcpy(out, tmp, sizeof(tmp));
}

static void mat3_adjoint(
    qcd_complex_16 out[3][3],
    const qcd_complex_16 in[3][3]
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = qcd_CONJ(
                in[j][i]
            );
        }
    }
}

static qcd_complex_16 mat3_trace(
    const qcd_complex_16 a[3][3]
) {
    qcd_complex_16 tr = {
        0.0,
        0.0
    };

    for (int i = 0; i < 3; i++) {
        tr = qcd_CADD(
            tr,
            a[i][i]
        );
    }

    return tr;
}

/*
  Exact rho_i = 1 factorization for R=0, T=1:

    sum_ij |tau_ij(x)|^2

      = Tr_3[
          U_0(x)
          K_t1(x,x)
          U_0(x)^\dagger
          K_t0(x,x)
        ]

  where

    K_t(x,x) = sum_n v_n(x,t) v_n(x,t)^\dagger.
*/
static qcd_complex_16 factorized_r0_one_link(
    const qcd_spinorComponent3d *v0,
    const qcd_spinorComponent3d *v1,
    const qcd_gaugeField *u,
    const long site4,
    const long site3
) {
    qcd_complex_16 k0[3][3];
    qcd_complex_16 k1[3][3];
    qcd_complex_16 temporal[3][3];
    qcd_complex_16 temporal_dag[3][3];

    qcd_complex_16 tmp1[3][3];
    qcd_complex_16 tmp2[3][3];
    qcd_complex_16 tmp3[3][3];

    mat3_zero(k0);
    mat3_zero(k1);

    for (int n = 0; n < NV; n++) {
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                k0[a][b] = qcd_CADD(
                    k0[a][b],
                    qcd_CMUL(
                        v0[n].D[site3][a],
                        qcd_CONJ(
                            v0[n].D[site3][b]
                        )
                    )
                );

                k1[a][b] = qcd_CADD(
                    k1[a][b],
                    qcd_CMUL(
                        v1[n].D[site3][a],
                        qcd_CONJ(
                            v1[n].D[site3][b]
                        )
                    )
                );
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            temporal[a][b] =
                u->D[site4][0][a][b];
        }
    }

    mat3_adjoint(
        temporal_dag,
        temporal
    );

    mat3_mul(
        tmp1,
        temporal,
        k1
    );

    mat3_mul(
        tmp2,
        tmp1,
        temporal_dag
    );

    mat3_mul(
        tmp3,
        tmp2,
        k0
    );

    return mat3_trace(tmp3);
}


static void build_laplace_kernel(
    qcd_complex_16 kernel[3][3],
    const qcd_spinorComponent3d *left,
    const long left_site,
    const qcd_spinorComponent3d *right,
    const long right_site
) {
    mat3_zero(kernel);

    for (int n = 0; n < NV; n++) {
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                kernel[a][b] = qcd_CADD(
                    kernel[a][b],
                    qcd_CMUL(
                        left[n].D[left_site][a],
                        qcd_CONJ(
                            right[n].D[right_site][b]
                        )
                    )
                );
            }
        }
    }
}


/*
  Exact rho_i = 1 factorization for a separated pair:

    L(x,y)
      = sum_ij tau_ij(y) tau_ij(x)^*

      = Tr_3[
          P_y
          K_t1(y,x)
          P_x^\dagger
          K_t0(x,y)
        ].
*/
static qcd_complex_16 factorized_pair_one_link(
    const qcd_spinorComponent3d *v0,
    const qcd_spinorComponent3d *v1,
    const qcd_gaugeField *u,
    const long x4,
    const long y4,
    const long x3,
    const long y3
) {
    qcd_complex_16 k1_yx[3][3];
    qcd_complex_16 k0_xy[3][3];

    qcd_complex_16 px[3][3];
    qcd_complex_16 py[3][3];
    qcd_complex_16 px_dag[3][3];

    qcd_complex_16 tmp1[3][3];
    qcd_complex_16 tmp2[3][3];
    qcd_complex_16 tmp3[3][3];

    build_laplace_kernel(
        k1_yx,
        v1,
        y3,
        v1,
        x3
    );

    build_laplace_kernel(
        k0_xy,
        v0,
        x3,
        v0,
        y3
    );

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            px[a][b] =
                u->D[x4][0][a][b];

            py[a][b] =
                u->D[y4][0][a][b];
        }
    }

    mat3_adjoint(
        px_dag,
        px
    );

    mat3_mul(
        tmp1,
        py,
        k1_yx
    );

    mat3_mul(
        tmp2,
        tmp1,
        px_dag
    );

    mat3_mul(
        tmp3,
        tmp2,
        k0_xy
    );

    return mat3_trace(tmp3);
}


static void check_r1_pair(
    const char *label,
    const qcd_spinorComponent3d *v0,
    const qcd_spinorComponent3d *v1,
    const qcd_gaugeField *u,
    const long x4,
    const long y4,
    const long x3,
    const long y3,
    const int gx,
    const int gy,
    const int gz,
    const qcd_geometry *geo
) {
    qcd_complex_16 direct = {
        0.0,
        0.0
    };

    for (int i = 0; i < NV; i++) {
        for (int j = 0; j < NV; j++) {
            const qcd_complex_16 tau_x =
                tau_one_link(
                    v0,
                    v1,
                    u,
                    i,
                    j,
                    x4,
                    x3
                );

            const qcd_complex_16 tau_y =
                tau_one_link(
                    v0,
                    v1,
                    u,
                    i,
                    j,
                    y4,
                    y3
                );

            direct = qcd_CADD(
                direct,
                qcd_CMUL(
                    tau_y,
                    qcd_CONJ(tau_x)
                )
            );
        }
    }

    const qcd_complex_16 factorized =
        factorized_pair_one_link(
            v0,
            v1,
            u,
            x4,
            y4,
            x3,
            y3
        );

    const double delta_re =
        factorized.re -
        direct.re;

    const double delta_im =
        factorized.im -
        direct.im;

    const double absdiff =
        sqrt(
            delta_re * delta_re +
            delta_im * delta_im
        );

    const double direct_abs =
        sqrt(
            direct.re * direct.re +
            direct.im * direct.im
        );

    const double scale =
        direct_abs > 0.0
        ? direct_abs
        : 1.0;

    const double reldiff =
        absdiff / scale;

    printf(
        "PAIR_R1 "
        "kind=%s "
        "x=(%d,%d,%d) "
        "y=(%d,%d,%d) "
        "x3=%ld y3=%ld "
        "x4=%ld y4=%ld "
        "y3_halo=%d "
        "y4_halo=%d "
        "direct=%+.17e %+.17e "
        "factorized=%+.17e %+.17e "
        "absdiff=%.17e "
        "reldiff=%.17e "
        "status=%s\n",
        label,
        gx,
        gy,
        gz,
        gx + 1,
        gy,
        gz,
        x3,
        y3,
        x4,
        y4,
        y3 >= geo->lV3,
        y4 >= geo->lV,
        direct.re,
        direct.im,
        factorized.re,
        factorized.im,
        absdiff,
        reldiff,
        reldiff < 1.0e-10
            ? "PASS"
            : "FAIL"
    );
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int myid = 0;
    int numprocs = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    const int L[4] = {
        NT,
        NS,
        NS,
        NS
    };

    /*
      Keep the complete temporal direction local.

      Spatial decomposition:
          64^3 -> 4 x 4 x 4 MPI blocks
                -> local 16^3 spatial volume.
    */
    const int P[4] = {
        1,
        4,
        4,
        4
    };

    const int expected_procs =
        P[0] * P[1] * P[2] * P[3];

    if (numprocs != expected_procs) {
        if (myid == 0) {
            fprintf(
                stderr,
                "ERROR: expected %d MPI ranks for process grid "
                "%d x %d x %d x %d, got %d\n",
                expected_procs,
                P[0],
                P[1],
                P[2],
                P[3],
                numprocs
            );
        }

        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    const int t0 = read_env_int(
        "T0",
        DEFAULT_T0,
        0,
        NT - 2
    );

    /*
      First qX smoke test deliberately uses exactly one forward
      temporal link.  No modulo Nt and no temporal wrapping:
      this is compatible with the open temporal boundary.
    */
    const int t1 = t0 + 1;

    const int gx = read_env_int(
        "GX",
        DEFAULT_GX,
        0,
        NS - 1
    );

    const int gy = read_env_int(
        "GY",
        DEFAULT_GY,
        0,
        NS - 1
    );

    const int gz = read_env_int(
        "GZ",
        DEFAULT_GZ,
        0,
        NS - 1
    );

    const int print_all = read_env_int(
        "PRINT_ALL",
        0,
        0,
        1
    );

    const char *qX_root = env_or_default(
        "QX_ROOT",
        "/hppfs/scratch/00/di98qes3/qX"
    );

    const char *evec_root = env_or_default(
        "EVROOT",
        "/hppfs/scratch/0E/di24gus/Masterarbeit/qX/"
        "mental_runs/"
        "n100_evec_Nv100_deg10_nf7_mr30_up015_p12x4x4x4/"
        "eigenvectors"
    );

    const char *nbase = "qX_p2";

    double theta[4] = {
        0.0,
        0.0,
        0.0,
        0.0
    };

    qcd_geometry geo;
    qcd_gaugeField u;

    if (qcd_initGeometry(
            &geo,
            L,
            P,
            theta,
            myid,
            numprocs
        )) {
        fprintf(
            stderr,
            "process %d: ERROR: qcd_initGeometry failed\n",
            myid
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    if (qcd_initEO(&geo)) {
        fprintf(
            stderr,
            "process %d: ERROR: qcd_initEO failed\n",
            myid
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    if (qcd_initGaugeField(&u, &geo)) {
        fprintf(
            stderr,
            "process %d: ERROR: qcd_initGaugeField failed\n",
            myid
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    qcd_spinorComponent3d v0[NV];
    qcd_spinorComponent3d v1[NV];

    for (int i = 0; i < NV; i++) {
        if (qcd_initSpinorComponent3d(
                &v0[i],
                &geo
            ) ||
            qcd_initSpinorComponent3d(
                &v1[i],
                &geo
            )) {
            fprintf(
                stderr,
                "process %d: ERROR: "
                "qcd_initSpinorComponent3d failed at mode %d\n",
                myid,
                i
            );
            MPI_Abort(
                MPI_COMM_WORLD,
                EXIT_FAILURE
            );
        }
    }

    char gauge_file[PATH_SIZE];
    char evec_t0_file[PATH_SIZE];
    char evec_t1_file[PATH_SIZE];

    snprintf(
        gauge_file,
        sizeof(gauge_file),
        "%s/qX_p2n%d",
        qX_root,
        CFG_ID
    );

    snprintf(
        evec_t0_file,
        sizeof(evec_t0_file),
        "%s/qX_p2n%d_evec_t%d.h5",
        evec_root,
        CFG_ID,
        t0
    );

    snprintf(
        evec_t1_file,
        sizeof(evec_t1_file),
        "%s/qX_p2n%d_evec_t%d.h5",
        evec_root,
        CFG_ID,
        t1
    );

    if (myid == 0) {
        printf(
            "qX one-point tau smoke test\n"
        );
        printf(
            "L=%d x %d x %d x %d\n",
            L[0],
            L[1],
            L[2],
            L[3]
        );
        printf(
            "P=%d x %d x %d x %d\n",
            P[0],
            P[1],
            P[2],
            P[3]
        );
        printf(
            "Nv=%d cfg=%d temporal_bc=open\n",
            NV,
            CFG_ID
        );
        printf(
            "t0=%d t1=%d global_xyz=(%d,%d,%d)\n",
            t0,
            t1,
            gx,
            gy,
            gz
        );
        printf(
            "Gauge field: %s\n",
            gauge_file
        );
        printf(
            "Eigenvectors t0: %s\n",
            evec_t0_file
        );
        printf(
            "Eigenvectors t1: %s\n",
            evec_t1_file
        );
        fflush(stdout);
    }

    /*
      qcd_getGaugeField performs collective parallel I/O.
      All MPI ranks must enter this call.
    */
    if (qcd_getGaugeField(
            gauge_file,
            qcd_GF_OPENQCD,
            &u
        )) {
        fprintf(
            stderr,
            "process %d: ERROR: "
            "qcd_getGaugeField failed for %s\n",
            myid,
            gauge_file
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    read_evec_time(
        evec_t0_file,
        CFG_ID,
        t0,
        nbase,
        &geo,
        v0
    );

    read_evec_time(
        evec_t1_file,
        CFG_ID,
        t1,
        nbase,
        &geo,
        v1
    );

    /*
      Populate the + spatial halo.

      qcd_communicateGaugeM:
          send local - face
          receive + neighbour into Bplus.

      Likewise for qcd_spinorComponent3d.
    */
    if (qcd_communicateGaugeM(&u)) {
        fprintf(
            stderr,
            "process %d: ERROR: qcd_communicateGaugeM failed\n",
            myid
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    qcd_waitall(&geo);

    for (int n = 0; n < NV; n++) {
        if (qcd_communicateSpinorComponent3dM(
                &v0[n]
            )) {
            fprintf(
                stderr,
                "process %d: ERROR: "
                "v0 halo exchange failed for mode %d\n",
                myid,
                n
            );
            MPI_Abort(
                MPI_COMM_WORLD,
                EXIT_FAILURE
            );
        }

        qcd_waitall(&geo);

        if (qcd_communicateSpinorComponent3dM(
                &v1[n]
            )) {
            fprintf(
                stderr,
                "process %d: ERROR: "
                "v1 halo exchange failed for mode %d\n",
                myid,
                n
            );
            MPI_Abort(
                MPI_COMM_WORLD,
                EXIT_FAILURE
            );
        }

        qcd_waitall(&geo);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
      Determine which spatial MPI block owns the requested
      global spatial point.
    */
    const int owner_px = gx / geo.lL[1];
    const int owner_py = gy / geo.lL[2];
    const int owner_pz = gz / geo.lL[3];

    const int owns_point =
        geo.Pos[1] == owner_px &&
        geo.Pos[2] == owner_py &&
        geo.Pos[3] == owner_pz;

    if (owns_point) {
        const int lx =
            gx - geo.Pos[1] * geo.lL[1];

        const int ly =
            gy - geo.Pos[2] * geo.lL[2];

        const int lz =
            gz - geo.Pos[3] * geo.lL[3];

        /*
          P[0]=1, therefore the local temporal coordinate is
          identical to the global temporal coordinate.
        */
        const int lt0 = t0;

        const long site3 =
            qcd_LEXIC0(
                lx,
                ly,
                lz,
                geo.lL
            );

        const long site4 =
            qcd_LEXIC(
                lt0,
                lx,
                ly,
                lz,
                geo.lL
            );

        qcd_complex_16 trace = {
            0.0,
            0.0
        };

        double frobenius_sq = 0.0;
        double maximum_abs = -1.0;
        int maximum_i = -1;
        int maximum_j = -1;

        qcd_complex_16 selected[4][4];

        memset(
            selected,
            0,
            sizeof(selected)
        );

        if (myid >= 0) {
            printf(
                "OWNER rank=%d Pos=(%d,%d,%d,%d) "
                "local_xyz=(%d,%d,%d) "
                "site3=%ld site4=%ld\n",
                myid,
                geo.Pos[0],
                geo.Pos[1],
                geo.Pos[2],
                geo.Pos[3],
                lx,
                ly,
                lz,
                site3,
                site4
            );

            if (print_all) {
                printf(
                    "# TAU i j Re Im\n"
                );
            }
        }

        for (int i = 0; i < NV; i++) {
            for (int j = 0; j < NV; j++) {
                const qcd_complex_16 tau =
                    tau_one_link(
                        v0,
                        v1,
                        &u,
                        i,
                        j,
                        site4,
                        site3
                    );

                if (i == j) {
                    trace = qcd_CADD(
                        trace,
                        tau
                    );
                }

                const double abs_sq =
                    tau.re * tau.re +
                    tau.im * tau.im;

                frobenius_sq += abs_sq;

                const double absolute =
                    sqrt(abs_sq);

                if (absolute > maximum_abs) {
                    maximum_abs = absolute;
                    maximum_i = i;
                    maximum_j = j;
                }

                if (i < 4 && j < 4) {
                    selected[i][j] = tau;
                }

                if (print_all) {
                    printf(
                        "TAU %d %d %+.17e %+.17e\n",
                        i,
                        j,
                        tau.re,
                        tau.im
                    );
                }
            }
        }

        const qcd_complex_16 factorized_r0 =
            factorized_r0_one_link(
                v0,
                v1,
                &u,
                site4,
                site3
            );

        const double direct_r0 =
            frobenius_sq;

        const double delta_re =
            factorized_r0.re -
            direct_r0;

        const double absdiff =
            sqrt(
                delta_re * delta_re +
                factorized_r0.im *
                factorized_r0.im
            );

        const double scale =
            direct_r0 > 0.0
            ? direct_r0
            : 1.0;

        const double reldiff =
            absdiff / scale;

        printf(
            "FACTOR_R0 "
            "direct=%+.17e "
            "factorized=%+.17e %+.17e "
            "absdiff=%.17e "
            "reldiff=%.17e "
            "status=%s\n",
            direct_r0,
            factorized_r0.re,
            factorized_r0.im,
            absdiff,
            reldiff,
            reldiff < 1.0e-10
                ? "PASS"
                : "FAIL"
        );

        printf(
            "SUMMARY trace=%+.17e %+.17e "
            "frobenius=%.17e "
            "maxabs=%.17e max_i=%d max_j=%d\n",
            trace.re,
            trace.im,
            sqrt(frobenius_sq),
            maximum_abs,
            maximum_i,
            maximum_j
        );

        printf(
            "# top-left 4 x 4 block\n"
        );

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                printf(
                    "SELECTED %d %d %+.17e %+.17e\n",
                    i,
                    j,
                    selected[i][j].re,
                    selected[i][j].im
                );
            }
        }

        fflush(stdout);
    }

    /*
      Controlled R=1 tests on the spatial block at
      process position (0,0,0).

        14 -> 15 : both endpoints local
        15 -> 16 : y endpoint lies in +x halo

      For both cases use qcd's own plus/plus3d maps.
    */
    if (
        geo.Pos[0] == 0 &&
        geo.Pos[1] == 0 &&
        geo.Pos[2] == 0 &&
        geo.Pos[3] == 0
    ) {
        const int ly = 0;
        const int lz = 0;
        const int lt = t0;

        const int pair_lx[2] = {
            14,
            15
        };

        const char *pair_label[2] = {
            "local",
            "cross_boundary"
        };

        for (int pcase = 0; pcase < 2; pcase++) {
            const int lx =
                pair_lx[pcase];

            const long x3 =
                qcd_LEXIC0(
                    lx,
                    ly,
                    lz,
                    geo.lL
                );

            const long y3 =
                geo.plus3d[x3][1];

            const long x4 =
                qcd_LEXIC(
                    lt,
                    lx,
                    ly,
                    lz,
                    geo.lL
                );

            const long y4 =
                geo.plus[x4][1];

            check_r1_pair(
                pair_label[pcase],
                v0,
                v1,
                &u,
                x4,
                y4,
                x3,
                y3,
                lx,
                0,
                0,
                &geo
            );
        }

        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < NV; i++) {
        qcd_destroySpinorComponent3d(
            &v0[i]
        );

        qcd_destroySpinorComponent3d(
            &v1[i]
        );
    }

    qcd_destroyGaugeField(&u);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();

    return 0;
}
