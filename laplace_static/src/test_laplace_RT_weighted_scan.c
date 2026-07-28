#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>

#include <qcd.h>
#include <DistEigvecsHdf5Reader.h>

#define NV 10
#define MAX_PROFILES NV
#define PAPER_NPROFILES 7
#define NT 48
#define DEFAULT_RMAX 3
#define DEFAULT_TMAX 4

static const double paper_profile_sigmas[PAPER_NPROFILES] = {
    0.0500,
    0.0894,
    0.1289,
    0.1683,
    0.2078,
    0.2472,
    0.2867
};

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


typedef enum {
    T0_MODE_FIXED = 0,
    T0_MODE_AVERAGE = 1
} t0_mode_t;

static t0_mode_t read_t0_mode(void) {
    const char *text = getenv("T0_MODE");

    if (text == NULL || *text == '\0' ||
        strcmp(text, "fixed") == 0) {
        return T0_MODE_FIXED;
    }

    if (strcmp(text, "average") == 0) {
        return T0_MODE_AVERAGE;
    }

    fprintf(
        stderr,
        "ERROR: T0_MODE must be 'fixed' or 'average', got '%s'\n",
        text
    );
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return T0_MODE_FIXED;
}

static const char *t0_mode_name(const t0_mode_t mode) {
    return (
        mode == T0_MODE_AVERAGE
        ? "average"
        : "fixed"
    );
}

typedef enum {
    PROFILE_BASIS_PAPER = 0,
    PROFILE_BASIS_DELTA = 1
} profile_basis_t;

static profile_basis_t read_profile_basis(void) {
    const char *text = getenv("PROFILE_BASIS");

    if (text == NULL || *text == '\0' ||
        strcmp(text, "paper") == 0) {
        return PROFILE_BASIS_PAPER;
    }

    if (strcmp(text, "delta") == 0) {
        return PROFILE_BASIS_DELTA;
    }

    fprintf(
        stderr,
        "ERROR: PROFILE_BASIS must be 'paper' or 'delta', "
        "got '%s'\n",
        text
    );
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

    return PROFILE_BASIS_PAPER;
}

static const char *profile_basis_name(
    const profile_basis_t basis
) {
    return (
        basis == PROFILE_BASIS_DELTA
        ? "delta"
        : "paper"
    );
}

static const char *profile_rho_name(
    const profile_basis_t basis
) {
    return (
        basis == PROFILE_BASIS_DELTA
        ? "constant+delta"
        : "constant+Gaussian"
    );
}

static const char *profile_weight_name(
    const profile_basis_t basis
) {
    return (
        basis == PROFILE_BASIS_DELTA
        ? "delta_ki"
        : "exp(-lambda^2/(2*sigma^2))"
    );
}

static void abort_hdf5(
    const char *operation,
    const char *fname
) {
    fprintf(
        stderr,
        "ERROR: %s failed for %s\n",
        operation,
        fname
    );
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static void read_scalar_string(
    const hid_t file,
    const char *dataset_name,
    const char *fname,
    char *value,
    const size_t value_size
) {
    hid_t dataset = H5Dopen2(
        file,
        dataset_name,
        H5P_DEFAULT
    );

    if (dataset < 0) {
        abort_hdf5("H5Dopen2 string dataset", fname);
    }

    hid_t datatype = H5Dget_type(dataset);

    if (datatype < 0 ||
        H5Tget_class(datatype) != H5T_STRING) {
        abort_hdf5("H5Dget_type string dataset", fname);
    }

    if (H5Tis_variable_str(datatype) > 0) {
        char *temporary = NULL;

        if (
            H5Dread(
                dataset,
                datatype,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &temporary
            ) < 0
        ) {
            abort_hdf5("H5Dread variable string", fname);
        }

        snprintf(
            value,
            value_size,
            "%s",
            temporary != NULL ? temporary : ""
        );

        if (temporary != NULL) {
            H5free_memory(temporary);
        }
    } else {
        const size_t stored_size = H5Tget_size(datatype);
        char *temporary = calloc(
            stored_size + 1,
            sizeof(*temporary)
        );

        if (temporary == NULL) {
            fprintf(
                stderr,
                "ERROR: allocation for %s in %s failed\n",
                dataset_name,
                fname
            );
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (
            H5Dread(
                dataset,
                datatype,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                temporary
            ) < 0
        ) {
            abort_hdf5("H5Dread fixed string", fname);
        }

        snprintf(
            value,
            value_size,
            "%s",
            temporary
        );
        free(temporary);
    }

    if (H5Tclose(datatype) < 0) {
        abort_hdf5("H5Tclose string datatype", fname);
    }

    if (H5Dclose(dataset) < 0) {
        abort_hdf5("H5Dclose string dataset", fname);
    }
}

static void read_eigenvalues(
    const char *fname,
    const int expected_cfg,
    const char *expected_run_name,
    double lambda[NT][NV]
) {
    hid_t file = H5Fopen(
        fname,
        H5F_ACC_RDONLY,
        H5P_DEFAULT
    );

    if (file < 0) {
        abort_hdf5("H5Fopen", fname);
    }

    hid_t dataset = H5Dopen2(
        file,
        "/disteigvals",
        H5P_DEFAULT
    );

    if (dataset < 0) {
        abort_hdf5("H5Dopen2 /disteigvals", fname);
    }

    hid_t filespace = H5Dget_space(dataset);

    if (filespace < 0) {
        abort_hdf5("H5Dget_space /disteigvals", fname);
    }

    hsize_t dimensions[2] = {0, 0};

    if (
        H5Sget_simple_extent_ndims(filespace) != 2 ||
        H5Sget_simple_extent_dims(
            filespace,
            dimensions,
            NULL
        ) < 0 ||
        dimensions[0] != NT ||
        dimensions[1] != NV
    ) {
        fprintf(
            stderr,
            "ERROR: /disteigvals in %s has shape "
            "(%llu,%llu), expected (%d,%d)\n",
            fname,
            (unsigned long long)dimensions[0],
            (unsigned long long)dimensions[1],
            NT,
            NV
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (
        H5Dread(
            dataset,
            H5T_NATIVE_DOUBLE,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            lambda
        ) < 0
    ) {
        abort_hdf5("H5Dread /disteigvals", fname);
    }

    if (H5Sclose(filespace) < 0) {
        abort_hdf5("H5Sclose /disteigvals", fname);
    }

    if (H5Dclose(dataset) < 0) {
        abort_hdf5("H5Dclose /disteigvals", fname);
    }

    int times[NT];

    dataset = H5Dopen2(
        file,
        "/times",
        H5P_DEFAULT
    );

    if (dataset < 0) {
        abort_hdf5("H5Dopen2 /times", fname);
    }

    filespace = H5Dget_space(dataset);

    hsize_t time_dimensions[1] = {0};

    if (
        filespace < 0 ||
        H5Sget_simple_extent_ndims(filespace) != 1 ||
        H5Sget_simple_extent_dims(
            filespace,
            time_dimensions,
            NULL
        ) < 0 ||
        time_dimensions[0] != NT
    ) {
        fprintf(
            stderr,
            "ERROR: /times in %s has length %llu, "
            "expected %d\n",
            fname,
            (unsigned long long)time_dimensions[0],
            NT
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (
        H5Dread(
            dataset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            times
        ) < 0) {
        abort_hdf5("H5Dread /times", fname);
    }

    if (H5Sclose(filespace) < 0) {
        abort_hdf5("H5Sclose /times", fname);
    }

    if (H5Dclose(dataset) < 0) {
        abort_hdf5("H5Dclose /times", fname);
    }

    for (int t = 0; t < NT; t++) {
        if (times[t] != t) {
            fprintf(
                stderr,
                "ERROR: /times[%d]=%d in %s\n",
                t,
                times[t],
                fname
            );
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    int cfg_id = -1;

    dataset = H5Dopen2(
        file,
        "/cnfg_id",
        H5P_DEFAULT
    );

    if (dataset < 0 ||
        H5Dread(
            dataset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            &cfg_id
        ) < 0) {
        abort_hdf5("read /cnfg_id", fname);
    }

    if (H5Dclose(dataset) < 0) {
        abort_hdf5("H5Dclose /cnfg_id", fname);
    }

    if (cfg_id != expected_cfg) {
        fprintf(
            stderr,
            "ERROR: /cnfg_id=%d in %s, expected %d\n",
            cfg_id,
            fname,
            expected_cfg
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char run_name[128];

    read_scalar_string(
        file,
        "/run_name",
        fname,
        run_name,
        sizeof(run_name)
    );

    if (strcmp(run_name, expected_run_name) != 0) {
        fprintf(
            stderr,
            "ERROR: /run_name='%s' in %s, expected '%s'\n",
            run_name,
            fname,
            expected_run_name
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (H5Fclose(file) < 0) {
        abort_hdf5("H5Fclose", fname);
    }

    for (int t = 0; t < NT; t++) {
        for (int i = 0; i < NV; i++) {
            if (!isfinite(lambda[t][i])) {
                fprintf(
                    stderr,
                    "ERROR: non-finite eigenvalue at "
                    "t=%d i=%d in %s\n",
                    t,
                    i,
                    fname
                );
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            if (
                i > 0 &&
                lambda[t][i] <
                    lambda[t][i - 1] - 1.0e-13
            ) {
                fprintf(
                    stderr,
                    "ERROR: unordered eigenvalues at "
                    "t=%d i=%d in %s\n",
                    t,
                    i,
                    fname
                );
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
    }
}

static void add_scaled_complex(
    qcd_complex_16 *sum,
    const qcd_complex_16 value,
    const double weight
) {
    sum->re += weight * value.re;
    sum->im += weight * value.im;
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
    const qcd_spinorComponent3d *v_src,
    const qcd_spinorComponent3d *v_sink,
    const qcd_gaugeField *u,
    const int i,
    const int j,
    const int t0,
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
        vec[a] = v_sink[j].D[site3][a];
    }

    /*
      Compute the forward temporal transporter

        U_0(x,t0) U_0(x,t0+1) ... U_0(x,t0+T-1)
        v_j(x,t0+T),

      with periodic temporal wrapping.

      The matrices are applied backwards to the sink vector.
    */
    for (int step = T - 1; step >= 0; step--) {
        const int tt = (t0 + step) % L[0];
        const int site4 = site4_index(
            tt,
            x,
            y,
            z,
            L
        );

        for (int a = 0; a < 3; a++) {
            tmp[a].re = 0.0;
            tmp[a].im = 0.0;

            for (int b = 0; b < 3; b++) {
                tmp[a] = qcd_CADD(
                    tmp[a],
                    qcd_CMUL(
                        u->D[site4][0][a][b],
                        vec[b]
                    )
                );
            }
        }

        for (int a = 0; a < 3; a++) {
            vec[a] = tmp[a];
        }
    }

    qcd_complex_16 tau = {0.0, 0.0};

    for (int a = 0; a < 3; a++) {
        tau = qcd_CADD(
            tau,
            qcd_CMUL(
                qcd_CONJ(v_src[i].D[site3][a]),
                vec[a]
            )
        );
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
            fprintf(
                stderr,
                "ERROR: this test is intended for exactly "
                "1 MPI rank.\n"
            );
        }

        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    const char *nbase = "Em1p4";

    int L[4] = {NT, 24, 24, 24};
    int P[4] = {1, 1, 1, 1};

    const int nc = read_env_int(
        "NCFG",
        1,
        1,
        INT_MAX
    );

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

    const t0_mode_t t0_mode = read_t0_mode();
    const profile_basis_t profile_basis =
        read_profile_basis();
    const int nprofiles = (
        profile_basis == PROFILE_BASIS_DELTA
        ? NV
        : PAPER_NPROFILES
    );

    const int t0_fixed = read_env_int(
        "T0_FIXED",
        0,
        0,
        L[0] - 1
    );

    const int t0_start = read_env_int(
        "T0_START",
        0,
        0,
        L[0] - 1
    );

    const int t0_stride = read_env_int(
        "T0_STRIDE",
        1,
        1,
        L[0]
    );

    int t0_values[L[0]];
    int nt0 = 0;

    if (t0_mode == T0_MODE_FIXED) {
        t0_values[nt0++] = t0_fixed;
    } else {
        for (
            int t0 = t0_start;
            t0 < L[0];
            t0 += t0_stride
        ) {
            t0_values[nt0++] = t0;
        }
    }

    if (nt0 < 1) {
        fprintf(
            stderr,
            "ERROR: no temporal source times selected\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    unsigned char needed_time[L[0]];
    unsigned char initialized_time[L[0]];

    memset(
        needed_time,
        0,
        sizeof(needed_time)
    );

    memset(
        initialized_time,
        0,
        sizeof(initialized_time)
    );

    for (int source = 0; source < nt0; source++) {
        const int t0 = t0_values[source];

        needed_time[t0] = 1;

        for (int T = 1; T <= tmax; T++) {
            const int t1 = (t0 + T) % L[0];
            needed_time[t1] = 1;
        }
    }

    int n_needed_times = 0;

    for (int t = 0; t < L[0]; t++) {
        if (needed_time[t]) {
            n_needed_times++;
        }
    }

    double theta[3] = {0.0, 0.0, 0.0};

    qcd_geometry geo;
    qcd_gaugeField u;

    if (
        qcd_initGeometry(
            &geo,
            L,
            P,
            theta,
            myid,
            numprocs
        )
    ) {
        fprintf(
            stderr,
            "ERROR: qcd_initGeometry failed\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    if (qcd_initEO(&geo)) {
        fprintf(
            stderr,
            "ERROR: qcd_initEO failed\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    if (qcd_initGaugeField(&u, &geo)) {
        fprintf(
            stderr,
            "ERROR: qcd_initGaugeField failed\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    qcd_spinorComponent3d *v_by_time = calloc(
        (size_t)L[0] * NV,
        sizeof(*v_by_time)
    );

    if (v_by_time == NULL) {
        fprintf(
            stderr,
            "ERROR: allocation of temporal eigenvectors failed\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    char gauge_file[1024];
    char evec_file[1024];
    char eval_file[1024];

    double lambda[NT][NV];
    double weights[MAX_PROFILES][NT][NV];

    snprintf(
        eval_file,
        sizeof(eval_file),
        "/home/m2130292/Masterarbeit/mental/"
        "runs_Em1p4_Nv10_qcdnew_full/n%d/"
        "eigenvalues/"
        "Em1p4_mental_n%d_Nv10_qcdnew_n%d_eval.h5",
        nc,
        nc,
        nc
    );

    if (myid == 0) {
        printf(
            "Reading eigenvalues: %s\n",
            eval_file
        );
        fflush(stdout);
    }

    read_eigenvalues(
        eval_file,
        nc,
        nbase,
        lambda
    );

    for (int k = 0; k < nprofiles; k++) {
        for (int t = 0; t < NT; t++) {
            for (int i = 0; i < NV; i++) {
                if (profile_basis == PROFILE_BASIS_DELTA) {
                    weights[k][t][i] = (
                        k == i
                        ? 1.0
                        : 0.0
                    );
                } else {
                    const double sigma =
                        paper_profile_sigmas[k];
                    const double scaled =
                        lambda[t][i] / sigma;

                    /*
                      The Gaussian quark profile is

                        rho_i =
                          exp[-lambda_i^2/(4 sigma^2)].

                      Equation (8) contains rho_i^2 at
                      each temporal endpoint, hence the
                      direct correlator weight below.
                    */
                    weights[k][t][i] = exp(
                        -0.5 * scaled * scaled
                    );
                }
            }
        }
    }

    snprintf(
        gauge_file,
        sizeof(gauge_file),
        "/home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n%d",
        nc
    );

    if (myid == 0) {
        printf(
            "Reading gauge field: %s\n",
            gauge_file
        );
        fflush(stdout);
    }

    if (
        qcd_getGaugeField(
            gauge_file,
            qcd_GF_OPENQCD,
            &u
        )
    ) {
        fprintf(
            stderr,
            "ERROR: qcd_getGaugeField failed\n"
        );
        MPI_Abort(
            MPI_COMM_WORLD,
            EXIT_FAILURE
        );
    }

    for (int t = 0; t < L[0]; t++) {
        if (!needed_time[t]) {
            continue;
        }

        qcd_spinorComponent3d *v_time = (
            &v_by_time[(size_t)t * NV]
        );

        for (int i = 0; i < NV; i++) {
            if (
                qcd_initSpinorComponent3d(
                    &v_time[i],
                    &geo
                )
            ) {
                fprintf(
                    stderr,
                    "ERROR: qcd_initSpinorComponent3d "
                    "failed at t=%d i=%d\n",
                    t,
                    i
                );
                MPI_Abort(
                    MPI_COMM_WORLD,
                    EXIT_FAILURE
                );
            }
        }

        initialized_time[t] = 1;

        snprintf(
            evec_file,
            sizeof(evec_file),
            "/home/m2130292/Masterarbeit/mental/"
            "runs_Em1p4_Nv10_qcdnew_full/n%d/"
            "eigenvectors/Em1p4n%d_evec_t%d.h5",
            nc,
            nc,
            t
        );

        if (myid == 0) {
            printf(
                "Reading temporal eigenvectors t=%d: %s\n",
                t,
                evec_file
            );
            fflush(stdout);
        }

        read_evec_time(
            evec_file,
            nc,
            t,
            nbase,
            &geo,
            v_time
        );
    }

    if (myid == 0) {
        printf(
            "\nRT source-average scan: "
            "cfg=%d Nv=%d R/a=0..%d T/a=1..%d "
            "T0Mode=%s Nt0=%d "
            "rho=%s NProfiles=%d\n",
            nc,
            NV,
            rmax,
            tmax,
            t0_mode_name(t0_mode),
            nt0,
            profile_rho_name(profile_basis),
            nprofiles
        );

        printf(
            "META cfg=%d Nv=%d Rmin=0 Rmax=%d "
            "Tmin=1 Tmax=%d T0Mode=%s Nt0=%d "
            "T0Fixed=%d T0Start=%d T0Stride=%d "
            "NeededTimes=%d temporal_bc=periodic "
            "rho=%s ProfileBasis=%s "
            "NProfiles=%d "
            "weight=%s\n",
            nc,
            NV,
            rmax,
            tmax,
            t0_mode_name(t0_mode),
            nt0,
            t0_fixed,
            t0_start,
            t0_stride,
            n_needed_times,
            profile_rho_name(profile_basis),
            profile_basis_name(profile_basis),
            nprofiles,
            profile_weight_name(profile_basis)
        );

        for (int k = 0; k < nprofiles; k++) {
            double minimum = weights[k][0][0];
            double maximum = weights[k][0][0];

            for (int t = 0; t < NT; t++) {
                for (int i = 0; i < NV; i++) {
                    const double value = weights[k][t][i];

                    if (value < minimum) {
                        minimum = value;
                    }

                    if (value > maximum) {
                        maximum = value;
                    }
                }
            }

            if (profile_basis == PROFILE_BASIS_DELTA) {
                printf(
                    "PROFILE %d delta_mode=%d "
                    "%.16e %.16e\n",
                    k,
                    k,
                    minimum,
                    maximum
                );
            } else {
                printf(
                    "PROFILE %d %.4f %.16e %.16e\n",
                    k,
                    paper_profile_sigmas[k],
                    minimum,
                    maximum
                );
            }
        }

        printf("T0_LIST");

        for (int source = 0; source < nt0; source++) {
            printf(" %d", t0_values[source]);
        }

        printf("\n");

        printf(
            "# DATA cfg T R Nsrc "
            "Re[avg L_constant(R,T)] "
            "Im[avg L_constant(R,T)]\n"
        );

        printf(
            "# WDATA cfg T R k l Nsrc "
            "Re[avg L_kl(R,T)] "
            "Im[avg L_kl(R,T)]\n"
        );

        fflush(stdout);
    }

    for (int T = 1; T <= tmax; T++) {
        for (int R = 0; R <= rmax; R++) {
            qcd_complex_16 Lsum = {0.0, 0.0};
            qcd_complex_16 Lsum_weighted
                [MAX_PROFILES][MAX_PROFILES];

            memset(
                Lsum_weighted,
                0,
                sizeof(Lsum_weighted)
            );

            long nsrc = 0;

            for (
                int source = 0;
                source < nt0;
                source++
            ) {
                const int t0 = t0_values[source];
                const int t1 = (t0 + T) % L[0];

                const qcd_spinorComponent3d *v_src = (
                    &v_by_time[(size_t)t0 * NV]
                );

                const qcd_spinorComponent3d *v_sink = (
                    &v_by_time[(size_t)t1 * NV]
                );

                for (int xx = 0; xx < L[1]; xx++) {
                    for (int yy = 0; yy < L[2]; yy++) {
                        for (int zz = 0; zz < L[3]; zz++) {
                            const int yx = (
                                xx + R
                            ) % L[1];

                            const int x_site3 = site3_index(
                                xx,
                                yy,
                                zz,
                                L
                            );

                            const int y_site3 = site3_index(
                                yx,
                                yy,
                                zz,
                                L
                            );

                            qcd_complex_16 Lxy = {
                                0.0,
                                0.0
                            };

                            qcd_complex_16 pair_matrix[NV][NV];

                            for (int i = 0; i < NV; i++) {
                                for (int j = 0; j < NV; j++) {
                                    const qcd_complex_16 tau_y =
                                        tau_time_path(
                                            v_src,
                                            v_sink,
                                            &u,
                                            i,
                                            j,
                                            t0,
                                            T,
                                            y_site3,
                                            yx,
                                            yy,
                                            zz,
                                            L
                                        );

                                    const qcd_complex_16 tau_x =
                                        tau_time_path(
                                            v_src,
                                            v_sink,
                                            &u,
                                            i,
                                            j,
                                            t0,
                                            T,
                                            x_site3,
                                            xx,
                                            yy,
                                            zz,
                                            L
                                        );

                                    const qcd_complex_16 pair =
                                        qcd_CMUL(
                                            tau_y,
                                            qcd_CONJ(tau_x)
                                        );

                                    pair_matrix[i][j] = pair;

                                    Lxy = qcd_CADD(
                                        Lxy,
                                        pair
                                    );
                                }
                            }

                            /*
                              Equation (8):

                                L_kl = sum_ij
                                  w_i^(k)(t0)
                                  A_ij
                                  w_j^(l)(t1),

                              with A_ij = tau_ij(y)
                              tau_ji(x) and
                              w = rho^2.  The two-stage
                              contraction uses the separable
                              endpoint weights and avoids a
                              profile-matrix update inside
                              every i,j pair.
                            */
                            qcd_complex_16 left
                                [MAX_PROFILES][NV];
                            qcd_complex_16 Lxy_weighted
                                [MAX_PROFILES][MAX_PROFILES];

                            memset(
                                left,
                                0,
                                sizeof(left)
                            );

                            memset(
                                Lxy_weighted,
                                0,
                                sizeof(Lxy_weighted)
                            );

                            for (
                                int k = 0;
                                k < nprofiles;
                                k++
                            ) {
                                for (int j = 0; j < NV; j++) {
                                    for (
                                        int i = 0;
                                        i < NV;
                                        i++
                                    ) {
                                        add_scaled_complex(
                                            &left[k][j],
                                            pair_matrix[i][j],
                                            weights[k][t0][i]
                                        );
                                    }
                                }
                            }

                            for (
                                int k = 0;
                                k < nprofiles;
                                k++
                            ) {
                                for (
                                    int l = 0;
                                    l < nprofiles;
                                    l++
                                ) {
                                    for (
                                        int j = 0;
                                        j < NV;
                                        j++
                                    ) {
                                        add_scaled_complex(
                                            &Lxy_weighted[k][l],
                                            left[k][j],
                                            weights[l][t1][j]
                                        );
                                    }
                                }
                            }

                            Lsum = qcd_CADD(
                                Lsum,
                                Lxy
                            );

                            for (
                                int k = 0;
                                k < nprofiles;
                                k++
                            ) {
                                for (
                                    int l = 0;
                                    l < nprofiles;
                                    l++
                                ) {
                                    Lsum_weighted[k][l] =
                                        qcd_CADD(
                                            Lsum_weighted[k][l],
                                            Lxy_weighted[k][l]
                                        );
                                }
                            }

                            nsrc++;
                        }
                    }
                }
            }

            const qcd_complex_16 Lavg = {
                Lsum.re / (double)nsrc,
                Lsum.im / (double)nsrc
            };

            if (myid == 0) {
                printf(
                    "DATA %d %d %d %ld %+.16e %+.16e\n",
                    nc,
                    T,
                    R,
                    nsrc,
                    Lavg.re,
                    Lavg.im
                );

                for (int k = 0; k < nprofiles; k++) {
                    for (
                        int l = 0;
                        l < nprofiles;
                        l++
                    ) {
                        const qcd_complex_16 Lavg_weighted = {
                            Lsum_weighted[k][l].re /
                                (double)nsrc,
                            Lsum_weighted[k][l].im /
                                (double)nsrc
                        };

                        printf(
                            "WDATA %d %d %d %d %d %ld "
                            "%+.16e %+.16e\n",
                            nc,
                            T,
                            R,
                            k,
                            l,
                            nsrc,
                            Lavg_weighted.re,
                            Lavg_weighted.im
                        );
                    }
                }

                fflush(stdout);
            }
        }
    }

    for (int t = 0; t < L[0]; t++) {
        if (!initialized_time[t]) {
            continue;
        }

        qcd_spinorComponent3d *v_time = (
            &v_by_time[(size_t)t * NV]
        );

        for (int i = 0; i < NV; i++) {
            qcd_destroySpinorComponent3d(
                &v_time[i]
            );
        }
    }

    free(v_by_time);

    qcd_destroyGaugeField(&u);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();

    return 0;
}
