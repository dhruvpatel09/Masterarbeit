#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <hdf5.h>
#include <qcd.h>

#define TINY_NT 8
#define TINY_NS 4
#define TINY_CFG 1
#define TINY_AVAILABLE_NV 16
#define DEFAULT_NV 16
#define DEFAULT_T_SRC 2
#define DEFAULT_T_SINK 6
#define DEFAULT_R_SEP 2
#define PATH_BUFFER_SIZE 2048

static int site3_index(const int x, const int y, const int z, const int L[4]) {
    return (x * L[2] + y) * L[3] + z;
}

static int site4_index(
    const int t,
    const int x,
    const int y,
    const int z,
    const int L[4]
) {
    return ((t * L[1] + x) * L[2] + y) * L[3] + z;
}

static size_t tau_index(
    const int site,
    const int i,
    const int j,
    const int nvecs
) {
    return (((size_t)site * (size_t)nvecs + (size_t)i) *
            (size_t)nvecs + (size_t)j);
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

static void check_readable_file(const char *path) {
    if (access(path, R_OK) != 0) {
        fprintf(
            stderr,
            "ERROR: required input is not readable: %s (%s)\n",
            path,
            strerror(errno)
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

static void format_path(
    char *output,
    const size_t output_size,
    const char *format,
    const char *root,
    const int number
) {
    const int written = snprintf(output, output_size, format, root, number);

    if (written < 0 || (size_t)written >= output_size) {
        fprintf(stderr, "ERROR: input path exceeds the path buffer\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

typedef struct {
    double r;
    double i;
} tiny_h5_complex;

static void hdf5_abort(const char *message, const char *fname) {
    fprintf(stderr, "ERROR: %s: %s\n", message, fname);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    exit(EXIT_FAILURE);
}

static int hdf5_read_scalar_int(
    const hid_t file,
    const char *dataset_name,
    const char *fname
) {
    const hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    if (dataset < 0) {
        hdf5_abort("failed to open integer dataset", fname);
    }

    const hid_t type = H5Dget_type(dataset);
    const hid_t space = H5Dget_space(dataset);
    if (type < 0 ||
        H5Tget_class(type) != H5T_INTEGER ||
        space < 0 ||
        H5Sget_simple_extent_type(space) != H5S_SCALAR) {
        hdf5_abort("invalid scalar integer dataset", fname);
    }

    int value = 0;
    if (H5Dread(
            dataset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            &value
        ) < 0) {
        hdf5_abort("failed to read scalar integer dataset", fname);
    }

    H5Sclose(space);
    H5Tclose(type);
    H5Dclose(dataset);
    return value;
}

static unsigned int hdf5_read_scalar_uint(
    const hid_t file,
    const char *dataset_name,
    const char *fname
) {
    const hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    if (dataset < 0) {
        hdf5_abort("failed to open unsigned-integer dataset", fname);
    }

    const hid_t type = H5Dget_type(dataset);
    const hid_t space = H5Dget_space(dataset);
    if (type < 0 ||
        H5Tget_class(type) != H5T_INTEGER ||
        space < 0 ||
        H5Sget_simple_extent_type(space) != H5S_SCALAR) {
        hdf5_abort("invalid scalar unsigned-integer dataset", fname);
    }

    unsigned int value = 0;
    if (H5Dread(
            dataset,
            H5T_NATIVE_UINT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            &value
        ) < 0) {
        hdf5_abort("failed to read scalar unsigned-integer dataset", fname);
    }

    H5Sclose(space);
    H5Tclose(type);
    H5Dclose(dataset);
    return value;
}

static void hdf5_read_spatial_size(
    const hid_t file,
    int spatial_size[3],
    const char *fname
) {
    const hid_t dataset = H5Dopen2(
        file, "/spatial_lat_size", H5P_DEFAULT
    );
    if (dataset < 0) {
        hdf5_abort("failed to open spatial-size dataset", fname);
    }

    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[1] = {0};
    if (space < 0 ||
        H5Sget_simple_extent_ndims(space) != 1 ||
        H5Sget_simple_extent_dims(space, dims, NULL) < 0 ||
        dims[0] != 3) {
        hdf5_abort("invalid spatial-size dataset", fname);
    }

    if (H5Dread(
            dataset,
            H5T_NATIVE_INT,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            spatial_size
        ) < 0) {
        hdf5_abort("failed to read spatial-size dataset", fname);
    }

    H5Sclose(space);
    H5Dclose(dataset);
}

static char *hdf5_read_scalar_string(
    const hid_t file,
    const char *dataset_name,
    const char *fname
) {
    const hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    if (dataset < 0) {
        hdf5_abort("failed to open string dataset", fname);
    }

    const hid_t file_type = H5Dget_type(dataset);
    const hid_t space = H5Dget_space(dataset);
    if (file_type < 0 ||
        H5Tget_class(file_type) != H5T_STRING ||
        H5Tis_variable_str(file_type) > 0 ||
        space < 0 ||
        H5Sget_simple_extent_type(space) != H5S_SCALAR) {
        hdf5_abort("invalid fixed-length scalar string dataset", fname);
    }

    const size_t size = H5Tget_size(file_type);
    char *value = calloc(size + 1U, sizeof(*value));
    if (value == NULL) {
        hdf5_abort("string allocation failed", fname);
    }

    const hid_t memory_type = H5Tcopy(H5T_C_S1);
    if (memory_type < 0 ||
        H5Tset_size(memory_type, size) < 0 ||
        H5Tset_strpad(memory_type, H5Tget_strpad(file_type)) < 0) {
        hdf5_abort("failed to create string memory type", fname);
    }

    if (H5Dread(
            dataset,
            memory_type,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            value
        ) < 0) {
        hdf5_abort("failed to read string dataset", fname);
    }
    value[size] = '\0';

    H5Tclose(memory_type);
    H5Sclose(space);
    H5Tclose(file_type);
    H5Dclose(dataset);
    return value;
}

static void read_evec_time_serial_hdf5(
    const char *fname,
    const int nc,
    const int time,
    const int nvecs,
    const char *nbase,
    qcd_geometry *geo,
    qcd_spinorComponent3d *v
) {
    const hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        hdf5_abort("failed to open HDF5 eigenvector file", fname);
    }

    const int tmp_nc = hdf5_read_scalar_int(file, "/cnfg_id", fname);
    const int tmp_time = hdf5_read_scalar_int(file, "/time", fname);
    const unsigned int tmp_nv = hdf5_read_scalar_uint(
        file, "/disteigpairs_num", fname
    );
    char *tmp_run_name = hdf5_read_scalar_string(
        file, "/run_name", fname
    );
    int spatial_size[3] = {0, 0, 0};
    hdf5_read_spatial_size(file, spatial_size, fname);

    if (tmp_nc != nc ||
        tmp_time != time ||
        tmp_nv < (unsigned int)nvecs ||
        strcmp(tmp_run_name, nbase) != 0 ||
        spatial_size[0] != geo->L[1] ||
        spatial_size[1] != geo->L[2] ||
        spatial_size[2] != geo->L[3]) {
        fprintf(
            stderr,
            "ERROR: eigenvector metadata mismatch for %s: "
            "tmp_nc=%d tmp_time=%d tmp_nv=%u tmp_run_name=%s "
            "spatial_size=%d,%d,%d\n",
            fname,
            tmp_nc,
            tmp_time,
            tmp_nv,
            tmp_run_name,
            spatial_size[0],
            spatial_size[1],
            spatial_size[2]
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    free(tmp_run_name);

    const hid_t dataset = H5Dopen2(file, "/disteigvecs", H5P_DEFAULT);
    if (dataset < 0) {
        hdf5_abort("failed to open eigenvector dataset", fname);
    }

    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[5] = {0, 0, 0, 0, 0};
    if (space < 0 ||
        H5Sget_simple_extent_ndims(space) != 5 ||
        H5Sget_simple_extent_dims(space, dims, NULL) < 0 ||
        dims[0] != tmp_nv ||
        dims[1] != (hsize_t)geo->L[1] ||
        dims[2] != (hsize_t)geo->L[2] ||
        dims[3] != (hsize_t)geo->L[3] ||
        dims[4] != 3) {
        hdf5_abort("unexpected eigenvector dataset shape", fname);
    }

    const size_t available_count = (
        (size_t)tmp_nv * (size_t)geo->lV3 * 3U
    );
    tiny_h5_complex *raw = calloc(available_count, sizeof(*raw));
    if (raw == NULL) {
        hdf5_abort("eigenvector allocation failed", fname);
    }

    const hid_t memory_type = H5Tcreate(H5T_COMPOUND, sizeof(*raw));
    if (memory_type < 0 ||
        H5Tinsert(
            memory_type,
            "r",
            HOFFSET(tiny_h5_complex, r),
            H5T_NATIVE_DOUBLE
        ) < 0 ||
        H5Tinsert(
            memory_type,
            "i",
            HOFFSET(tiny_h5_complex, i),
            H5T_NATIVE_DOUBLE
        ) < 0) {
        hdf5_abort("failed to create compound complex memory type", fname);
    }

    if (H5Dread(
            dataset,
            memory_type,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            raw
        ) < 0) {
        hdf5_abort("failed to read eigenvector dataset", fname);
    }

    for (int mode = 0; mode < nvecs; mode++) {
        for (long site = 0; site < geo->lV3; site++) {
            for (int color = 0; color < 3; color++) {
                const size_t index = (
                    ((size_t)mode * (size_t)geo->lV3 + (size_t)site) *
                    3U + (size_t)color
                );
                v[mode].D[site][color].re = raw[index].r;
                v[mode].D[site][color].im = raw[index].i;
            }
        }
    }

    free(raw);
    H5Tclose(memory_type);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
}

static qcd_complex_16 tau_time_path_open(
    const qcd_spinorComponent3d *v_src,
    const qcd_spinorComponent3d *v_sink,
    const qcd_gaugeField *u,
    const int i,
    const int j,
    const int t_src,
    const int t_sink,
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
      Open temporal path with no modulo operation:

        U_0(x,t_src) ... U_0(x,t_sink-1) v_j(x,t_sink).

      Applying links from the sink backwards produces the ordered product.
    */
    for (int tt = t_sink - 1; tt >= t_src; tt--) {
        const int site4 = site4_index(tt, x, y, z, L);

        for (int a = 0; a < 3; a++) {
            tmp[a].re = 0.0;
            tmp[a].im = 0.0;

            for (int b = 0; b < 3; b++) {
                tmp[a] = qcd_CADD(
                    tmp[a],
                    qcd_CMUL(u->D[site4][0][a][b], vec[b])
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
            qcd_CMUL(qcd_CONJ(v_src[i].D[site3][a]), vec[a])
        );
    }

    return tau;
}

static qcd_complex_16 pair_correlator_origin_adjointed(
    const qcd_complex_16 *tau,
    const int origin_site,
    const int displaced_site,
    const int nvecs
) {
    qcd_complex_16 value = {0.0, 0.0};

    /*
      The spatial line at local (0,0,0) is the adjointed line and the
      line at local (R,0,0) is unadjointed:

        L = sum_ij tau_ij(displaced) tau_ij(origin)^*.

      This naming makes the endpoint convention explicit.  It is
      mathematically the same ordering used by the validated L(R,tau)
      calculation.
    */
    for (int i = 0; i < nvecs; i++) {
        for (int j = 0; j < nvecs; j++) {
            const qcd_complex_16 tau_displaced = (
                tau[tau_index(displaced_site, i, j, nvecs)]
            );
            const qcd_complex_16 tau_origin = (
                tau[tau_index(origin_site, i, j, nvecs)]
            );

            value = qcd_CADD(
                value,
                qcd_CMUL(tau_displaced, qcd_CONJ(tau_origin))
            );
        }
    }

    return value;
}


static int wrap_mod(const int value, const int extent) {
    const int remainder = value % extent;
    return (remainder < 0) ? remainder + extent : remainder;
}

typedef struct {
    double plane_action[6];
    double S_raw;
} tiny_clover_density;

static double raw_plane_action_from_Q(qcd_complex_16 Q[3][3]) {
    qcd_complex_16 F[3][3];

    /*
      qcdlib clover normalization used in this cross-check:

        F_munu = (Q_munu - Q_munu^\dagger) / 8,

      with no color-trace subtraction.  Since F is anti-Hermitian,

        -Tr(F F) = Tr(F^\dagger F) = sum_ab |F_ab|^2.
    */
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            F[a][b].re = 0.125 * (Q[a][b].re - Q[b][a].re);
            F[a][b].im = 0.125 * (Q[a][b].im + Q[b][a].im);
        }
    }

    double plane_action = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            plane_action += (
                F[a][b].re * F[a][b].re +
                F[a][b].im * F[a][b].im
            );
        }
    }

    return plane_action;
}

static tiny_clover_density local_clover_density(
    const qcd_cloverField *clover,
    const int site4
) {
    tiny_clover_density density;
    memset(&density, 0, sizeof(density));

    /*
      qcd_cloverField plane order:
        0: 01, 1: 02, 2: 03, 3: 12, 4: 13, 5: 23.

      The observable is the complete six-plane action density;
      no electric/magnetic separation is made here.
    */
    for (int plane = 0; plane < 6; plane++) {
        density.plane_action[plane] = raw_plane_action_from_Q(
            clover->D[site4][plane]
        );
        density.S_raw += density.plane_action[plane];
    }

    return density;
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
                "ERROR: tiny Laplace-clover validation requires exactly one MPI rank\n"
            );
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int L[4] = {TINY_NT, TINY_NS, TINY_NS, TINY_NS};
    int P[4] = {1, 1, 1, 1};
    double theta[3] = {0.0, 0.0, 0.0};

    const int nvecs = read_env_int(
        "NVECS", DEFAULT_NV, 1, TINY_AVAILABLE_NV
    );
    const int t_src = read_env_int(
        "T_SRC", DEFAULT_T_SRC, 1, L[0] - 2
    );
    const int t_sink = read_env_int(
        "T_SINK", DEFAULT_T_SINK, 1, L[0] - 2
    );
    const int r_sep = read_env_int(
        "R_SEP", DEFAULT_R_SEP, 0, L[1] / 2
    );
    const int point_x = read_env_int("POINT_X", 0, 0, L[1] - 1);
    const int point_y = read_env_int("POINT_Y", 1, 0, L[2] - 1);
    const int point_z = read_env_int("POINT_Z", 2, 0, L[3] - 1);
    const int probe_long = read_env_int(
        "PROBE_LONG", 0, -L[1], L[1]
    );
    const int probe_tr1 = read_env_int(
        "PROBE_TR1", 1, -L[2], L[2]
    );
    const int probe_tr2 = read_env_int(
        "PROBE_TR2", 2, -L[3], L[3]
    );
    const int print_site_values = read_env_int(
        "PRINT_SITE_VALUES", 0, 0, 1
    );

    if (t_sink <= t_src) {
        fprintf(
            stderr,
            "ERROR: require 0 < T_SRC < T_SINK < %d for open boundaries\n",
            L[0] - 1
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if ((t_sink - t_src) % 2 != 0) {
        fprintf(
            stderr,
            "ERROR: an even temporal separation is required for an integer midpoint\n"
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const int tau_extent = t_sink - t_src;
    const int t_mid = t_src + tau_extent / 2;
    const char *nbase = "open8x4x4x4";
    const char *data_root = getenv("TINY_DIST_ROOT");

    if (data_root == NULL || *data_root == '\0') {
        data_root = (
            "/home/m2130292/Masterarbeit/laplace_static/validation/"
            "tiny_open8x4x4x4/data"
        );
    }

    char gauge_file[PATH_BUFFER_SIZE];
    char evec_src_file[PATH_BUFFER_SIZE];
    char evec_sink_file[PATH_BUFFER_SIZE];

    format_path(
        gauge_file,
        sizeof(gauge_file),
        "%s/cnfg/open8x4x4x4n%d",
        data_root,
        TINY_CFG
    );
    format_path(
        evec_src_file,
        sizeof(evec_src_file),
        "%s/disteigvecs/disteigvecs_open8x4x4x4n1_t%d.h5",
        data_root,
        t_src
    );
    format_path(
        evec_sink_file,
        sizeof(evec_sink_file),
        "%s/disteigvecs/disteigvecs_open8x4x4x4n1_t%d.h5",
        data_root,
        t_sink
    );

    check_readable_file(gauge_file);
    check_readable_file(evec_src_file);
    check_readable_file(evec_sink_file);

    qcd_geometry geo;
    qcd_gaugeField gauge;
    qcd_cloverField clover;

    if (qcd_initGeometry(&geo, L, P, theta, myid, numprocs)) {
        fprintf(stderr, "ERROR: qcd_initGeometry failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (qcd_initEO(&geo)) {
        fprintf(stderr, "ERROR: qcd_initEO failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (qcd_initGaugeField(&gauge, &geo)) {
        fprintf(stderr, "ERROR: qcd_initGaugeField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (qcd_initCloverField(&clover, &geo)) {
        fprintf(stderr, "ERROR: qcd_initCloverField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const long spatial_volume = (long)L[1] * L[2] * L[3];
    if (geo.lV3 != spatial_volume) {
        fprintf(
            stderr,
            "ERROR: serial validation expected lV3=%ld, got %ld\n",
            spatial_volume,
            geo.lV3
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (qcd_getGaugeField(gauge_file, qcd_GF_OPENQCD, &gauge)) {
        fprintf(
            stderr,
            "ERROR: qcd_getGaugeField failed for %s\n",
            gauge_file
        );
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const double plaquette = qcd_calculatePlaquette(&gauge);
    if (qcd_calculateCloverField(&clover, &gauge)) {
        fprintf(stderr, "ERROR: qcd_calculateCloverField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    qcd_spinorComponent3d *v_src = calloc(
        (size_t)nvecs, sizeof(*v_src)
    );
    qcd_spinorComponent3d *v_sink = calloc(
        (size_t)nvecs, sizeof(*v_sink)
    );
    if (v_src == NULL || v_sink == NULL) {
        fprintf(stderr, "ERROR: eigenvector-array allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int mode = 0; mode < nvecs; mode++) {
        if (qcd_initSpinorComponent3d(&v_src[mode], &geo) ||
            qcd_initSpinorComponent3d(&v_sink[mode], &geo)) {
            fprintf(
                stderr,
                "ERROR: qcd_initSpinorComponent3d failed at mode=%d\n",
                mode
            );
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    read_evec_time_serial_hdf5(
        evec_src_file,
        TINY_CFG,
        t_src,
        nvecs,
        nbase,
        &geo,
        v_src
    );
    read_evec_time_serial_hdf5(
        evec_sink_file,
        TINY_CFG,
        t_sink,
        nvecs,
        nbase,
        &geo,
        v_sink
    );

    const size_t tau_count = (
        (size_t)spatial_volume * (size_t)nvecs * (size_t)nvecs
    );
    qcd_complex_16 *tau = calloc(tau_count, sizeof(*tau));
    if (tau == NULL) {
        fprintf(stderr, "ERROR: static-perambulator allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int x = 0; x < L[1]; x++) {
        for (int y = 0; y < L[2]; y++) {
            for (int z = 0; z < L[3]; z++) {
                const int site = site3_index(x, y, z, L);

                for (int i = 0; i < nvecs; i++) {
                    for (int j = 0; j < nvecs; j++) {
                        tau[tau_index(site, i, j, nvecs)] =
                            tau_time_path_open(
                                v_src,
                                v_sink,
                                &gauge,
                                i,
                                j,
                                t_src,
                                t_sink,
                                site,
                                x,
                                y,
                                z,
                                L
                            );
                    }
                }
            }
        }
    }

    const int point_site4 = site4_index(
        t_mid, point_x, point_y, point_z, L
    );
    const tiny_clover_density point_density = local_clover_density(
        &clover, point_site4
    );
    const char *plane_name[6] = {
        "01", "02", "03", "12", "13", "23"
    };
    const double expected_point_plane[6] = {
        0.1421358974693371,
        0.2229753324656666,
        0.3360429657631333,
        0.3473193984043027,
        0.2150515455120535,
        0.2047295356401899
    };
    const double expected_point_S_raw = 1.4682546752546830;
    const int point_regression_applicable = (
        t_mid == 4 &&
        point_x == 0 &&
        point_y == 1 &&
        point_z == 2
    );
    double max_point_plane_absdiff = 0.0;
    for (int plane = 0; plane < 6; plane++) {
        const double absdiff = fabs(
            point_density.plane_action[plane] -
            expected_point_plane[plane]
        );
        if (absdiff > max_point_plane_absdiff) {
            max_point_plane_absdiff = absdiff;
        }
    }
    const double point_S_raw_absdiff = fabs(
        point_density.S_raw - expected_point_S_raw
    );
    const double point_regression_tolerance = 1.0e-14;
    const int point_regression_ok = (
        !point_regression_applicable ||
        (max_point_plane_absdiff <= point_regression_tolerance &&
         point_S_raw_absdiff <= point_regression_tolerance)
    );

    if (myid == 0) {
        printf(
            "META cfg=%d Nt=%d Ns=%d Nv=%d r=%d "
            "t_src=%d t_sink=%d tau=%d t_mid=%d "
            "temporal_bc=open spatial_bc=periodic rho=constant "
            "evec_reader=serial_hdf5 clover=qcdlib\n",
            TINY_CFG,
            L[0],
            L[1],
            nvecs,
            r_sep,
            t_src,
            t_sink,
            tau_extent,
            t_mid
        );
        printf("INPUT gauge=%s\n", gauge_file);
        printf("INPUT evec_src=%s\n", evec_src_file);
        printf("INPUT evec_sink=%s\n", evec_sink_file);
        printf("PATH temporal_links");
        for (int time = t_src; time < t_sink; time++) {
            printf(" %d", time);
        }
        printf("\n");
        printf(
            "CLOVER_CONVENTION F=(Q-Qdagger)/8 color=unprojected "
            "plane_action=-Tr(FF)=Tr(FdaggerF) "
            "S_raw=sum_mu_lt_nu_plane_action\n"
        );
        printf("PLAQUETTE value=%+.16e\n", plaquette);
        for (int plane = 0; plane < 6; plane++) {
            printf(
                "CLOVER_PLANE t=%d x=%d y=%d z=%d plane=%s "
                "minus_Tr_FF=%+.16e\n",
                t_mid,
                point_x,
                point_y,
                point_z,
                plane_name[plane],
                point_density.plane_action[plane]
            );
        }
        printf(
            "CLOVER_POINT t=%d x=%d y=%d z=%d "
            "S_raw=-sum_mu_lt_nu_Tr_FF=%+.16e\n",
            t_mid,
            point_x,
            point_y,
            point_z,
            point_density.S_raw
        );
        printf(
            "CLOVER_REFERENCE applicable=%d "
            "expected_S_raw=%+.16e S_raw_absdiff=%+.16e "
            "max_plane_absdiff=%+.16e tolerance=%+.16e ok=%d\n",
            point_regression_applicable,
            expected_point_S_raw,
            point_S_raw_absdiff,
            max_point_plane_absdiff,
            point_regression_tolerance,
            point_regression_ok
        );
        printf(
            "PROBE_DEF frame=source_relative "
            "longitudinal=%d transverse1=%d transverse2=%d "
            "rotation=cyclic_with_axis\n",
            probe_long,
            probe_tr1,
            probe_tr2
        );
        printf(
            "STATIC_LINE_DEF local_origin=0,0,0 origin_line=adjointed "
            "local_displaced=%d,0,0 displaced_line=unadjointed "
            "product=displaced_times_conjugate_origin\n",
            r_sep
        );
    }

    const char *axis_name[3] = {"x", "y", "z"};
    qcd_complex_16 combined_L_sum = {0.0, 0.0};
    qcd_complex_16 combined_LS_sum = {0.0, 0.0};
    double combined_S_raw_sum = 0.0;
    long combined_nsrc = 0;

    double axis_S_raw_avg[3] = {0.0, 0.0, 0.0};
    double max_pair_conjugacy_abs = 0.0;
    double max_axis_L_abs_im = 0.0;

    for (int axis = 0; axis < 3; axis++) {
        const int transverse1 = (axis + 1) % 3;
        const int transverse2 = (axis + 2) % 3;
        int global_probe_offset[3] = {0, 0, 0};
        global_probe_offset[axis] = probe_long;
        global_probe_offset[transverse1] = probe_tr1;
        global_probe_offset[transverse2] = probe_tr2;

        qcd_complex_16 L_sum = {0.0, 0.0};
        qcd_complex_16 LS_sum = {0.0, 0.0};
        qcd_complex_16 origin_pair = {0.0, 0.0};
        double S_raw_sum = 0.0;
        long axis_nsrc = 0;

        for (int x = 0; x < L[1]; x++) {
            for (int y = 0; y < L[2]; y++) {
                for (int z = 0; z < L[3]; z++) {
                    const int source_coord[3] = {x, y, z};
                    int sink_coord[3] = {x, y, z};
                    int probe_coord[3] = {x, y, z};

                    sink_coord[axis] = wrap_mod(
                        source_coord[axis] + r_sep,
                        L[axis + 1]
                    );
                    probe_coord[axis] = wrap_mod(
                        source_coord[axis] + probe_long,
                        L[axis + 1]
                    );
                    probe_coord[transverse1] = wrap_mod(
                        source_coord[transverse1] + probe_tr1,
                        L[transverse1 + 1]
                    );
                    probe_coord[transverse2] = wrap_mod(
                        source_coord[transverse2] + probe_tr2,
                        L[transverse2 + 1]
                    );

                    const int source_site = site3_index(x, y, z, L);
                    const int sink_site = site3_index(
                        sink_coord[0],
                        sink_coord[1],
                        sink_coord[2],
                        L
                    );
                    const int probe_site4 = site4_index(
                        t_mid,
                        probe_coord[0],
                        probe_coord[1],
                        probe_coord[2],
                        L
                    );

                    const qcd_complex_16 L_value =
                        pair_correlator_origin_adjointed(
                            tau, source_site, sink_site, nvecs
                        );
                    const qcd_complex_16 reverse =
                        pair_correlator_origin_adjointed(
                            tau, sink_site, source_site, nvecs
                        );
                    const tiny_clover_density density =
                        local_clover_density(&clover, probe_site4);

                    L_sum = qcd_CADD(L_sum, L_value);
                    S_raw_sum += density.S_raw;
                    LS_sum.re += L_value.re * density.S_raw;
                    LS_sum.im += L_value.im * density.S_raw;
                    axis_nsrc++;

                    const double residual_re = (
                        L_value.re - reverse.re
                    );
                    const double residual_im = (
                        L_value.im + reverse.im
                    );
                    const double residual_abs = hypot(
                        residual_re, residual_im
                    );
                    if (residual_abs > max_pair_conjugacy_abs) {
                        max_pair_conjugacy_abs = residual_abs;
                    }

                    if (x == 0 && y == 0 && z == 0) {
                        origin_pair = L_value;
                    }

                    if (print_site_values && myid == 0) {
                        printf(
                            "SITE axis=%s source=%d,%d,%d sink=%d,%d,%d "
                            "probe_t=%d probe=%d,%d,%d "
                            "L_Re=%+.16e L_Im=%+.16e "
                            "S_raw=%+.16e\n",
                            axis_name[axis],
                            x,
                            y,
                            z,
                            sink_coord[0],
                            sink_coord[1],
                            sink_coord[2],
                            t_mid,
                            probe_coord[0],
                            probe_coord[1],
                            probe_coord[2],
                            L_value.re,
                            L_value.im,
                            density.S_raw
                        );
                    }
                }
            }
        }

        const double inverse_nsrc = 1.0 / (double)axis_nsrc;
        const qcd_complex_16 L_avg = {
            L_sum.re * inverse_nsrc,
            L_sum.im * inverse_nsrc
        };
        const qcd_complex_16 LS_avg = {
            LS_sum.re * inverse_nsrc,
            LS_sum.im * inverse_nsrc
        };
        const double S_raw_avg = S_raw_sum * inverse_nsrc;

        axis_S_raw_avg[axis] = S_raw_avg;
        if (fabs(L_avg.im) > max_axis_L_abs_im) {
            max_axis_L_abs_im = fabs(L_avg.im);
        }

        combined_L_sum = qcd_CADD(combined_L_sum, L_sum);
        combined_LS_sum = qcd_CADD(combined_LS_sum, LS_sum);
        combined_S_raw_sum += S_raw_sum;
        combined_nsrc += axis_nsrc;

        if (myid == 0) {
            printf(
                "PAIR spacetime_plane=(%d,0) axis=%s "
                "local_origin=0,0,0 origin_line=adjointed "
                "local_displaced=%d,0,0 displaced_line=unadjointed "
                "Re=%+.16e Im=%+.16e\n",
                axis + 1,
                axis_name[axis],
                r_sep,
                origin_pair.re,
                origin_pair.im
            );
            printf(
                "AXIS_LS spacetime_plane=(%d,0) axis=%s Nsrc=%ld "
                "local_probe=%d,%d,%d global_probe_offset=%d,%d,%d "
                "L_SUM_Re=%+.16e L_SUM_Im=%+.16e "
                "L_AVG_Re=%+.16e L_AVG_Im=%+.16e "
                "S_RAW_SUM=%+.16e S_RAW_AVG=%+.16e "
                "LS_SUM_Re=%+.16e LS_SUM_Im=%+.16e "
                "LS_AVG_Re=%+.16e LS_AVG_Im=%+.16e\n",
                axis + 1,
                axis_name[axis],
                axis_nsrc,
                probe_long,
                probe_tr1,
                probe_tr2,
                global_probe_offset[0],
                global_probe_offset[1],
                global_probe_offset[2],
                L_sum.re,
                L_sum.im,
                L_avg.re,
                L_avg.im,
                S_raw_sum,
                S_raw_avg,
                LS_sum.re,
                LS_sum.im,
                LS_avg.re,
                LS_avg.im
            );
        }
    }

    const double inverse_combined_nsrc = (
        1.0 / (double)combined_nsrc
    );
    const qcd_complex_16 combined_L_avg = {
        combined_L_sum.re * inverse_combined_nsrc,
        combined_L_sum.im * inverse_combined_nsrc
    };
    const qcd_complex_16 combined_LS_avg = {
        combined_LS_sum.re * inverse_combined_nsrc,
        combined_LS_sum.im * inverse_combined_nsrc
    };
    const double combined_S_raw_avg = (
        combined_S_raw_sum * inverse_combined_nsrc
    );

    double max_axis_S_raw_spread = 0.0;
    for (int first = 0; first < 3; first++) {
        for (int second = first + 1; second < 3; second++) {
            const double spread = fabs(
                axis_S_raw_avg[first] - axis_S_raw_avg[second]
            );
            if (spread > max_axis_S_raw_spread) {
                max_axis_S_raw_spread = spread;
            }
        }
    }

    const double expected_L_nvec16 = 5.2853688927290924e-04;
    const int regression_applicable = (
        nvecs == 16 &&
        t_src == 2 &&
        t_sink == 6 &&
        r_sep == 2
    );
    const double L_reference_absdiff = fabs(
        combined_L_avg.re - expected_L_nvec16
    );
    const double L_reference_reldiff = (
        L_reference_absdiff / fabs(expected_L_nvec16)
    );
    const double L_regression_tolerance = 1.0e-15;
    const int L_regression_ok = (
        !regression_applicable ||
        L_reference_absdiff <= L_regression_tolerance
    );

    if (myid == 0) {
        printf(
            "COMBINED_DIAGNOSTIC Nsrc=%ld "
            "L_Re=%+.16e L_Im=%+.16e "
            "S_RAW_AVG=%+.16e "
            "LS_Re=%+.16e LS_Im=%+.16e\n",
            combined_nsrc,
            combined_L_avg.re,
            combined_L_avg.im,
            combined_S_raw_avg,
            combined_LS_avg.re,
            combined_LS_avg.im
        );
        printf(
            "CHECK max_pair_conjugacy_abs=%+.16e "
            "max_axis_L_abs_im=%+.16e "
            "max_axis_S_raw_spread=%+.16e "
            "L_reference_applicable=%d "
            "L_reference_expected=%+.16e "
            "L_reference_absdiff=%+.16e "
            "L_reference_reldiff=%+.16e "
            "L_regression_tolerance=%+.16e "
            "L_regression_ok=%d\n",
            max_pair_conjugacy_abs,
            max_axis_L_abs_im,
            max_axis_S_raw_spread,
            regression_applicable,
            expected_L_nvec16,
            L_reference_absdiff,
            L_reference_reldiff,
            L_regression_tolerance,
            L_regression_ok
        );
        fflush(stdout);
    }

    free(tau);
    for (int mode = 0; mode < nvecs; mode++) {
        qcd_destroySpinorComponent3d(&v_src[mode]);
        qcd_destroySpinorComponent3d(&v_sink[mode]);
    }
    free(v_src);
    free(v_sink);

    qcd_destroyCloverField(&clover);
    qcd_destroyGaugeField(&gauge);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();

    if (!point_regression_ok) {
        fprintf(
            stderr,
            "ERROR: raw clover point regression failed: "
            "S_absdiff=%+.16e plane_max_absdiff=%+.16e "
            "tolerance=%+.16e\n",
            point_S_raw_absdiff,
            max_point_plane_absdiff,
            point_regression_tolerance
        );
        return EXIT_FAILURE;
    }

    if (!L_regression_ok) {
        fprintf(
            stderr,
            "ERROR: validated L(R,tau) regression failed: "
            "absdiff=%+.16e tolerance=%+.16e\n",
            L_reference_absdiff,
            L_regression_tolerance
        );
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
