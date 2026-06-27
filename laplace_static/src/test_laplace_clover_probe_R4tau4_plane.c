#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <mpi.h>

#include <qcd.h>
#include <DistEigvecsHdf5Reader.h>

#define NV 10

static int site3_index(const int x, const int y, const int z, const int L[4]) {
    return (x * L[2] + y) * L[3] + z;
}

static int site4_index(const int t, const int x, const int y, const int z, const int L[4]) {
    return ((t * L[1] + x) * L[2] + y) * L[3] + z;
}

static int wrap_mod(const int a, const int n) {
    int r = a % n;
    return (r < 0) ? r + n : r;
}

static int read_env_int(const char *name, const int def, const int minv, const int maxv) {
    const char *text = getenv(name);
    if (text == NULL || *text == '\0') return def;

    errno = 0;
    char *end = NULL;
    long v = strtol(text, &end, 10);

    if (errno == ERANGE || end == text || *end != '\0' || v < minv || v > maxv) {
        fprintf(stderr, "ERROR: %s must be integer in [%d,%d], got '%s'\n",
                name, minv, maxv, text);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    return (int)v;
}

typedef enum {
    T0_MODE_FIXED = 0,
    T0_MODE_AVERAGE = 1
} t0_mode_t;

static t0_mode_t read_t0_mode(void) {
    const char *text = getenv("T0_MODE");

    if (text == NULL || *text == '\0' || strcmp(text, "average") == 0) {
        return T0_MODE_AVERAGE;
    }

    if (strcmp(text, "fixed") == 0) {
        return T0_MODE_FIXED;
    }

    fprintf(stderr, "ERROR: T0_MODE must be fixed or average, got '%s'\n", text);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return T0_MODE_AVERAGE;
}

static const char *t0_mode_name(const t0_mode_t mode) {
    return mode == T0_MODE_AVERAGE ? "average" : "fixed";
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

    char *tmp_run_name = NULL;
    int tmp_nc = -1;
    int tmp_time = -1;
    unsigned int tmp_nv = 0;

    int spatial_lat_size[3] = {geo->L[1], geo->L[2], geo->L[3]};
    int spatial_proclat_size[3] = {
        geo->L[1] / geo->lL[1],
        geo->L[2] / geo->lL[2],
        geo->L[3] / geo->lL[3]
    };
    int spatial_proclat_pos[3] = {geo->Pos[1], geo->Pos[2], geo->Pos[3]};

    DistEigvecsHdf5_Metadata metadata;
    DistEigvecsHdf5_Programdata programdata;

    int ret = DistEigvecsHdf5Reader_Init(&reader);
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

    for (int step = T - 1; step >= 0; step--) {
        const int tt = (t0 + step) % L[0];
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

static double su3_norm2_traceless_F_from_Q(qcd_complex_16 Q[3][3]) {
    qcd_complex_16 F[3][3];

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            F[a][b].re = 0.125 * (Q[a][b].re - Q[b][a].re);
            F[a][b].im = 0.125 * (Q[a][b].im + Q[b][a].im);
        }
    }

    qcd_complex_16 tr = {0.0, 0.0};

    for (int a = 0; a < 3; a++) {
        tr.re += F[a][a].re;
        tr.im += F[a][a].im;
    }

    tr.re /= 3.0;
    tr.im /= 3.0;

    for (int a = 0; a < 3; a++) {
        F[a][a].re -= tr.re;
        F[a][a].im -= tr.im;
    }

    double n2 = 0.0;

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            n2 += F[a][b].re * F[a][b].re + F[a][b].im * F[a][b].im;
        }
    }

    return n2;
}

static void local_E2_B2_S(
    const qcd_cloverField *clv,
    const int site4,
    double *E2,
    double *B2,
    double *S
) {
    double e = 0.0;
    double b = 0.0;

    for (int p = 0; p < 3; p++) {
        e += su3_norm2_traceless_F_from_Q(clv->D[site4][p]);
    }

    for (int p = 3; p < 6; p++) {
        b += su3_norm2_traceless_F_from_Q(clv->D[site4][p]);
    }

    *E2 = e;
    *B2 = b;
    *S = 0.5 * (e + b);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int myid = 0;
    int numprocs = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (numprocs != 1) {
        if (myid == 0) {
            fprintf(stderr, "ERROR: this prototype is intended for exactly 1 MPI rank.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const char *nbase = "Em1p4";

    int L[4] = {48, 24, 24, 24};
    int P[4] = {1, 1, 1, 1};
    double theta[3] = {0.0, 0.0, 0.0};

    const int nc = read_env_int("NCFG", 1, 1, INT_MAX);
    const int R = read_env_int("RSEP", 4, 0, L[1] / 2);
    const int Tau = read_env_int("TAU", 4, 1, L[0] - 1);

    const int dx_min = read_env_int("DX_MIN", -2, -L[1], L[1]);
    const int dx_max = read_env_int("DX_MAX", 6, -L[1], L[1]);
    const int dy_min = read_env_int("DY_MIN", -6, -L[2], L[2]);
    const int dy_max = read_env_int("DY_MAX", 6, -L[2], L[2]);
    const int dz_fixed = read_env_int("DZ_FIXED", 0, -L[3], L[3]);

    if (dx_max < dx_min || dy_max < dy_min) {
        fprintf(stderr, "ERROR: invalid probe plane bounds\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (Tau % 2 != 0 && myid == 0) {
        fprintf(stderr, "WARNING: Tau is odd; midpoint uses integer Tau/2 floor.\n");
    }

    const int t_ins_offset = Tau / 2;

    const t0_mode_t t0_mode = read_t0_mode();
    const int t0_fixed = read_env_int("T0_FIXED", 0, 0, L[0] - 1);
    const int t0_start = read_env_int("T0_START", 0, 0, L[0] - 1);
    const int t0_stride = read_env_int("T0_STRIDE", 1, 1, L[0]);

    int t0_values[48];
    int nt0 = 0;

    if (t0_mode == T0_MODE_FIXED) {
        t0_values[nt0++] = t0_fixed;
    } else {
        for (int t0 = t0_start; t0 < L[0]; t0 += t0_stride) {
            t0_values[nt0++] = t0;
        }
    }

    if (nt0 < 1) {
        fprintf(stderr, "ERROR: no t0 values selected\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    unsigned char needed_time[48];
    unsigned char initialized_time[48];

    memset(needed_time, 0, sizeof(needed_time));
    memset(initialized_time, 0, sizeof(initialized_time));

    for (int s = 0; s < nt0; s++) {
        const int t0 = t0_values[s];
        const int t1 = (t0 + Tau) % L[0];

        needed_time[t0] = 1;
        needed_time[t1] = 1;
    }

    int n_needed_times = 0;
    for (int t = 0; t < L[0]; t++) {
        if (needed_time[t]) n_needed_times++;
    }

    const int ndx = dx_max - dx_min + 1;
    const int ndy = dy_max - dy_min + 1;
    const int nprobe = ndx * ndy;

    qcd_complex_16 *LO_E2 = calloc((size_t)nprobe, sizeof(*LO_E2));
    qcd_complex_16 *LO_B2 = calloc((size_t)nprobe, sizeof(*LO_B2));
    qcd_complex_16 *LO_S  = calloc((size_t)nprobe, sizeof(*LO_S));

    double *sum_E2 = calloc((size_t)nprobe, sizeof(*sum_E2));
    double *sum_B2 = calloc((size_t)nprobe, sizeof(*sum_B2));
    double *sum_S  = calloc((size_t)nprobe, sizeof(*sum_S));

    if (!LO_E2 || !LO_B2 || !LO_S || !sum_E2 || !sum_B2 || !sum_S) {
        fprintf(stderr, "ERROR: probe accumulator allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    qcd_geometry geo;
    qcd_gaugeField u;
    qcd_cloverField clv;

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

    if (qcd_initCloverField(&clv, &geo)) {
        fprintf(stderr, "ERROR: qcd_initCloverField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    qcd_spinorComponent3d *v_by_time = calloc((size_t)L[0] * NV, sizeof(*v_by_time));
    if (!v_by_time) {
        fprintf(stderr, "ERROR: allocation of temporal eigenvectors failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char gauge_file[1024];
    char evec_file[1024];

    snprintf(gauge_file, sizeof(gauge_file),
             "/home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n%d", nc);

    if (myid == 0) {
        printf("Laplace-clover probe plane prototype\n");
        printf("NCFG=%d R=%d Tau=%d t_ins_offset=%d axis=x\n", nc, R, Tau, t_ins_offset);
        printf("T0Mode=%s T0Fixed=%d T0Start=%d T0Stride=%d Nt0=%d NeededTimes=%d\n",
               t0_mode_name(t0_mode), t0_fixed, t0_start, t0_stride, nt0, n_needed_times);
        printf("Probe plane: dx=%d..%d dy=%d..%d dz=%d nprobe=%d\n",
               dx_min, dx_max, dy_min, dy_max, dz_fixed, nprobe);
        printf("Reading gauge field: %s\n", gauge_file);
        fflush(stdout);
    }

    if (qcd_getGaugeField(gauge_file, qcd_GF_OPENQCD, &u)) {
        fprintf(stderr, "ERROR: qcd_getGaugeField failed for %s\n", gauge_file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const double plaq = qcd_calculatePlaquette(&u);

    if (myid == 0) {
        printf("Average plaquette = %.16e\n", plaq);
        printf("Calculating clover field...\n");
        fflush(stdout);
    }

    if (qcd_calculateCloverField(&clv, &u)) {
        fprintf(stderr, "ERROR: qcd_calculateCloverField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (myid == 0) {
        printf("Clover field ready.\n");
        fflush(stdout);
    }

    for (int t = 0; t < L[0]; t++) {
        if (!needed_time[t]) continue;

        qcd_spinorComponent3d *v_time = &v_by_time[(size_t)t * NV];

        for (int i = 0; i < NV; i++) {
            if (qcd_initSpinorComponent3d(&v_time[i], &geo)) {
                fprintf(stderr, "ERROR: qcd_initSpinorComponent3d failed at t=%d i=%d\n", t, i);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }

        initialized_time[t] = 1;

        snprintf(evec_file, sizeof(evec_file),
                 "/home/m2130292/Masterarbeit/mental/"
                 "runs_Em1p4_Nv10_qcdnew_full/n%d/"
                 "eigenvectors/Em1p4n%d_evec_t%d.h5",
                 nc, nc, t);

        if (myid == 0) {
            printf("Reading temporal eigenvectors t=%d: %s\n", t, evec_file);
            fflush(stdout);
        }

        read_evec_time(evec_file, nc, t, nbase, &geo, v_time);
    }

    qcd_complex_16 Lsum = {0.0, 0.0};
    long nsrc = 0;

    for (int s = 0; s < nt0; s++) {
        const int t0 = t0_values[s];
        const int t1 = (t0 + Tau) % L[0];
        const int t_ins = (t0 + t_ins_offset) % L[0];

        const qcd_spinorComponent3d *v_src = &v_by_time[(size_t)t0 * NV];
        const qcd_spinorComponent3d *v_sink = &v_by_time[(size_t)t1 * NV];

        for (int xx = 0; xx < L[1]; xx++) {
            for (int yy = 0; yy < L[2]; yy++) {
                for (int zz = 0; zz < L[3]; zz++) {
                    const int yx = (xx + R) % L[1];

                    const int x_site3 = site3_index(xx, yy, zz, L);
                    const int y_site3 = site3_index(yx, yy, zz, L);

                    qcd_complex_16 Lxy = {0.0, 0.0};

                    for (int i = 0; i < NV; i++) {
                        for (int j = 0; j < NV; j++) {
                            const qcd_complex_16 tau_y =
                                tau_time_path(v_src, v_sink, &u, i, j,
                                              t0, Tau, y_site3, yx, yy, zz, L);

                            const qcd_complex_16 tau_x =
                                tau_time_path(v_src, v_sink, &u, i, j,
                                              t0, Tau, x_site3, xx, yy, zz, L);

                            Lxy = qcd_CADD(Lxy, qcd_CMUL(tau_y, qcd_CONJ(tau_x)));
                        }
                    }

                    Lsum = qcd_CADD(Lsum, Lxy);

                    for (int idx = 0; idx < ndx; idx++) {
                        const int dx = dx_min + idx;
                        const int xp = wrap_mod(xx + dx, L[1]);

                        for (int idy = 0; idy < ndy; idy++) {
                            const int dy = dy_min + idy;
                            const int yp = wrap_mod(yy + dy, L[2]);
                            const int zp = wrap_mod(zz + dz_fixed, L[3]);

                            const int pidx = idx * ndy + idy;
                            const int site4p = site4_index(t_ins, xp, yp, zp, L);

                            double E2 = 0.0;
                            double B2 = 0.0;
                            double Sval = 0.0;

                            local_E2_B2_S(&clv, site4p, &E2, &B2, &Sval);

                            sum_E2[pidx] += E2;
                            sum_B2[pidx] += B2;
                            sum_S[pidx]  += Sval;

                            LO_E2[pidx].re += Lxy.re * E2;
                            LO_E2[pidx].im += Lxy.im * E2;

                            LO_B2[pidx].re += Lxy.re * B2;
                            LO_B2[pidx].im += Lxy.im * B2;

                            LO_S[pidx].re += Lxy.re * Sval;
                            LO_S[pidx].im += Lxy.im * Sval;
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
        printf("LAVG cfg r tau nsrc L_Re L_Im plaquette\n");
        printf("LAVG %d %d %d %ld %+.16e %+.16e %.16e\n",
               nc, R, Tau, nsrc, Lavg.re, Lavg.im, plaq);

        printf("# PROBE cfg r tau axis dx dy dz Nsrc L_Re L_Im ");
        printf("E2_vac B2_vac S_vac ");
        printf("LE2_Re LE2_Im LB2_Re LB2_Im LS_Re LS_Im ");
        printf("rho_E2 rho_B2 rho_S\n");

        for (int idx = 0; idx < ndx; idx++) {
            const int dx = dx_min + idx;

            for (int idy = 0; idy < ndy; idy++) {
                const int dy = dy_min + idy;
                const int pidx = idx * ndy + idy;

                const double E2avg = sum_E2[pidx] / (double)nsrc;
                const double B2avg = sum_B2[pidx] / (double)nsrc;
                const double Savg  = sum_S[pidx]  / (double)nsrc;

                const qcd_complex_16 LE2avg = {
                    LO_E2[pidx].re / (double)nsrc,
                    LO_E2[pidx].im / (double)nsrc
                };

                const qcd_complex_16 LB2avg = {
                    LO_B2[pidx].re / (double)nsrc,
                    LO_B2[pidx].im / (double)nsrc
                };

                const qcd_complex_16 LSavg = {
                    LO_S[pidx].re / (double)nsrc,
                    LO_S[pidx].im / (double)nsrc
                };

                const double rho_E2 = LE2avg.re / Lavg.re - E2avg;
                const double rho_B2 = LB2avg.re / Lavg.re - B2avg;
                const double rho_S  = LSavg.re  / Lavg.re - Savg;

                printf(
                    "PROBE %d %d %d 0 %d %d %d %ld "
                    "%+.16e %+.16e "
                    "%.16e %.16e %.16e "
                    "%+.16e %+.16e %+.16e %+.16e %+.16e %+.16e "
                    "%+.16e %+.16e %+.16e\n",
                    nc, R, Tau, dx, dy, dz_fixed, nsrc,
                    Lavg.re, Lavg.im,
                    E2avg, B2avg, Savg,
                    LE2avg.re, LE2avg.im,
                    LB2avg.re, LB2avg.im,
                    LSavg.re, LSavg.im,
                    rho_E2, rho_B2, rho_S
                );
            }
        }

        fflush(stdout);
    }

    for (int t = 0; t < L[0]; t++) {
        if (!initialized_time[t]) continue;

        qcd_spinorComponent3d *v_time = &v_by_time[(size_t)t * NV];

        for (int i = 0; i < NV; i++) {
            qcd_destroySpinorComponent3d(&v_time[i]);
        }
    }

    free(v_by_time);

    free(LO_E2);
    free(LO_B2);
    free(LO_S);
    free(sum_E2);
    free(sum_B2);
    free(sum_S);

    qcd_destroyCloverField(&clv);
    qcd_destroyGaugeField(&u);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
