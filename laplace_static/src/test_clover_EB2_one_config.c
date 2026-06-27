#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <mpi.h>

#include <qcd.h>

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

    if (
        errno == ERANGE ||
        end == text ||
        *end != '\0' ||
        value < min_value ||
        value > max_value
    ) {
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

static double su3_norm2_traceless_F_from_Q(
    qcd_complex_16 Q[3][3]
) {
    qcd_complex_16 F[3][3];

    /*
      qcdlib convention copied from qcd_makeCloverDiagonal:

        F_munu = anti-Hermitian part of Q_munu with factor 1/8.

      For the density prototype we use the positive norm
        Tr(F^\dagger F) = sum_ab |F_ab|^2

      and subtract the color trace to project approximately to su(3).
      Absolute normalization is not yet the physics point; first we validate
      locality, finiteness, and stable averages.
    */

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

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int myid = 0;
    int numprocs = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if (numprocs != 1) {
        if (myid == 0) {
            fprintf(stderr, "ERROR: this first clover test is intended for exactly 1 MPI rank.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const int nc = read_env_int("NCFG", 1, 1, INT_MAX);

    int L[4] = {48, 24, 24, 24};
    int P[4] = {1, 1, 1, 1};
    double theta[3] = {0.0, 0.0, 0.0};

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

    char gauge_file[1024];

    snprintf(
        gauge_file,
        sizeof(gauge_file),
        "/home/m2130292/Masterarbeit/Em1/cnfg/Em1p4n%d",
        nc
    );

    if (myid == 0) {
        printf("Clover E2/B2 one-config test\n");
        printf("HOST=%s\n", getenv("HOSTNAME") ? getenv("HOSTNAME") : "(unknown)");
        printf("NCFG=%d\n", nc);
        printf("Gauge field: %s\n", gauge_file);
        fflush(stdout);
    }

    if (qcd_getGaugeField(gauge_file, qcd_GF_OPENQCD, &u)) {
        fprintf(stderr, "ERROR: qcd_getGaugeField failed for %s\n", gauge_file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const double plaq = qcd_calculatePlaquette(&u);

    if (myid == 0) {
        printf("Average plaquette = %.16e\n", plaq);
        fflush(stdout);
    }

    if (qcd_calculateCloverField(&clv, &u)) {
        fprintf(stderr, "ERROR: qcd_calculateCloverField failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /*
      qcd_cloverField plane order:
        0: 0/1
        1: 0/2
        2: 0/3
        3: 1/2
        4: 1/3
        5: 2/3

      Here 0 is Euclidean time.
      So E-like planes are 0,1,2 and B-like planes are 3,4,5.
    */

    double sum_E2 = 0.0;
    double sum_B2 = 0.0;
    double sum_S = 0.0;
    double sum_eps = 0.0;

    double per_t_E2[48];
    double per_t_B2[48];
    double per_t_count[48];

    for (int t = 0; t < 48; t++) {
        per_t_E2[t] = 0.0;
        per_t_B2[t] = 0.0;
        per_t_count[t] = 0.0;
    }

    for (long l = 0; l < geo.lV; l++) {
        double E2 = 0.0;
        double B2 = 0.0;

        for (int p = 0; p < 3; p++) {
            E2 += su3_norm2_traceless_F_from_Q(clv.D[l][p]);
        }

        for (int p = 3; p < 6; p++) {
            B2 += su3_norm2_traceless_F_from_Q(clv.D[l][p]);
        }

        const double S = 0.5 * (E2 + B2);
        const double eps = 0.5 * (E2 - B2);

        sum_E2 += E2;
        sum_B2 += B2;
        sum_S += S;
        sum_eps += eps;

        const int t = (int)(l / (24 * 24 * 24));

        if (t >= 0 && t < 48) {
            per_t_E2[t] += E2;
            per_t_B2[t] += B2;
            per_t_count[t] += 1.0;
        }
    }

    if (myid == 0) {
        const double norm = (double)geo.V;

        printf("# Global raw positive clover norms, traceless color projection\n");
        printf("CLOVER_AVG cfg E2_avg B2_avg S_avg eps_avg plaquette\n");
        printf(
            "CLOVER_AVG %d %.16e %.16e %.16e %.16e %.16e\n",
            nc,
            sum_E2 / norm,
            sum_B2 / norm,
            sum_S / norm,
            sum_eps / norm,
            plaq
        );

        printf("# TSLICE cfg t E2_avg B2_avg S_avg eps_avg count\n");

        for (int t = 0; t < 48; t++) {
            const double c = per_t_count[t];

            if (c > 0.0) {
                const double E2t = per_t_E2[t] / c;
                const double B2t = per_t_B2[t] / c;
                const double St = 0.5 * (E2t + B2t);
                const double epst = 0.5 * (E2t - B2t);

                printf(
                    "TSLICE %d %d %.16e %.16e %.16e %.16e %.0f\n",
                    nc,
                    t,
                    E2t,
                    B2t,
                    St,
                    epst,
                    c
                );
            }
        }

        fflush(stdout);
    }

    qcd_destroyCloverField(&clv);
    qcd_destroyGaugeField(&u);
    qcd_destroyGeometry(&geo);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
