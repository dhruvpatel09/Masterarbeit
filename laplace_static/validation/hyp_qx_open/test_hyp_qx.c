#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "qcd.h"
#include "solver.h"

static int rank_id;

static void check(int failed, const char *message)
{
    int any;
    MPI_Allreduce(&failed, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (any) {
        if (rank_id == 0) fprintf(stderr, "FAIL: %s\n", message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

/* Maximum component differences: temporal physical links, spatial links,
   and absent outgoing temporal links at the upper open boundary. */
static void differences(qcd_gaugeField *a, qcd_gaugeField *b, double out[3])
{
    qcd_geometry *g = a->geo;
    double local[3] = {0, 0, 0};
    int bad = 0;
    for (int t=0; t<g->lL[0]; t++)
    for (int x=0; x<g->lL[1]; x++)
    for (int y=0; y<g->lL[2]; y++)
    for (int z=0; z<g->lL[3]; z++) {
        int v = qcd_LEXIC(t,x,y,z,g->lL);
        int gt = t + g->Pos[0]*g->lL[0];
        for (int mu=0; mu<4; mu++)
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            int k = mu ? 1 : (gt == g->L[0]-1 ? 2 : 0);
            double dr = fabs(a->D[v][mu][i][j].re -
                             b->D[v][mu][i][j].re);
            double di = fabs(a->D[v][mu][i][j].im -
                             b->D[v][mu][i][j].im);
            if (!isfinite(dr) || !isfinite(di)) bad = 1;
            if (dr > local[k]) local[k] = dr;
            if (di > local[k]) local[k] = di;
        }
    }
    check(bad, "non-finite gauge-field component");
    MPI_Allreduce(local, out, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

static double temporal_unitarity(qcd_gaugeField *u)
{
    qcd_geometry *g = u->geo;
    double local = 0, global;
    for (int t=0; t<g->lL[0]; t++)
    for (int x=0; x<g->lL[1]; x++)
    for (int y=0; y<g->lL[2]; y++)
    for (int z=0; z<g->lL[3]; z++) {
        if (t + g->Pos[0]*g->lL[0] == g->L[0]-1) continue;
        int v = qcd_LEXIC(t,x,y,z,g->lL);
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            double re = -(i == j), im = 0;
            for (int k=0; k<3; k++) {
                double ar=u->D[v][0][k][i].re;
                double ai=u->D[v][0][k][i].im;
                double br=u->D[v][0][k][j].re;
                double bi=u->D[v][0][k][j].im;
                re += ar*br + ai*bi;
                im += ar*bi - ai*br;
            }
            double d = hypot(re, im);
            if (!isfinite(d)) local = INFINITY;
            else if (d > local) local = d;
        }
    }
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global;
}


/* Maximum |det(U)-1| over physical temporal links.
   Excludes absent outgoing temporal links at the upper open boundary. */
static double temporal_determinant(qcd_gaugeField *u)
{
    qcd_geometry *g = u->geo;
    const int perm[6][3] = {
        {0,1,2}, {1,2,0}, {2,0,1},
        {0,2,1}, {1,0,2}, {2,1,0}
    };
    double local = 0, global;

    for (int t=0; t<g->lL[0]; t++)
    for (int x=0; x<g->lL[1]; x++)
    for (int y=0; y<g->lL[2]; y++)
    for (int z=0; z<g->lL[3]; z++) {
        if (t + g->Pos[0]*g->lL[0] == g->L[0]-1) continue;
        int v = qcd_LEXIC(t,x,y,z,g->lL);
        double det_re = 0, det_im = 0;

        for (int p=0; p<6; p++) {
            double re = 1, im = 0;
            for (int row=0; row<3; row++) {
                int col = perm[p][row];
                double ar = u->D[v][0][row][col].re;
                double ai = u->D[v][0][row][col].im;
                double next_re = re*ar - im*ai;
                im = re*ai + im*ar;
                re = next_re;
            }
            double sign = p < 3 ? 1.0 : -1.0;
            det_re += sign*re;
            det_im += sign*im;
        }

        double defect = hypot(det_re-1.0, det_im);
        if (!isfinite(defect)) local = INFINITY;
        else if (defect > local) local = defect;
    }

    MPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    return global;
}

int main(int argc, char **argv)
{
    setvbuf(stdout, NULL, _IONBF, 0);
    solver_mpiInit(argc, argv);
    int ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    check(ranks != 384 || argc != 2, "require 384 ranks and a gauge-file path");

    int L[4]={192,64,64,64}, P[4]={6,4,4,4};
    double theta[4]={0,0,0,0};
    qcd_geometry geo;
    qcd_gaugeField original, output;

    check(qcd_initGeometry(&geo,L,P,theta,rank_id,ranks), "geometry");
    check(qcd_initEO(&geo), "even/odd geometry");
    check(qcd_initGaugeField(&original,&geo), "original field allocation");
    check(qcd_initGaugeField(&output,&geo), "output field allocation");
    check(qcd_getGaugeField(argv[1],qcd_GF_OPENQCD,&original), "gauge read");

    /* No Dirac solves: all sol[].solver entries remain solver_NONE. */
    solver_parameters par = {0};
    par.kappa = 0.125;
    par.bc = 0;
    par.cF = par.cF_prime = 1.0;
    check(solver_init(&geo,&par), "solver initialization");

    double d[3];
    check(solver_updateGaugeField(&original), "import gauge field");
    check(solver_getOpenQCDGaugeField(&output), "export gauge field");
    differences(&original,&output,d);
    if (!rank_id)
        printf("ROUNDTRIP temporal=%.16e spatial=%.16e boundary=%.16e\n",
               d[0],d[1],d[2]);
    check(fmax(d[0],fmax(d[1],d[2])) > 1e-13, "round-trip mismatch");

    check(solver_hypTimeLinks(&output,0.0,1.0,0.5,10), "zero-alpha HYP");
    differences(&original,&output,d);
    if (!rank_id)
        printf("ZERO_ALPHA temporal=%.16e spatial=%.16e boundary=%.16e\n",
               d[0],d[1],d[2]);
    check(fmax(d[0],fmax(d[1],d[2])) > 1e-13, "zero-alpha mismatch");

    /* Reset before the nonzero-alpha test. */
    check(solver_updateGaugeField(&original), "reset gauge field");
    check(solver_hypTimeLinks(&output,1.0,1.0,0.5,10), "temporal HYP");
    differences(&original,&output,d);
    double before=temporal_unitarity(&original);
    double after=temporal_unitarity(&output);
    if (!rank_id) {
        printf("HYP steps=1 piter=10 alpha=1,1,0.5\n");
        printf("HYP_CHANGE temporal=%.16e spatial=%.16e boundary=%.16e\n",
               d[0],d[1],d[2]);
        printf("TEMPORAL_UNITARITY before=%.16e after=%.16e\n",
               before,after);
    }
    check(d[1] > 1e-13 || d[2] > 1e-13, "spatial/boundary links changed");
    check(d[0] <= 1e-13, "nonzero-alpha HYP had no temporal effect");
    check(!isfinite(after), "non-finite unitarity diagnostic");
    /* output holds the piter=10 result. Keep it for comparison. */
    double det_before = temporal_determinant(&original);
    double det10 = temporal_determinant(&output);

    qcd_gaugeField output20;
    check(qcd_initGaugeField(&output20,&geo), "piter=20 field allocation");

    /* Independent one-step run: reset to the original unsmeared field. */
    check(solver_updateGaugeField(&original), "reset before piter=20");
    check(solver_hypTimeLinks(&output20,1.0,1.0,0.5,20), "piter=20 HYP");

    double change20[3], comparison[3];
    differences(&original,&output20,change20);
    differences(&output,&output20,comparison);

    double unit20 = temporal_unitarity(&output20);
    double det20 = temporal_determinant(&output20);

    if (!rank_id) {
        printf("DETERMINANT original=%.16e piter10=%.16e piter20=%.16e\n",
               det_before,det10,det20);
        printf("UNITARITY piter10=%.16e piter20=%.16e\n",after,unit20);
        printf("HYP20_CHANGE temporal=%.16e spatial=%.16e boundary=%.16e\n",
               change20[0],change20[1],change20[2]);
        printf("PITER_10_VS_20 temporal=%.16e spatial=%.16e boundary=%.16e\n",
               comparison[0],comparison[1],comparison[2]);
        fflush(stdout);
    }

    /* Absolute validation targets, carried over from the tiny test. */
    const double projection_tol = 1e-12;
    const double comparison_tol = 1e-12;

    check(change20[1] > 1e-13 || change20[2] > 1e-13,
          "piter=20 changed spatial/boundary links");
    check(change20[0] <= 1e-13, "piter=20 had no temporal effect");
    check(!isfinite(after) || !isfinite(unit20) ||
          after > projection_tol || unit20 > projection_tol,
          "unitarity tolerance");
    check(!isfinite(det10) || !isfinite(det20) ||
          det10 > projection_tol || det20 > projection_tol,
          "determinant tolerance");
    check(fmax(comparison[0],fmax(comparison[1],comparison[2])) >
          comparison_tol, "piter=10 versus 20 sensitivity");

    if (!rank_id) {
        printf("PASS qX temporal-HYP structural, SU(3), and "
               "piter-sensitivity checks.\n");
        fflush(stdout);
    }

    qcd_destroyGaugeField(&output20);

    qcd_destroyGaugeField(&output);
    qcd_destroyGaugeField(&original);
    qcd_destroyEO(&geo);
    qcd_destroyGeometry(&geo);
    MPI_Finalize();
    return 0;
}
