/* Historical HYP2 (piter=4) versus a converged diagnostic (piter=40).
   SU(3) defects of piter=4 are measurements, not acceptance thresholds.
   Geometry and import/export controls follow the validated qX test. */
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
        if (!rank_id) fprintf(stderr, "FAIL: %s\n", message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

/* Real/imaginary component maxima: physical temporal, spatial, absent
   outgoing temporal links at the upper open boundary. */
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
            double dr = fabs(a->D[v][mu][i][j].re-b->D[v][mu][i][j].re);
            double di = fabs(a->D[v][mu][i][j].im-b->D[v][mu][i][j].im);
            if (!isfinite(dr) || !isfinite(di)) bad = 1;
            if (dr > local[k]) local[k] = dr;
            if (di > local[k]) local[k] = di;
        }
    }
    check(bad, "non-finite gauge-field component");
    MPI_Allreduce(local, out, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

typedef struct {
    double max_entry, max_frob, rms_frob, max_det, rms_det;
} su3_report;

static su3_report su3_statistics(qcd_gaugeField *u, const char *label)
{
    qcd_geometry *g = u->geo;
    const int perm[6][3] = {
        {0,1,2}, {1,2,0}, {2,0,1}, {0,2,1}, {1,0,2}, {2,1,0}
    };
    const double threshold[4] = {1e-12,1e-8,1e-4,1e-2};
    double maxima[3]={0,0,0}, sums[2]={0,0};
    double global_max[3], global_sums[2];
    /* Count, four unitarity exceedances, four determinant exceedances. */
    unsigned long long counts[9]={0}, global_counts[9];
    int bad=0;
    for (int t=0; t<g->lL[0]; t++)
    for (int x=0; x<g->lL[1]; x++)
    for (int y=0; y<g->lL[2]; y++)
    for (int z=0; z<g->lL[3]; z++) {
        if (t + g->Pos[0]*g->lL[0] == g->L[0]-1) continue;
        int v=qcd_LEXIC(t,x,y,z,g->lL);
        double frob2=0, entry_max=0;
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            double re=-(i==j), im=0;
            for (int k=0; k<3; k++) {
                double ar=u->D[v][0][k][i].re, ai=u->D[v][0][k][i].im;
                double br=u->D[v][0][k][j].re, bi=u->D[v][0][k][j].im;
                re += ar*br+ai*bi;
                im += ar*bi-ai*br;
            }
            double entry=hypot(re,im);
            if (entry>entry_max) entry_max=entry;
            frob2 += re*re+im*im;
        }
        double det_re=0, det_im=0;
        for (int p=0; p<6; p++) {
            double re=1, im=0;
            for (int row=0; row<3; row++) {
                int col=perm[p][row];
                double ar=u->D[v][0][row][col].re;
                double ai=u->D[v][0][row][col].im;
                double next_re=re*ar-im*ai;
                im=re*ai+im*ar;
                re=next_re;
            }
            double sign=p<3 ? 1.0 : -1.0;
            det_re += sign*re;
            det_im += sign*im;
        }
        double frob=sqrt(frob2), det=hypot(det_re-1,det_im);
        if (!isfinite(frob) || !isfinite(det)) bad=1;
        if (entry_max>maxima[0]) maxima[0]=entry_max;
        if (frob>maxima[1]) maxima[1]=frob;
        if (det>maxima[2]) maxima[2]=det;
        sums[0] += frob2;
        sums[1] += det*det;
        counts[0]++;
        for (int k=0; k<4; k++) {
            if (frob>threshold[k]) counts[1+k]++;
            if (det>threshold[k]) counts[5+k]++;
        }
    }
    check(bad, "non-finite SU(3) diagnostic");
    MPI_Allreduce(maxima,global_max,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(sums,global_sums,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(counts,global_counts,9,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    unsigned long long expected=(unsigned long long)(g->L[0]-1)*g->L[1]*g->L[2]*g->L[3];
    check(global_counts[0]!=expected, "physical temporal-link count");
    double n=(double)global_counts[0];
    su3_report report={global_max[0],global_max[1],sqrt(global_sums[0]/n),
                       global_max[2],sqrt(global_sums[1]/n)};
    check(!isfinite(report.rms_frob) || !isfinite(report.rms_det),
          "non-finite global SU(3) summary");
    if (!rank_id) {
        printf("SU3 %s N=%llu UNIT_MAX_ENTRY=%.16e UNIT_MAX_FROB=%.16e "
               "UNIT_RMS_FROB=%.16e DET_MAX=%.16e DET_RMS=%.16e\n",
               label,global_counts[0],report.max_entry,report.max_frob,
               report.rms_frob,report.max_det,report.rms_det);
        for (int k=0; k<4; k++)
            printf("SU3_EXCEED %s threshold=%.1e UNIT_COUNT=%llu UNIT_FRACTION=%.16e "
                   "DET_COUNT=%llu DET_FRACTION=%.16e\n",label,threshold[k],
                   global_counts[1+k],global_counts[1+k]/n,
                   global_counts[5+k],global_counts[5+k]/n);
    }
    return report;
}

static void print_difference(const char *label, const double d[3])
{
    if (!rank_id)
        printf("%s temporal=%.16e spatial=%.16e boundary=%.16e\n",
               label,d[0],d[1],d[2]);
}

int main(int argc, char **argv)
{
    setvbuf(stdout,NULL,_IONBF,0);
    solver_mpiInit(argc,argv);
    int ranks;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    check(ranks!=384 || argc!=2,"require 384 ranks and a gauge-file path");
    check(solver_getNThread()!=1,"require one compiled thread per MPI rank");
    int L[4]={192,64,64,64}, P[4]={6,4,4,4};
    double theta[4]={0,0,0,0};
    qcd_geometry geo;
    qcd_gaugeField original,output4,output40;
    check(qcd_initGeometry(&geo,L,P,theta,rank_id,ranks),"geometry");
    check(qcd_initEO(&geo),"even/odd geometry");
    check(qcd_initGaugeField(&original,&geo),"original allocation");
    check(qcd_initGaugeField(&output4,&geo),"piter4 allocation");
    check(qcd_initGaugeField(&output40,&geo),"piter40 allocation");
    check(qcd_getGaugeField(argv[1],qcd_GF_OPENQCD,&original),"gauge read");
    solver_parameters par={0};
    par.kappa=0.125;par.bc=0;par.cF=par.cF_prime=1.0;
    check(solver_init(&geo,&par),"solver initialization");
    double d[3];
    (void)su3_statistics(&original,"ORIGINAL");
    check(solver_updateGaugeField(&original),"import gauge");
    check(solver_getOpenQCDGaugeField(&output4),"export gauge");
    differences(&original,&output4,d);print_difference("ROUNDTRIP",d);
    check(fmax(d[0],fmax(d[1],d[2]))>1e-13,"round-trip mismatch");
    check(solver_hypTimeLinks(&output4,0,1,0.5,4),"zero-alpha piter4");
    differences(&original,&output4,d);print_difference("ZERO_ALPHA_P4",d);
    check(fmax(d[0],fmax(d[1],d[2]))>1e-13,"zero-alpha mismatch");
    check(solver_updateGaugeField(&original),"reset before piter4");
    check(solver_hypTimeLinks(&output4,1,1,0.5,4),"piter4 HYP");
    differences(&original,&output4,d);print_difference("HYP4_CHANGE",d);
    check(d[1]>1e-13 || d[2]>1e-13,"piter4 spatial/boundary change");
    check(d[0]<=1e-13,"piter4 no temporal effect");
    (void)su3_statistics(&output4,"PITER4");
    check(solver_updateGaugeField(&original),"reset before repeat4");
    check(solver_hypTimeLinks(&output40,1,1,0.5,4),"repeat4 HYP");
    differences(&output4,&output40,d);print_difference("REPEAT4",d);
    check(fmax(d[0],fmax(d[1],d[2]))>1e-13,"piter4 repeatability");
    check(solver_updateGaugeField(&original),"reset before piter40");
    check(solver_hypTimeLinks(&output40,1,1,0.5,40),"piter40 HYP");
    differences(&original,&output40,d);print_difference("HYP40_CHANGE",d);
    check(d[1]>1e-13 || d[2]>1e-13,"piter40 spatial/boundary change");
    check(d[0]<=1e-13,"piter40 no temporal effect");
    su3_report reference=su3_statistics(&output40,"PITER40");
    check(reference.max_frob>1e-12 || reference.max_det>1e-12,
          "piter40 converged-reference SU(3) tolerance");
    differences(&output4,&output40,d);print_difference("PITER_4_VS_40",d);
    if (!rank_id)
        printf("PASS structural and repeatability checks; piter4 SU(3) defects "
               "reported without an acceptance threshold. No observable comparison yet.\n");
    qcd_destroyGaugeField(&output40);qcd_destroyGaugeField(&output4);
    qcd_destroyGaugeField(&original);qcd_destroyEO(&geo);qcd_destroyGeometry(&geo);
    MPI_Finalize();
    return 0;
}
