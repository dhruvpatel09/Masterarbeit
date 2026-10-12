/**
 * @file solver.c
 * @brief Interface to openQCD's solvers.
 * @details The functions in this file provide an interface to the openQCD solvers, 
 *          such that they can be used from programs based on the qcdlib library.
 *          Version based on openQCD-2.4.2 and Bjoern Leder's Wilson loop program
 *          For usage instructions see the file README.solver
 * @version 2.4.2
 * @author Tomasz Korzec
 * @date 2015-2026
 * @copyright GNU Public License.
 *
 */


#define MAIN_PROGRAM  ///< forces openQCD's global variables to be defined (not just declared)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "mpi.h"
#include "su3.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "archive.h"
#include "sflds.h"
#include "linalg.h"
#include "dirac.h"
#include "sap.h"
#include "dfl.h"
#include "version.h"
#include "uflds.h"
#include "sw_term.h"
#include "forces.h"
#include "update.h"
#include "tcharge.h"
#include "wflow.h"
#include "random.h"
#include "global.h"
#include "qcd.h"
#include "lnk.h"
#include "smear.h"
#include "solver.h"


/** @def N0
 *  @brief global lattice size in 0 direction
 */
#define N0 (NPROC0*L0)
/** @def N1
 *  @brief global lattice size in 1 direction
 */
#define N1 (NPROC1*L1)
/** @def N2
 *  @brief global lattice size in 2 direction
 */
#define N2 (NPROC2*L2)
/** @def N3
 *  @brief global lattice size in 3 direction
 */
#define N3 (NPROC3*L3)


/* private global variables */
static int my_rank;
static int *status=NULL;  // status of last solve or dfl generation
static int initialized=0; // set to 1 after initialization
static int hasU=0;        // set to 1 after gaugeField update
static int hasDFL=0;      // set to 1 after dfl generation
static int flowinit=0;    // set to 1 after first flow call

static bc_parms_t bcp;    // boundary conditions
static lat_parms_t lat;   // other lattice parameters


/* private helper functions */

static void maxn(int *n,int m)
{
   if ((*n)<m)
      (*n)=m;
}//end maxn


static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpr;

   dp=dfl_parms();
   dpr=dfl_pro_parms();

   maxn(nws,dp.Ns+2);
   maxn(nwv,2*dpr.nmx_gcr+3);
   maxn(nwvd,2*dpr.nkv+4);
}//end dfl_wsize

static void solver_wsize(int isp,int nsds,int *nws,int *nwv,int *nwvd)
{
   solver_parms_t sp;

   sp=solver_parms(isp);

   if (sp.solver==CGNE)
      maxn(nws,nsds+11);
   else if (sp.solver==SAP_GCR)
      maxn(nws,nsds+2*sp.nkv+5);
   else if (sp.solver==DFL_SAP_GCR)
   {
      maxn(nws,nsds+2*sp.nkv+6);
      dfl_wsize(nws,nwv,nwvd);
   }
   else if (sp.solver==SOLVERS)
   {
      // no space needed
   }
   else
   {
      error_root(1,1,"solver_wsize [solver.c]","Unknown or unsupported solver");  
   }
}//end solver_wsize


static void wsize(int *nws,int *nwv,int *nwvd)
{
   (*nws)=0;
   (*nwv)=0;
   (*nwvd)=0;
   int i;

   for(i=0; i<solver_ISPMAX; i++)
   {
      solver_wsize(i,4,nws,nwv,nwvd);
   }
}
//end wsize


static void qcd2oqcd(qcd_vector *v, spinor_dble *vo)
{
   int ix, iy;
   int x0,x1,x2,x3;
   qcd_geometry *geo=v->geo;
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      ix = qcd_LEXIC(x0,x1,x2,x3,geo->lL);      //qcdlib index
      iy = ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];  //openQCD index

      /* the dirac/color/complex components have the same memory layout */
      memcpy(vo+iy,&(v->D[ix][0][0].re),192);  // 192 = 4x3xcomplex
   }
//   printf("process %i: vo[3 2 1 1][2][1] = %f %+fi\n"
//          "             v[3 2 1 1][2][1] = %f %+fi\n",
//          my_rank,vo[ipt[1+L3*1+L2*L3*2+L1*L2*L3*3]].c3.c2.re, vo[ipt[1+L3*1+L2*L3*2+L1*L2*L3*3]].c3.c2.im,
//                  v->D[qcd_LEXIC(3,2,1,1,geo->lL)][2][1].re, v->D[qcd_LEXIC(3,2,1,1,geo->lL)][2][1].im);
}//end qcd2oqcd

static void oqcd2qcd(spinor_dble *vo, qcd_vector *v)
{
   int ix, iy;
   int x0,x1,x2,x3;
   qcd_geometry *geo=v->geo;
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      ix = qcd_LEXIC(x0,x1,x2,x3,geo->lL);      //qcdlib index
      iy = ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0];  //openQCD index

      /* the dirac/color/complex components have the same memory layout */
      memcpy(&(v->D[ix][0][0].re),vo+iy,192);  // 192 = 4x3xcomplex
   }
}//end oqcd2qcd


static double dbl(qflt q)
{
   double q0 = q.q[0];
   return q0;
}


/* externally visible functions */
#if (defined _OPENMP)

void solver_mpiInit(int argc,char **argv)
{

   int provided;

   MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);

   qcd_errorRoot(provided<MPI_THREAD_FUNNELED,1,"mpi_init_thread [testsolv_open-2.4.2.c]",
              "Thread-safety level MPI_THREAD_FUNNELED is not supported");

   omp_set_num_threads(solver_getNThread());
   //printf("nthread = %i\n",solver_getNThread());
   omp_set_dynamic(0);
   omp_set_nested(0);
}

int solver_getNThread(void)
{
   return NTHREAD;
}

#else

void mpi_init(int argc,char **argv)
{
   MPI_Init(&argc,&argv);
}

int solver_getNThread(void)
{
   return 1;
}
#endif





int solver_init(qcd_geometry *geo, solver_parameters *par)
{
   int i,err=0;
   su3 *us=NULL;         //single precision gauge field
   su3_dble *ud=NULL;    //double precision gauge field
   solver_t st=0;
   int needsSAP=0;
   int needsDFL=0;
   int nws,nwv,nwvd; // number of different types of spinors on workspace

   error(initialized,1,"solver_init [solver.c]",
                       "solver_init can be called only once");

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   /* check whether global lattice sizes match */
   if( (geo->L[0] != N0) || (geo->L[1] != N1) ||
       (geo->L[2] != N2) || (geo->L[3] != N3))
   {
      err=1;
   }
   error(err,1,"solver_init [solver.c]",
               "Global lattice size  missmatch between qcdlib and openQCD");

   /* check whether local lattice sizes match */
   if((geo->lL[0]!= L0) || (geo->lL[1]!= L1) ||
      (geo->lL[2]!= L2) || (geo->lL[3]!= L3))
   {
      err=1;
   }
   error(err,1,"solver_init [solver.c]",
               "Local lattice size  missmatch between qcdlib and openQCD");

   status = alloc_std_status();
   
   /* set up lattice parameters */
   lat=set_lat_parms(0.0, 1.0, 1, &(par->kappa), par->isw, par->csw);
   /* set up boundary conditions */
   bcp=set_bc_parms(par->bc, 1.0, 1.0, par->cF,par->cF_prime, par->phi, par->phi_prime, par->theta);
   /* create the openQCD neighborhood relations */
   check_machine();
   geometry();

   /* check whether the map "my_rank <--> cartesian grid" matches */
   if((cpr[0] != geo->Pos[0]) || (cpr[1] != geo->Pos[1]) ||
      (cpr[2] != geo->Pos[2]) || (cpr[3] != geo->Pos[3]))
   {
      err=1;
   }
   error(err,1,"solver_init [solver.c]",
                 "process grid missmatch between qcdlib and openQCD\n"
                 "(Try turning process blocks off)");


   /* get space for gauge field (single + double) */
   ud = udfld();
   us = ufld();
   error((ud==0) || (us==0), 1, "solver_init [solver.c]",
                         "could not allocate openQCD gauge fields");

   set_sw_parms(0.5/par->kappa-4.0);
   set_tm_parms(0);    //TM term lives on all sites

   /* loop through solvers and initialize each. */
   for(i=0; i<solver_ISPMAX; i++)
   {
      switch (par->sol[i].solver)
      {
         case solver_NONE: st=SOLVERS;
                           break;
         case solver_CGNE: st=CGNE;
                           break;
         case solver_MSCG: st=MSCG;
                           break;
         case solver_SAP_GCR: st=SAP_GCR;
                              needsSAP = 1;
                              break;
         case solver_DFL_SAP_GCR: st=DFL_SAP_GCR;
                                  needsSAP = 1;
                                  needsDFL = 1;
                                  break;
         default: error_root(err,1,"solver_init [solver.c]",
                                   "unsupported solver\n");
      }
      if( st != SOLVERS )
         set_solver_parms(i, st, par->sol[i].nkv, par->sol[i].isolv,
                          par->sol[i].nmr, par->sol[i].ncy, par->sol[i].nmx,
                          par->sol[i].istop, par->sol[i].res);
   }
   /* if necessary: set sap parameters */
   if( needsSAP )
   {
      set_sap_parms(par->sap.bs,1,4,5);
   }

   /* if necessary: set dfl parameters */
   if( needsDFL )
   {
      start_ranlux(0,137); // not so important, used for deflation subspace only
      set_dfl_parms(par->dfl.bs,par->dfl.Ns);
      set_dfl_gen_parms(par->dfl.kappa, par->dfl.mu, par->dfl.ninv, par->dfl.nmr, par->dfl.ncy);
      set_dfl_pro_parms(par->dfl.nkv,par->dfl.nmx,par->dfl.res, 
                        par->dfl.nmx_gcr, par->dfl.res_gcr);
   }

   /* reserve space for single and double precision spinors */
   wsize(&nws,&nwv,&nwvd);
   alloc_ws(nws);
   wsd_uses_ws();
   alloc_wv(nwv);
   alloc_wvd(nwvd);
   
   initialized = 1;
   return(0);
}//end solver_init


int solver_hypLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter)
{
   int ret=0;
   error(!initialized,1,"solver_hypLinks [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_hypLinks [solver.c]",
                 "solver_updateGaugeField must be called first");
   
   unset_ud_phase();
   hyp_links(alpha1, alpha2, alpha3, 1, piter, 0);
   set_ud_phase();
   
   ret = solver_getOpenQCDGaugeField(qcdU);
   return ret;
}// end solver_hypLinks

int solver_hypSpaceLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter)
{
   int ret=0;
   error(!initialized,1,"solver_hypSpaceLinks [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_hypSpaceLinks [solver.c]",
                 "solver_updateGaugeField must be called first");
   
   unset_ud_phase();
   hyp_space_links(alpha1, alpha2, alpha3, 1, piter, 0);
   set_ud_phase();
   
   ret = solver_getOpenQCDGaugeField(qcdU);
   return ret;
}// end solver_hypSpaceLinks


int solver_hypTimeLinks(qcd_gaugeField *qcdU, double alpha1, double alpha2, double alpha3, int piter)
{
   int ret=0;
   error(!initialized,1,"solver_hypTimeLinks [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_hypTimeLinks [solver.c]",
                 "solver_updateGaugeField must be called first");
   
   unset_ud_phase();
   hyp_time_links(alpha1, alpha2, alpha3, 1, piter, 0);
   set_ud_phase();
   
   ret = solver_getOpenQCDGaugeField(qcdU);
   return ret;
}// end solver_hypTimeLinks


int solver_hypLinks3d(qcd_gaugeField *qcdU, double alpha1, double alpha2, int piter)
{
   int ret=0;
   error(!initialized,1,"solver_hypLinks3d [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_hypLinks3d [solver.c]",
                 "solver_updateGaugeField must be called first");
   
   unset_ud_phase();
   hyp_spatial_links(udfld(), alpha1, alpha2, 1, piter, 0);
   set_ud_phase();
   
   ret = solver_getOpenQCDGaugeField(qcdU);
   return ret;
}// end solver_hypLinks3d


void solver_wilsonFlow(int dnms, double eps)
{
   error(!initialized,1,"solver_wilsonFlow [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_wilsonFlow [solver.c]",
                 "solver_updateGaugeField must be called first");
   if (flowinit==0)
   {
      alloc_wfd(1);
      flowinit=1;
   }
   fwd_rk3(dnms, eps);
   return;
}//end solver_wilsonFlow

double solver_computeWsl(double *wsl)
{
   double wslsum;
   error(!initialized,1,"solver_computeWsl [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_computeWsl [solver.c]",
                 "solver_updateGaugeField must be called first");
   wslsum = plaq_action_slices(wsl);
   return(wslsum);
}

double solver_computeYsl(double *ysl)
{
   double yslsum;
   error(!initialized,1,"solver_computeYsl [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_computeYsl [solver.c]",
                 "solver_updateGaugeField must be called first");
   yslsum = ym_action_slices(ysl);
   return(yslsum);
}

double solver_computeQsl(double *qsl)
{
   double qslsum;
   error(!initialized,1,"solver_computeQsl [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_computeQsl [solver.c]",
                 "solver_updateGaugeField must be called first");
   qslsum = tcharge_slices(qsl);
   return(qslsum);
}

double **solver_allocFlowObs(int nn)
{
   double **pp;
   double *p;
   int i;
   error(!initialized,1,"solver_allocFlowObs [solver.c]",
                        "solver_init must be called first");
   pp=amalloc((nn+1)*sizeof(*pp),3);
   p=amalloc((nn+1)*(N0)*sizeof(*p),4);

   error((pp==NULL)||(p==NULL),1,"solver_allocFlowObs [solver.c]",
         "Unable to allocate data arrays");

   for (i=0;i<nn+1;i++)
   {
      pp[i] = p;
      p+=N0;
   }
   return(pp);
}

void solver_freeFlowObs(double **obs)
{
    afree(obs[0]);
    afree(obs);
}

int solver_getOpenQCDGaugeField(qcd_gaugeField *u)
{
   su3_dble *ud;
   int ix,iy,iu;
   int x0,x1,x2,x3,mu;
   int ip[4];
   qcd_geometry *geo;

   error(!initialized,1,"solver_getOpenQCDGaugeField [solver.c]",
                        "solver_init must be called first");

   error(!(u->initialized),1,"solver_getOpenQCDGaugeField [solver.c]",
                             "Gauge field was not initialized");
   geo = u->geo;
   ud = udfld();
   error(ud==NULL,1,"solver_updateGaugeField [solver.c]",
                    "openQCD gauge field was not initialized");
   unset_ud_phase();

   copy_bnd_ud();
   
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      ix = qcd_LEXIC(x0,x1,x2,x3,geo->lL);     // qcdlib index
      iy = ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0]; // openQCD index
      if((x0+x1+x2+x3)%2) //odd points
      {
         iu = 8*(iy-VOLUME/2);
         for(mu=0; mu<4; mu++)
            memcpy(&(u->D[ix][mu][0][0].re),ud+(iu+2*mu), 144); // 144 = 3x3xcomplex
      }else //even points
      {
         for(mu=0; mu<4; mu++)
         {
            iu = iup[iy][mu];
            if(iu<VOLUME) //boundary points are treated later
            {
               iu = 8*(iu-VOLUME/2)+2*mu+1;
               memcpy(&(u->D[ix][mu][0][0].re),ud+iu, 144);
            }else
            {
                //figure out index in the outer boundary
                //TODO: try to do this more efficiently
                if(mu==0)
                {
                   plaq_uidx(0,iy,ip);
                   iu = ip[0];
                }
                else
                {
                   plaq_uidx(mu-1,iy,ip);
                   iu = ip[2];
                }
                memcpy(&(u->D[ix][mu][0][0].re),ud+iu, 144);
            }
         }
      }
   }
   u->BPup2date = 0;
   u->BMup2date = 0;
   set_ud_phase();
   return(0);

}//end solver_getOpenQCDGaugeField

int solver_updateGaugeField(qcd_gaugeField *u)
{
   su3_dble *ud;
   int ix,ix2,iy,iu,ib;
   int x0,x1,x2,x3,mu;
   int x[4];
   int ibc;
   qcd_geometry *geo;

   error(!initialized,1,"solver_updateGaugeField [solver.c]",
                        "solver_init must be called first");

   error(!(u->initialized),1,"solver_updateGaugeField [solver.c]",
                             "Gauge field was not initialized");
   geo = u->geo;
   ud = udfld();
   error(ud==NULL,1,"solver_updateGaugeField [solver.c]",
                    "openQCD gauge field was not initialized");
   unset_ud_phase();

   //TODO copy only the necessary mu directions of u
   error(qcd_communicateGaugeP(u),1,"solver_updateGaugeField [solver.c]", 
                                    "Communication error");
   //copy available links while communication goes on
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      ix = qcd_LEXIC(x0,x1,x2,x3,geo->lL);     // qcdlib index
      iy = ipt[x3+L3*x2+L2*L3*x1+L1*L2*L3*x0]; // openQCD index
      if((x0+x1+x2+x3)%2) //odd points
      {
         iu = 8*(iy-VOLUME/2);
         for(mu=0; mu<4; mu++)
            memcpy(ud+(iu+2*mu), &(u->D[ix][mu][0][0].re), 144); // 144 = 3x3xcomplex
      }else //even points
      {
         for(mu=0; mu<4; mu++)
         {
            iu = iup[iy][mu];
            if(iu<VOLUME) //boundary points are treated later
            {
               iu = 8*(iu-VOLUME/2)+2*mu+1;
               memcpy(ud+iu,  &(u->D[ix][mu][0][0].re), 144);
            }
         }
      }
   }
   qcd_waitall(geo);
   //some odd inner-boundary points are still incomplete
   for(ib=0; ib<geo->edgePoints; ib++)   //TODO (maybe) loop only over odd points
   {
      ix = geo->edge[ib]; // qcdlib index of an inner boundary point
      qcd_antilexic(x,ix,geo->lL);
      if( (x[0]+x[1]+x[2]+x[3])%2 )
      {
         for(mu=0; mu<4; mu++)
         {
            ix2 = geo->minus[ix][mu];
            if( ix2>geo->lV ) //ix2 is an even point on the newly copied boundary
            {
               iy = ipt[x[3]+L3*x[2]+L2*L3*x[1]+L1*L2*L3*x[0]]; //openQCD index of an odd inner boundary point
               memcpy(ud+(8*(iy-VOLUME/2)+2*mu+1), &(u->D[ix2][mu][0][0].re), 144);
            }
         }
      }
   }

   set_flags(UPDATED_UD);

   ibc=check_bc(64.0*DBL_EPSILON);
   error_root(ibc!=1,1,"solver_updateGaugeField [solver.c]",
              "Incompatible boundary conditions");

   //for testing:
   //export_cnfg("/media/data/testfld");

   set_bc();
   set_ud_phase();
   hasU=1;
   return(0);
}//end solver_updateGaugeField


int solver_applyDw(qcd_vector *Dv, qcd_vector *v, double muTM)
{
   spinor_dble **wsd, *vo, *Dvo; //openQCD versions of v and Dv

   error(!initialized,1,"solver_applyDw [solver.c]",
                        "solver_init must be called first");

   error(!hasU,1,"solver_applyDw [solver.c]",
                 "solver_updateGaugeField must be called first");

   error(!(v->initialized) || !(Dv->initialized),1,
                             "solver_applyDw [solver.c]",
                             "qcdlib spinors not initialized");

   wsd=reserve_wsd(2);
   vo =wsd[0];
   Dvo=wsd[1];
   error((wsd==NULL) || (vo==NULL) || (Dvo==NULL),1,
                             "solver_applyDw [solver.c]",
                             "openQCD spinors not initialized");

   /* translate v to openQCD format */
   qcd2oqcd(v,vo);

   /* apply the Wilson (clover) Dirac opertor with twisted mass muTM */
   sw_term(NO_PTS);      // compute sw term, do not invert
   Dw_dble(muTM,vo,Dvo); 

   /* translate Dvo to libqcd format */
   oqcd2qcd(Dvo,Dv);

   release_wsd();
   return 0;
}//end solver_applyDw


int solver_solveDw(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp)
{
   spinor_dble **wsd, *vo, *Dinvvo; //openQCD versions of v and Dv
   solver_parms_t sp;
   sap_parms_t sap;
   int ifail0[2];
   int *stat0;

   error(!initialized,1,"solver_solveDw [solver.c]",
                        "solver_init must be called first");
   error(!hasU,1,"solver_solveDw [solver.c]",
                 "solver_updateGaugeField must be called first");
   error(!(v->initialized) || !(Dinvv->initialized),1,
                             "solver_solveDw [solver.c]",
                             "qcdlib spinors not initialized");

   stat0=alloc_std_status();
   ifail0[0]=0;
   ifail0[1]=0;
   reset_std_status(stat0);
   reset_std_status(status);
   
   wsd=reserve_wsd(2);
   vo =wsd[0];
   Dinvvo=wsd[1];
   error((vo==NULL) || (Dinvvo==NULL),1,
                             "solver_solveDw [solver.c]",
                             "openQCD spinors not initialized");

   /* translate v to openQCD format */
   qcd2oqcd(v,vo);

   /* solve the Dirac equation using solver isp */
   sp=solver_parms(isp);
   //set_sw_parms(0.5/kappa-4.0);
   if (sp.solver==CGNE)
   {
      mulg5_dble(VOLUME_TRD,2,vo);
      tmcg(sp.nmx,sp.istop,sp.res,muTM,vo,vo,ifail0,stat0);
      acc_std_status("tmcg",ifail0,stat0,0,status);

      if (ifail0[0]<0)
      {
         print_status("tmcg",ifail0,stat0);
         error_root(1,1,"solver_solveDW [solver.c]","CGNE solver failed "
                    "(mu = %.2e, parameter set no 0)",muTM);
      }
      //TODO: better mechanism for ^

      Dw_dble(-muTM,vo,Dinvvo);
      mulg5_dble(VOLUME_TRD,2,Dinvvo);
   }
   else if (sp.solver==SAP_GCR)
   {
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      sap_gcr(sp.nkv,sp.nmx,sp.istop,sp.res,muTM,vo,Dinvvo,ifail0,stat0);
      acc_std_status("sap_gcr",ifail0,stat0,0,status);

      if (ifail0[0]<0)
      {
         print_status("sap_gcr",ifail0,stat0);
         error_root(1,1,"solver_solveDW [solver.c]","SAP_GCR solver failed "
                    "(mu = %.2e, parameter set no 0)",muTM);
      }
      //TODO: better mechanism for ^
   }
   else if (sp.solver==DFL_SAP_GCR)
   {
      error(!hasDFL,1,"solver_solveDw [solver.c]",
                      "solver_dflGeneration must be called first");
      sap=sap_parms();
      set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

      dfl_sap_gcr2(sp.nkv,sp.nmx,sp.istop,sp.res,muTM,vo,Dinvvo,ifail0,stat0);
      acc_std_status("dfl_sap_gcr2",ifail0,stat0,0,status);

      if ((ifail0[0]<-2)||(ifail0[1]<0))
      {
         print_status("dfl_sap_gcr2",ifail0,stat0);
         error_root(1,1,"solver_solveDW [solver.c]","DFL_SAP_GCR solver failed "
                    "(mu = %.2e, parameter set no 0)",muTM);
      }

   }
   else
      error_root(1,1,"solver_solveDw [solver.c]",
                 "Unknown or unsupported solver");


   /* translate Dvo to libqcd format */
   oqcd2qcd(Dinvvo,Dinvv);

   release_wsd();
   free(stat0);
   return 0;
}//end solver_solveDw


int solver_solveDwGuessed(qcd_vector *Dinvv, qcd_vector *v, double muTM, int isp)
{
   spinor_dble **wsd, *vo, *Dinvvo; //openQCD versions of v and Dv
   spinor_dble *rho; //residuum
   solver_parms_t sp;
   sap_parms_t sap;
   double res1;
   int ifail0[2];
   int *stat0;

   error(!initialized,1,"solver_solveDwFrom [solver.c]",
                        "solver_init must be called first");
   error(!hasU,1,"solver_solveDwFrom [solver.c]",
                 "solver_updateGaugeField must be called first");
   error(!(v->initialized) || !(Dinvv->initialized),1,
                             "solver_solveDwFrom [solver.c]",
                             "qcdlib spinors not initialized");

   stat0=alloc_std_status();
   ifail0[0]=0;
   ifail0[1]=0;
   reset_std_status(stat0);
   reset_std_status(status);

   wsd=reserve_wsd(3);
   vo =wsd[0];
   Dinvvo=wsd[1];
   rho=wsd[2];
   error((vo==NULL) || (Dinvvo==NULL) || (rho==NULL),1,
                             "solver_solveDwFrom [solver.c]",
                             "openQCD spinors not initialized");

   /* translate v to openQCD format */
   qcd2oqcd(v,vo);
   
   /* translate initial guess Dinvvo to openQCD format */
   qcd2oqcd(Dinvv,Dinvvo);
   
   /* compute the relative residuum of initial guess */
   sw_term(NO_PTS);      // compute sw term, do not invert
   Dw_dble(muTM,Dinvvo,rho);
   mulr_spinor_add_dble(VOLUME_TRD,2,rho,vo,-1.0); // residuum of guess
   res1=sqrt(dbl(norm_square_dble(VOLUME_TRD,2,rho)) / dbl(norm_square_dble(VOLUME_TRD,2,vo)));

   if (res1<1.0)
   {
      /* solve the D Dinvo = rho  using solver isp */
      sp=solver_parms(isp);
      if (sp.solver==CGNE)
      {
         mulg5_dble(VOLUME_TRD,2,rho);
         tmcg(sp.nmx,sp.istop,sp.res/res1,muTM,rho,rho,ifail0,stat0);
         acc_std_status("tmcg",ifail0,stat0,0,status);
         if (ifail0[0]<0)
         {
            print_status("tmcg",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","CGNE solver failed "
                    "(mu = %.2e, parameter set no 0)",muTM);
         }
         //TODO: better mechanism for ^

         Dw_dble(-muTM,rho,vo);
         mulg5_dble(VOLUME_TRD,2,vo);
      }
      else if (sp.solver==SAP_GCR)
      {
         sap=sap_parms();
         set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

         sap_gcr(sp.nkv,sp.nmx,sp.istop,sp.res/res1,muTM,rho,vo,ifail0,stat0);
         acc_std_status("sap_gcr",ifail0,stat0,0,status);
         if (ifail0[0]<0)
         {
            print_status("sap_gcr",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","SAP_GCR solver failed "
                       "(mu = %.2e, parameter set no 0)",muTM);
         }
         //TODO: better mechanism for ^
      }
      else if (sp.solver==DFL_SAP_GCR)
      {
         error(!hasDFL,1,"solver_solveDwGuessed [solver.c]",
                         "solver_dflGeneration must be called first");
         sap=sap_parms();
         set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.istop,sp.res/res1,muTM,rho,vo,ifail0,stat0);
         acc_std_status("dfl_sap_gcr2",ifail0,stat0,0,status);
         if ((ifail0[0]<-2)||(ifail0[1]<0))
         {
            print_status("dfl_sap_gcr2",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","DFL_SAP_GCR solver failed "
                       "(mu = %.2e, parameter set no 0)",muTM);
         }      
      }
      else
         error_root(1,1,"solver_solveDwGuessed [solver.c]",
                    "Unknown or unsupported solver");
      /* update the solution */
      mulr_spinor_add_dble(VOLUME_TRD,2,Dinvvo,vo,-1.0);
   }
   else
   {
      /* first guess is useless, make standard solve */
      sp=solver_parms(isp);
      if (sp.solver==CGNE)
      {
         mulg5_dble(VOLUME_TRD,2,vo);
         tmcg(sp.nmx,sp.istop,sp.res,muTM,vo,vo,ifail0,stat0);
         acc_std_status("tmcg",ifail0,stat0,0,status);
         if (ifail0[0]<0)
         {
            print_status("tmcg",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","CGNE solver failed "
                    "(mu = %.2e, parameter set no 0)",muTM);
         }
         //TODO: better mechanism for ^

         Dw_dble(-muTM,vo,Dinvvo);
         mulg5_dble(VOLUME_TRD,2,Dinvvo);
      }
      else if (sp.solver==SAP_GCR)
      {
         sap=sap_parms();
         set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

         sap_gcr(sp.nkv,sp.nmx,sp.istop,sp.res,muTM,vo,Dinvvo,ifail0,stat0);
         acc_std_status("sap_gcr",ifail0,stat0,0,status);
         if (ifail0[0]<0)
         {
            print_status("sap_gcr",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","SAP_GCR solver failed "
                       "(mu = %.2e, parameter set no 0)",muTM);
         }
         //TODO: better mechanism for ^
      }
      else if (sp.solver==DFL_SAP_GCR)
      {
         error(!hasDFL,1,"solver_solveDwGuessed [solver.c]",
                         "solver_dflGeneration must be called first");
         sap=sap_parms();
         set_sap_parms(sap.bs,sp.isolv,sp.nmr,sp.ncy);

         dfl_sap_gcr2(sp.nkv,sp.nmx,sp.istop,sp.res,muTM,vo,Dinvvo,ifail0,stat0);
         acc_std_status("dfl_sap_gcr2",ifail0,stat0,0,status);
         if ((ifail0[0]<-2)||(ifail0[1]<0))
         {
            print_status("dfl_sap_gcr2",ifail0,stat0);
            error_root(1,1,"solver_solveDWGuessed [solver.c]","DFL_SAP_GCR solver failed "
                       "(mu = %.2e, parameter set no 0)",muTM);
         }     
      }
      else
         error_root(1,1,"solver_solveDwGuessed [solver.c]",
                    "Unknown or unsupported solver");
   }

   /* translate Dvo to libqcd format */
   oqcd2qcd(Dinvvo,Dinvv);

   release_wsd();
   free(stat0);
   return 0;
}//end solver_solveDwGuessed


int solver_dflGeneration(void)
{
   dfl_parms_t dfl;
   int ifail0[2];

   error(!initialized,1,"solver_dflGeneration [solver.c]",
                        "solver_init must be called first");
   error(!hasU,1,"solver_dflGeneration [solver.c]",
                        "solver_updateGaugeField must be called first");

   reset_std_status(status);
   dfl=dfl_parms();
   if(dfl.Ns)
   {
      dfl_modes2(ifail0,status);
      /*
      if(my_rank==0)
      {
         printf("Generation of the deflation subspace\n");
         print_status("dfl_modes2",ifail0,status);
      }*/
      error_root((ifail0[0]<-2)||(ifail0[1]<0),1,"solver_dflGeneration [solver.c]",
                    "Deflation subspace generation failed");
   }
   hasDFL=1;
   return 0;
}//end solver_dflGeneration


int solver_setKappa(double kappa)
{
   error(!initialized,1,"solver_setKappa [solver.c]",
                        "solver_init must be called first");

   set_sw_parms(0.5/kappa-4.0);
   return 0;
}//end solver_setKappa



int solver_readDFL(solver_parameters *spar)
{
   int bs[4],Ns;
   int ninv,nmr,ncy,nkv,nmx;
   double kappa,mudfl,res;
   int my_rank;
   int err=0;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if (my_rank==0)
   {
      if (qcd_findOpenQCDSection("Deflation subspace")==-1)
         err=1;
      if (qcd_readOpenQCDLine("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3)==-1)
         err=1;
      if (qcd_readOpenQCDLine("Ns","%d",&Ns)==-1)
         err=1;
      
      if (qcd_findOpenQCDSection("Deflation subspace generation")==-1)
         err=1;
      if (qcd_readOpenQCDLine("kappa","%lf",&kappa)==-1)
         err=1;
      if (qcd_readOpenQCDLine("mu","%lf",&mudfl)==-1)
         err=1;
      if (qcd_readOpenQCDLine("ninv","%d",&ninv)==-1)
         err=1;
      if (qcd_readOpenQCDLine("nmr","%d",&nmr)==-1)
         err=1;
      if (qcd_readOpenQCDLine("ncy","%d",&ncy)==-1)
         err=1;
      
      if (qcd_findOpenQCDSection("Deflation projection")==-1)
         err=1;
      if (qcd_readOpenQCDLine("nkv","%d",&nkv)==-1)
         err=1;
      if (qcd_readOpenQCDLine("nmx","%d",&nmx)==-1)
         err=1;
      if (qcd_readOpenQCDLine("res","%lf",&res)==-1)
         err=1;
   }

   MPI_Bcast(&err,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Ns,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&mudfl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&ninv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if (err==1)
      return 1;
   else
   {
      spar->dfl.bs[0] = bs[0];
      spar->dfl.bs[1] = bs[1];
      spar->dfl.bs[2] = bs[2];
      spar->dfl.bs[3] = bs[3];
      spar->dfl.Ns    = Ns;
      spar->dfl.nkv   = nkv;
      spar->dfl.nmx   = nmx;
      spar->dfl.res   = res;
      spar->dfl.ninv  = ninv;
      spar->dfl.nmr   = nmr;
      spar->dfl.ncy   = ncy;
      spar->dfl.kappa = kappa;
      spar->dfl.mu    = mudfl;
      return 0;
   }
}



int solver_readSAP(solver_parameters *spar)
{
   int bs[4];
   int my_rank;
   int err = 0;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   if (my_rank==0)
   {
      if( qcd_findOpenQCDSection("SAP")==-1 )
         err = 1;
      if( qcd_readOpenQCDLine("bs","%d %d %d %d",bs,bs+1,bs+2,bs+3)==-1)
         err = 1;
   }

   MPI_Bcast(&err,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(bs,4,MPI_INT,0,MPI_COMM_WORLD);
   if(err==1)
      return 1;
   else
   {
      spar->sap.bs[0] = bs[0];
      spar->sap.bs[1] = bs[1];
      spar->sap.bs[2] = bs[2];
      spar->sap.bs[3] = bs[3];
      return 0;
   }
}



int solver_readSolver(solver_parameters *spar, int isp)
{
   int my_rank,ids;
   int nkv,isolv,nmr,ncy,nmx;
   double res;
   char line[NAME_SIZE];
   int err=0;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   ids=solver_NONE;
   nkv=0;
   isolv=0;
   nmr=0;
   ncy=0;
   nmx=0;
   int istop;
   res=1.0;
   
   if (my_rank==0)
   {
      if (isp<0 || isp>31)
         err=1;

      sprintf(line,"Solver %d",isp);
      if(qcd_findOpenQCDSection(line)==-1)
         err = 1;
      
      if(qcd_readOpenQCDLine("solver","%s",line)==-1)
         err = 1;

      if (strcmp(line,"MSCG")==0)
         ids=solver_MSCG;
      else if (strcmp(line,"SAP_GCR")==0)
      {
         ids=solver_SAP_GCR;
         if (qcd_readOpenQCDLine("nkv","%d",&nkv)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("isolv","%d",&isolv)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("nmr","%d",&nmr)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("ncy","%d",&ncy)==-1)
            err = 1;
      }
      else if (strcmp(line,"DFL_SAP_GCR")==0)
      {
         ids=solver_DFL_SAP_GCR;
         if (qcd_readOpenQCDLine("nkv","%d",&nkv)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("isolv","%d",&isolv)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("nmr","%d",&nmr)==-1)
            err = 1;
         if (qcd_readOpenQCDLine("ncy","%d",&ncy)==-1)
            err = 1;
      }
      else if (strcmp(line,"CGNE")==0)
      {
         ids=solver_CGNE;  
      }
      else
      {
         fprintf(stderr,"Error in solver_readParameters! Unknown solver %s\n"
                 ,line);
         err = 1;
      }

      if (qcd_readOpenQCDLine("nmx","%d",&nmx)==-1)
         err = 1;
      if (qcd_readOpenQCDLine("istop","%d",&istop)==-1)
         err = 1;
      if (qcd_readOpenQCDLine("res","%lf",&res)==-1)
         err = 1;
      
   }

   if (NPROC>1)
   {     
      MPI_Bcast(&err,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ids,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nkv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&isolv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nmr,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&ncy,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&nmx,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&istop,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
   }
   if (err==1)
      return 1;
   else
   {
      spar->sol[isp].isolv = isolv;
      spar->sol[isp].ncy   = ncy;
      spar->sol[isp].nkv   = nkv;
      spar->sol[isp].nmr   = nmr;
      spar->sol[isp].nmx   = nmx;
      spar->sol[isp].res   = res;
      spar->sol[isp].istop = istop;
      spar->sol[isp].solver= ids;
      return 0;
   }
}

void solver_printSolvers(solver_parameters *par)
{
   int i;
   int ps=0, pd=0;
   for (i=0; i<solver_ISPMAX; i++)
   {
      switch( par->sol[i].solver )
      {
         case solver_CGNE:        
            printf("\nSolver %i\n",i);
            printf("   type:   CGNE\n");
            printf("   nmx   = %i\n",par->sol[i].nmx);
            printf("   istop = %i\n",par->sol[i].istop);
            printf("   res   = %e\n",par->sol[i].res);
            break;
         case solver_SAP_GCR:
            printf("\nSolver %i\n",i);
            printf("   type:   SAP_GCR\n");
            printf("   nkv   = %i\n",par->sol[i].nkv);
            printf("   isolv = %i\n",par->sol[i].isolv);
            printf("   nmr   = %i\n",par->sol[i].nmr);
            printf("   ncy   = %i\n",par->sol[i].ncy);
            printf("   nmx   = %i\n",par->sol[i].nmx);
            printf("   istop = %i\n",par->sol[i].istop);
            printf("   res   = %e\n",par->sol[i].res);
            ps=1;
            break;
         case solver_DFL_SAP_GCR:
            printf("\nSolver %i\n",i);
            printf("   type:   DFL_SAP_GCR\n");
            printf("   nkv   = %i\n",par->sol[i].nkv);
            printf("   isolv = %i\n",par->sol[i].isolv);
            printf("   nmr   = %i\n",par->sol[i].nmr);
            printf("   ncy   = %i\n",par->sol[i].ncy);
            printf("   nmx   = %i\n",par->sol[i].nmx);
            printf("   istop = %i\n",par->sol[i].istop);
            printf("   res   = %e\n",par->sol[i].res);
            ps=1;
            pd=1;
            break;
         case solver_MSCG:        
            printf("\nSolver %i\n",i);
            printf("   type:   MSCG\n");
            printf("   nmx   = %i\n",par->sol[i].nmx);
            printf("   istop = %i\n",par->sol[i].istop);
            printf("   res   = %e\n",par->sol[i].res);
            break;
      }
   }
   if (ps)
   {
      printf("\nSAP\n   bs= %i x %i x %i x %i\n",par->sap.bs[0],par->sap.bs[1],
                                                 par->sap.bs[2],par->sap.bs[3]);
   }
   if (pd)
   {
      printf("\nDeflation subspace\n");
      printf("   bs = %i x %i x %i x %i\n",par->dfl.bs[0], par->dfl.bs[1],
                                           par->dfl.bs[2], par->dfl.bs[3]);
      printf("   Ns = %i\n",par->dfl.Ns);
     
      printf("\nDeflation subspace generation\n");
      printf("   kappa = %0.15f\n",par->dfl.kappa);
      printf("   mu    = %0.15f\n",par->dfl.mu);
      printf("   ninv  = %i\n",par->dfl.ninv);
      printf("   nmr   = %i\n",par->dfl.nmr);
      printf("   ncy   = %i\n",par->dfl.ncy);
     
      printf("\nDeflation projection\n");
      printf("   nkv     = %i\n",par->dfl.nkv);
      printf("   nmx     = %i\n",par->dfl.nmx);
      printf("   res     = %e\n",par->dfl.res);
      printf("   nmx_gcr = %d\n",par->dfl.nmx_gcr);
      printf("   res_gcr = %.1e\n\n",par->dfl.res_gcr);

   }
}// end solver_printSolvers


void solver_getProcessGrid(int *p)
{
   p[0]=NPROC0;
   p[1]=NPROC1;
   p[2]=NPROC2;
   p[3]=NPROC3;
   return;
}

void solver_getLatticeSize(int *l)
{
   l[0]=N0;
   l[1]=N1;
   l[2]=N2;
   l[3]=N3;
   return;
}

void solver_getStatus(int *stat)
{
   memcpy(stat,status,4*sizeof(int));
   return;
}


int solver_free()
{
   /* free memory allocated by openQCD 
      TODO -->modules/utils/wspace.c ??*/
   free(status);
   return 0;   
}//end solver_free
