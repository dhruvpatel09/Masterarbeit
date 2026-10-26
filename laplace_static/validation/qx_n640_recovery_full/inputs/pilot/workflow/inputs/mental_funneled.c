/* mental.c
 *
 * Creates elementals with momenta and/or derivatives
 *
 *
 * Tomasz Korzec 2025, 
 *    distredelemshdf5 by Andreas Risch 2025,
 *    disteigvecshdf5 by Andreas Risch 2025,
 *    disteigvalshdf5 by Andreas Risch 2025
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <mpi.h>
#include <math.h>
#include <getopt.h>
#include <qcd.h>
#include <DistRedelemsHdf5Writer.h>
#include <DistEigvecsHdf5Writer.h>
#include <DistEigvecsHdf5Reader.h>
#include <DistEigvalsHdf5Writer.h>

#define mental_RELEASE "mental 1.0"

static struct
{
   int Nv;      /* number of eigenvectors */
   int nset;    /* number of elemental sets */
   int *nder;   /* nder[i] is the number of derivatives of elemental i */
   int **dermu; /* dermu[i][j] is the direction of the j'th derivative in elemental i */
   int *p2max;  /* max (integer) momentum squared for each set */
   
   /* Fields related to the file writer (up to 5 different dat files)*/
   DistRedelemsHdf5Writer writer[5];
   int wr_nd[5];
   int *wr_muscombs[5];
   int wr_np[5];
   int *wr_moms[5];
   int *wr_times;
   DistRedelemsHdf5_Metadata wr_mdat;
   DistRedelemsHdf5_Programdata wr_pdat;
   
   /* Fields related to the eigenvector reader  / writer */
   DistEigvecsHdf5Writer evec_writer;
   DistEigvecsHdf5Reader evec_reader;
   DistEigvalsHdf5Writer eval_writer;
   
   /* File names and folders */
   char nbase[qcd_MAX_STRING_LENGTH];
   char outbase[qcd_MAX_STRING_LENGTH];
   char in_file[qcd_MAX_STRING_LENGTH];
   char dat_files[5][qcd_MAX_STRING_LENGTH];
   char out_file[qcd_MAX_STRING_LENGTH];
   char end_file[qcd_MAX_STRING_LENGTH];
   char cnfg_dir[qcd_MAX_STRING_LENGTH];
   char dat_dir[qcd_MAX_STRING_LENGTH];
   char evec_dir[qcd_MAX_STRING_LENGTH];
   char eval_dir[qcd_MAX_STRING_LENGTH];
   char log_dir[qcd_MAX_STRING_LENGTH];
   
   
   int L[4]; /* global lattice size */
   int P[4]; /* process grid */
   
   /* cnfgs to measure on */
   int nc_first;
   int nc_last;
   int nc_step;
   
   /* flags */
   char compute_ev; // 0: load eigenvectors, 1: compute eigenvectors
   char write_ev;   // 1: write eigenvectors (only together with comp_ev=1)
   char hasDer;
   char needsAPE;
   char differentAPE;

   /* smearing parameters for 3D Laplace */
   int  laplace_nAPE;
   double laplace_alphaAPE;
   /* smearing parameters for cov. deriv. */
   int der_nAPE;
   double der_alphaAPE;
   
   /* parameters for Laplace eigenvalue problem */
} params;

static struct
{
   int np;
   int *p1;
   int *p2;
   int *p3;
   int *pp;
} tabs;

static struct
{
   qcd_spinorComponent3d *evec; /* the first params.Nv 3d Laplace eigenvectors */
   double *eval;                /* the first params.Nv 3d Laplace eigenvalues */
   qcd_complex_16 *elem;        /* elemental */
   char cnfg[qcd_MAX_STRING_LENGTH];
   int nc;
} data;

//MPI related
static int myid,numprocs;

static qcd_lanczosParams lpar;

static void qcd_errorLoc(int test,int no,const char * const name,const char * const format,...)
{
   va_list args;
   int myid;
   
   if (test != 0)
   {
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      fprintf(stderr,"\nError in %s on process %i:\n",name,myid);
      va_start(args,format);
      vfprintf(stderr,format,args);
      va_end(args);
      fprintf(stderr,"\nProgram aborted\n\n");
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD,no);
   }
}


static int cmp_meta(DistEigvecsHdf5_Metadata *m1, DistEigvecsHdf5_Metadata *m2)
{
   // return 0 if m1 and m2 are equal, 1 otherwise
   if( (m1->bc.type == m2->bc.type) &&
       (m1->bc.theta[0] == m2->bc.theta[0]) &&
       (m1->bc.theta[1] == m2->bc.theta[1]) &&
       (m1->bc.theta[2] == m2->bc.theta[2]) &&
       (m1->laplacesmear.itercnt == m2->laplacesmear.itercnt) &&
       (fabs(m1->laplacesmear.alpha - m2->laplacesmear.alpha)<1e-10) &&
       (strcmp(m1->laplacesmear.type, m2->laplacesmear.type)==0) )
   {
      return 0;
   } else
   {
      return 1;
   }
}

static void read_dirs(void)
{
   if (myid==0)
   {
      qcd_findOpenQCDSection("Run name");
      qcd_readOpenQCDLine("name","%s",params.nbase);
      qcd_readOpenQCDLineOpt("output",params.nbase,"%s",params.outbase);

      qcd_findOpenQCDSection("Directories");
      qcd_readOpenQCDLine("log_dir","%s",params.log_dir);
      qcd_readOpenQCDLine("cnfg_dir","%s",params.cnfg_dir);
      qcd_readOpenQCDLine("dat_dir","%s",params.dat_dir);
      qcd_readOpenQCDLine("evec_dir","%s",params.evec_dir);
      qcd_readOpenQCDLine("eval_dir","%s",params.eval_dir);
   }

   MPI_Bcast(params.nbase,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.outbase,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(params.log_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.cnfg_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.dat_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.evec_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.eval_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
}

static void read_lat_parms(void)
{
   if (myid==0)
   {
      qcd_findOpenQCDSection("Lattice");
      qcd_readOpenQCDLine("glob_size", "%d %d %d %d", &(params.L[0]),&(params.L[1]), &(params.L[2]), &(params.L[3]));
      qcd_readOpenQCDLine("proc_grid", "%d %d %d %d", &(params.P[0]),&(params.P[1]), &(params.P[2]), &(params.P[3]));
   }
   MPI_Bcast(params.L, 4, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(params.P, 4, MPI_INT, 0, MPI_COMM_WORLD);
}


static void setup_files(void)
{
   int i;
   int ret,nder;

   qcd_errorLoc(!qcd_isDirWritable(params.dat_dir), 1, "setup_files [mental.c]",
                                                 "dat directory is not writable");
      
   if( (params.compute_ev==1) && (params.write_ev==1))
   {
         qcd_errorRoot(!qcd_isDirWritable(params.evec_dir), 1, "setup_files [mental.c]",
                                                    "eigenvector directory is not writable");
      
         qcd_errorRoot(!qcd_isDirWritable(params.eval_dir), 1, "setup_files [mental.c]",
                                                    "eigenvalue directory is not writable");
         
         ret = DistEigvecsHdf5Writer_Init(&(params.evec_writer));
         qcd_errorLoc(ret!=0,1,"setup_files [mental.c]",
                 "DistEigvecsHdf5Writer_Init failed");
         
         ret = DistEigvalsHdf5Writer_Init(&(params.eval_writer));
         qcd_errorLoc(ret!=0,1,"setup_files [mental.c]",
                 "DistEigvalsHdf5Writer_Init failed");
   }
   
   if( params.compute_ev==0 )
   {
      //qcd_errorRoot(!qcd_isDirReadable(params.evec_dir), 1, "setup_files [mental.c]",
      //                                              "eigenvector directory is not readable");
      params.evec_reader.initialised_ = 0; // TODO: THAT SHOULD NOT BE NECESSARY!
      ret = DistEigvecsHdf5Reader_Init(&(params.evec_reader));
      qcd_errorLoc(ret!=0,1,"setup_files [mental.c]",
                 "DistEigvecsHdf5Reader_Init failed");
   }
   
   for(nder=0; nder<5; nder++)
   {
      if( params.wr_nd[nder]>0 )
      {
         ret = DistRedelemsHdf5Writer_Init(&(params.writer[nder]));
            qcd_errorLoc(ret!=0,1,"setup_files [mental.c]",
                 "DistRedelemsHdf5Writer_Init failed");
      }
   }
}


static void initDatFiles(qcd_geometry *geo, int ic)
{
   int ret=0, i,nder;
   int p2,np, nd;
   char fname[qcd_MAX_STRING_LENGTH];
   
   // if eigenvectors from file: get meta-data from one of the files
   if( params.compute_ev == 0 )
   {
      ret = snprintf(fname, qcd_MAX_STRING_LENGTH,"%s/%sn%i_evec_t%i.h5",params.evec_dir,params.nbase,params.nc_first,0+geo->Pos[0]*geo->lL[0]);
      qcd_errorLoc(ret>=qcd_MAX_STRING_LENGTH,1,"initDatFiles [mental.c]",
               "file name for eigenvectors too long");
      int P3[3];
      P3[0] = geo->L[1]/geo->lL[1];
      P3[1] = geo->L[2]/geo->lL[2];
      P3[2] = geo->L[3]/geo->lL[3];
      char *tmpstr;
      int tmpnc;
      int tmpt;
      int tmpNv;
      DistEigvecsHdf5_Metadata tmpmdata;
      DistEigvecsHdf5_Programdata tmppdata;
      ret = DistEigvecsHdf5Reader_OpenFile(&(params.evec_reader), geo->comm3d,
                                           fname, &tmpstr, &tmpnc, &(geo->L[1]),
                                           P3, &tmpt, &tmpNv, 
                                           &tmpmdata,
                                           &tmppdata);
      qcd_errorLoc(ret != 0,1,"initDatFiles [mental.c]",
              "DistEigvecsHdf5Reader_OpenFile failed.");
      qcd_errorLoc((tmpnc != params.nc_first) || (tmpt != geo->Pos[0]*geo->lL[0]) || (tmpNv < params.Nv)
                || (strcmp(tmpstr,params.nbase)!=0), 1, "initDatFiles [mental.c]",
                "Eigenvector file does not match run");
            
      params.wr_mdat.bc.type = tmpmdata.bc.type;
      params.wr_mdat.bc.theta[0] = tmpmdata.bc.theta[0];
      params.wr_mdat.bc.theta[1] = tmpmdata.bc.theta[1];
      params.wr_mdat.bc.theta[2] = tmpmdata.bc.theta[2];
      //TODO: do something about "laplacesmear.type"
      params.wr_mdat.laplacesmear.type = "APE";
      params.wr_mdat.laplacesmear.itercnt = tmpmdata.laplacesmear.itercnt;
      params.wr_mdat.laplacesmear.alpha = tmpmdata.laplacesmear.alpha;
      MPI_Bcast(&(params.wr_mdat), sizeof(DistEigvecsHdf5_Metadata), MPI_CHAR, 0, MPI_COMM_WORLD);
       //TODO: do something about "laplacesmear.type"
      params.wr_mdat.laplacesmear.type = "APE";
      if(myid==0)
      {
         printf("Eigenvectors were created by %s with input files like %s.\n",tmppdata.name, tmppdata.infile);
         printf("   bc=%i, theta=[%f, %f, %f]\n",tmpmdata.bc.type, tmpmdata.bc.theta[0],tmpmdata.bc.theta[1],tmpmdata.bc.theta[2]);
         printf("   Smearing: %s with n_iter=%i, alpha=%f\n",tmpmdata.laplacesmear.type,tmpmdata.laplacesmear.itercnt, tmpmdata.laplacesmear.alpha);
      }
      ret =DistEigvecsHdf5Reader_CloseFile(&(params.evec_reader));
            qcd_errorLoc(ret != 0,1,"initDatFiles [mental.c]",
              "DistEigvecsHdf5Reader_CloseFile failed.");
   }
   
   for(nder=0; nder<5; nder++)
   {
      if( params.wr_nd[nder]>0 )
      {
         ret=snprintf(fname,qcd_MAX_STRING_LENGTH,"%s/%s_der%i_n%i.h5",params.dat_dir,params.outbase,nder,ic);
         qcd_errorLoc(ret>=qcd_MAX_STRING_LENGTH,1,"initDatFiles [mental.c]",
              "file name for hdf5 file too long");
         //check if exists
         
         ret = DistRedelemsHdf5Writer_CreateFile(&(params.writer[nder]),
                                                 fname, params.nbase, ic,
                                                 nder, params.wr_nd[nder], 
                                                 params.wr_muscombs[nder], params.wr_np[nder],
                                                 params.wr_moms[nder], geo->L[0], params.wr_times, params.Nv,
                                                 &(params.wr_mdat), &(params.wr_pdat));
         qcd_errorLoc(ret!=0,1,"initDatFiles [mental.c]",
              "DistRedelemsHdf5Writer_CreateFile failed");
         
         if( geo->myid == 0)
         {
            printf("initialized dat file: %s\n",fname);
         }
      }
   }
   
   if( params.compute_ev && params.write_ev)
   {
      ret=snprintf(fname,qcd_MAX_STRING_LENGTH,"%s/%s_n%i_eval.h5",params.eval_dir,params.outbase,ic);
         qcd_errorLoc(ret>=qcd_MAX_STRING_LENGTH,1,"initDatFiles [mental.c]",
              "file name for eigenvalues file too long");
      ret = DistEigvalsHdf5Writer_CreateFile(&(params.eval_writer),
                                             fname, params.nbase, ic,
                                             geo->L[0], params.wr_times, params.Nv,
                                             (DistEigvalsHdf5_Metadata *) &(params.wr_mdat), (DistEigvalsHdf5_Programdata *) &(params.wr_pdat));
      qcd_errorLoc(ret!=0,1,"initDatFiles [mental.c]",
              "DistEigvalsHdf5Writer_CreateFile failed");
   }
}

static void closeDatFiles()
{
   int nder;
   int ret;
   for(nder=0; nder<5; nder++)
   {
      if( params.wr_nd[nder]>0 )
      {
         ret = DistRedelemsHdf5Writer_CloseFile(&(params.writer[nder]));
         qcd_errorLoc(ret!=0,1,"closeDatFiles [mental.c]",
              "DistRedelemsHdf5Writer_CloseFile failed");
      }
   }
   if( params.compute_ev && params.write_ev)
   {
      ret = DistEigvalsHdf5Writer_CloseFile(&(params.eval_writer));
      qcd_errorLoc(ret!=0,1,"closeDatFiles [mental.c]",
              "DistEigvalsHdf5Writer_CloseFile failed");
   }
}


static void printhelp(void)
{
   printf("usage:\n");
   printf("mental -i <input-file> -o <output-file>\n");
   printf("   run with input file <input-file> and write output to <output-file>\n\n");
   printf("mental -i <input-file>\n");
   printf("   run with input file <input-file> and write output to stdout\n\n");
   printf("mental -t\n");
   printf("   write an input-file template to stdout\n\n");
   printf("mental <anything else>\n");
   printf("   print this help message\n");
}

static void printversion(void)
{
   printf("mental version <%s>\n\n",mental_RELEASE);
}

static void printtemplate(void)
{
   //print out a template input-file
   
}

void print_p_tab()
{
   int i=0;
   if (myid==0)
   {
      printf("\nMomenta table for set with largest pmax\n");
      printf("p1   p2   p3   p^2\n");
      for(i=0; i<tabs.np; i++)
         printf("%+3d %+3d %+3d %4d\n",tabs.p1[i],tabs.p2[i],tabs.p3[i],tabs.pp[i]);
   }
}

static void setup_out()
{
   FILE *flog;
   flog=fopen(params.out_file,"r");
   qcd_errorRoot(flog!=NULL,1,"setup_out [mental.c]",
              "output file %s already exists!",params.out_file);

   flog = freopen(params.out_file,"a",stdout);
   qcd_errorRoot(flog==NULL,1,"setup_out [mental.c]",
              "could not open %s for appending",params.out_file);
   flog = freopen(params.out_file,"a",stderr);
   qcd_errorRoot(flog==NULL,1,"check_out [mental.c]",
              "could not open %s for appending",params.out_file);
   return;
}




static void make_tables(qcd_geometry *geo)
{
   int p1,p2,p3,pmax,p2max=0;
   int i,j,k,tmp;
   int is_in=0;
   int mu1,mu2,mu3,mu4;
   int *nder;
   int **dermu;
   int nset;
   int *p2max_ptr;
   int threepow[5] = {1, 3, 9, 27, 81};
   int pind;

   // find largest p2max
   for(i=0; i<params.nset; i++)
   {
      if(params.p2max[i] > p2max)
         p2max = params.p2max[i];
   }
   // unsorted momenta table
   pmax = (int) sqrt((double) p2max);
   tabs.np = 0;
   for( p1=pmax; p1>=-pmax; p1--)
   for( p2=pmax; p2>=-pmax; p2--)
   for( p3=pmax; p3>=-pmax; p3--)
   {
      if(p1*p1+p2*p2+p3*p3<=p2max)
      {
         tabs.np++;
      }
   }
   tabs.p1 = malloc(tabs.np*sizeof(int));
   tabs.p2 = malloc(tabs.np*sizeof(int));
   tabs.p3 = malloc(tabs.np*sizeof(int));
   tabs.pp = malloc(tabs.np*sizeof(int));
   tabs.np=0;
   for( p1=pmax; p1>=-pmax; p1--)
   for( p2=pmax; p2>=-pmax; p2--)
   for( p3=pmax; p3>=-pmax; p3--)
   {
      if(p1*p1+p2*p2+p3*p3<=p2max)
      {
         // check symmetry
         is_in=0;
         for(i=0; i<tabs.np; i++)
         {
            if( (p1 == -tabs.p1[i]) &&
                (p2 == -tabs.p2[i]) &&
                (p3 == -tabs.p3[i]) )
            {
               // this momentum will by realized by symmetry
               // no computation needed
               is_in =1;
            }
         }
         if(!is_in)
         {
            tabs.p1[tabs.np] = p1;
            tabs.p2[tabs.np] = p2;
            tabs.p3[tabs.np] = p3;
            tabs.pp[tabs.np] = p1*p1+p2*p2+p3*p3;
            tabs.np++;
         }
      }
   }
   // sort tables, increasing p2
   for(i=1; i<tabs.np; i++)
   {
      for(j=i-1; j>=0; j--)
      {
         if( tabs.pp[j+1]<tabs.pp[j] )
         {
            //swap j+1 and j
            tmp = tabs.p1[j+1]; tabs.p1[j+1] = tabs.p1[j]; tabs.p1[j] = tmp;
            tmp = tabs.p2[j+1]; tabs.p2[j+1] = tabs.p2[j]; tabs.p2[j] = tmp;
            tmp = tabs.p3[j+1]; tabs.p3[j+1] = tabs.p3[j]; tabs.p3[j] = tmp;
            tmp = tabs.pp[j+1]; tabs.pp[j+1] = tabs.pp[j]; tabs.pp[j] = tmp;
         }else
         {   
            break; //j-loop          
         }
      }
   }
   
   // TODO: pre-compute sin, cos up to pmax*2pi/L
   
   // expand sets: include all combinations up to nder
   nset = 0;
   for(i=0; i<params.nset; i++)
   {
      nset += threepow[params.nder[i]];
   }
   nder      = malloc(nset*sizeof(int));
   p2max_ptr = malloc(nset*sizeof(int));
   dermu     = malloc(nset*sizeof(int*));
   
   nset=0;
   for(i=0; i<params.nset; i++)
   {
      switch (params.nder[i])
      {
         case 0: nder[nset] = 0;
                 p2max_ptr[nset] = params.p2max[i];
                 nset++;
                 break;
         case 1: for(mu1=1; mu1<=3; mu1++)
                 {
                    nder[nset] = 1;
                    p2max_ptr[nset] = params.p2max[i];
                    dermu[nset] = malloc(nder[nset]*sizeof(int));
                    dermu[nset][0] = mu1;
                    nset++;
                 }
                 break;
         case 2: for(mu1=1; mu1<=3; mu1++)
                 for(mu2=1; mu2<=3; mu2++)
                 {
                    nder[nset] = 2;
                    p2max_ptr[nset] = params.p2max[i];
                    dermu[nset] = malloc(nder[nset]*sizeof(int));
                    dermu[nset][0] = mu1;
                    dermu[nset][1] = mu2;
                    nset++;
                 }
                 break;
         case 3: for(mu1=1; mu1<=3; mu1++)
                 for(mu2=1; mu2<=3; mu2++)
                 for(mu3=1; mu3<=3; mu3++)
                 {
                    nder[nset] = 3;
                    p2max_ptr[nset] = params.p2max[i];
                    dermu[nset] = malloc(nder[nset]*sizeof(int));
                    dermu[nset][0] = mu1;
                    dermu[nset][1] = mu2;
                    dermu[nset][2] = mu3;
                    nset++;
                 }
                 break;
         case 4: for(mu1=1; mu1<=3; mu1++)
                 for(mu2=1; mu2<=3; mu2++)
                 for(mu3=1; mu3<=3; mu3++)
                 for(mu4=1; mu4<=3; mu4++)
                 {
                    nder[nset] = 4;
                    p2max_ptr[nset] = params.p2max[i];
                    dermu[nset] = malloc(nder[nset]*sizeof(int));
                    dermu[nset][0] = mu1;
                    dermu[nset][1] = mu2;
                    dermu[nset][2] = mu3;
                    dermu[nset][3] = mu4;
                    nset++;
                 }
                 break;
      }//end switch
   }// end in-file-sets-loop
   //replace preliminary set-table by final one:
   free(params.dermu);
   params.dermu = dermu;
   free(params.nder);
   params.nder = nder;
   free(params.p2max);
   params.p2max = p2max_ptr;
   params.nset = nset;
   
   
   /* prepare tables for the writer */
   for(i=1; i<5; i++)
   {
      for(j=0; j<params.nset; j++)
      {
         if(params.nder[j]==i)
         {
            //j-derivative elemental exists. Allocate muscombs table
            params.wr_muscombs[i] = malloc(threepow[i]*i*sizeof(int));
            break;
         }
      }
   }
   
   
   for(i=0; i<5; i++)
   {
      params.wr_nd[i] = 0;
      for( j=0; j<params.nset; j++ )
      {
         if(params.nder[j]==i)
         {
            //add derivatives combo to list
            for(k=0; k<i; k++)
               params.wr_muscombs[i][i*params.wr_nd[i]+k] = params.dermu[j][k];
            params.wr_nd[i]++;
            p2max = params.p2max[j];
         }
      }
      params.wr_np[i] = 0;
      if(params.wr_nd[i]>0)
      {
         // count momenta in range
         for(j=0; j<tabs.np; j++)
            if(tabs.pp[j] <= p2max)
               params.wr_np[i]++;
         // allocate writer table
         params.wr_moms[i] = malloc(3*params.wr_np[i]*sizeof(int));
         
         // fill writer table
         pind=0;
         for(j=0; j<tabs.np; j++)
            if(tabs.pp[j] <= p2max)
            {
               params.wr_moms[i][3*pind]   = tabs.p1[j];
               params.wr_moms[i][3*pind+1] = tabs.p2[j];
               params.wr_moms[i][3*pind+2] = tabs.p3[j];
               pind++;
            }
      }
   }
   params.wr_times = malloc(geo->L[0]*sizeof(int));
   for(i=0; i<geo->L[0]; i++)
      params.wr_times[i] = i;
   // prepare meta-data for writer
   params.wr_mdat.bc.type = 0; //TODO: ??
   params.wr_mdat.bc.theta[0] = 0.0;
   params.wr_mdat.bc.theta[1] = 0.0;
   params.wr_mdat.bc.theta[2] = 0.0;
   params.wr_mdat.laplacesmear.type= "APE"; // TODO: ??
   params.wr_mdat.laplacesmear.itercnt= params.laplace_nAPE;
   params.wr_mdat.laplacesmear.alpha= params.laplace_alphaAPE;
   params.wr_mdat.derivativesmear.type= "APE"; // TODO: ??
   params.wr_mdat.derivativesmear.itercnt= params.der_nAPE;
   params.wr_mdat.derivativesmear.alpha= params.der_alphaAPE;
   
   params.wr_pdat.name=mental_RELEASE;
   params.wr_pdat.infile=params.in_file;
}


static void read_infile(void)
{
   int lf,i,j;
   char tmpstr[qcd_MAX_STRING_LENGTH];
   // command line option positions/flags
   int ifile;

   FILE *fin;

   if(myid==0)
   {
      qcd_errorRoot(qcd_isBigEndian()==-1,1,"read_infile [mental.c]",
                 "Machine has unknown endianness");

      fin=freopen(params.in_file,"r",stdin);
      qcd_errorRoot(fin==NULL,1,"read_infile [mental.c]",
                 "Unable to open input file");
   }
   read_dirs();
   setup_files();
   read_lat_parms();

   if(myid==0)
   {
      qcd_findOpenQCDSection("Configurations");
      qcd_readOpenQCDLine("first","%d",&(params.nc_first));
      qcd_readOpenQCDLine("last","%d",&(params.nc_last));
      qcd_readOpenQCDLine("step","%d",&(params.nc_step));

      qcd_errorRoot((params.nc_last<params.nc_first)||
                 (params.nc_step<1)||
                 (((params.nc_last-params.nc_first)%params.nc_step)!=0),1,
                 "read_dirs [mental.c]","Improper configuration range");
      
      qcd_findOpenQCDSection("Measurements");
      qcd_readOpenQCDLine("nset","%d",&(params.nset));
      qcd_readOpenQCDLine("Nv","%d",&(params.Nv));
      qcd_readOpenQCDLine("compute_ev","%d",&(params.compute_ev));
      qcd_readOpenQCDLine("write_ev","%d",&(params.write_ev));
      
      
      //TODO: sanity check for these parameters
   }
   MPI_Bcast(&(params.nc_first),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.nc_last),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.nc_step),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.nset),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.Nv),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.compute_ev),1,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.write_ev),1,MPI_CHAR,0,MPI_COMM_WORLD);

   //allocate memory for all set parameters
   params.p2max = malloc(params.nset*sizeof(int));
   params.nder  = malloc(params.nset*sizeof(int));
   params.hasDer=0;
   for(i=0; i<params.nset; i++)
   {
      if(myid==0)
      {
         sprintf(tmpstr, "set %d",i);
         qcd_findOpenQCDSection(tmpstr);
         qcd_readOpenQCDLine("p2max","%d",&(params.p2max[i]));
         qcd_readOpenQCDLine("nder","%d",&(params.nder[i]));
      }
      MPI_Bcast(&(params.p2max[i]), 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&(params.nder[i]), 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(params.nder[i]>0) params.hasDer = 1;
      
      qcd_errorRoot(params.nder[i]>4,1,"read_infile [mental.c]","maximum nder is 4");
      //check if the same nder was used in a different set (illegal)
      for(j=0; j<i; j++)
      {
         qcd_errorRoot(params.nder[i]==params.nder[j],1,"read_infile [mental.c]","Each nder may appear only once");
      }
   }
   
   qcd_defaultLanczosParams(&lpar);
   if( myid==0 )
   {
      qcd_findOpenQCDSection("Measurements");
      if(params.hasDer) //read APE parameters for derivatives
      {
         qcd_readOpenQCDLine("der_nAPE","%d",&(params.der_nAPE));
         qcd_readOpenQCDLine("der_alphaAPE","%lf",&(params.der_alphaAPE));
      }
      if(params.compute_ev) //read APE parameters for derivatives
      {
         qcd_readOpenQCDLine("laplace_nAPE","%d",&(params.laplace_nAPE));
         qcd_readOpenQCDLine("laplace_alphaAPE","%lf",&(params.laplace_alphaAPE));
         // read lanczos parameters
         qcd_findOpenQCDSection("Lanczos");
         qcd_readOpenQCDLineOpt("res", "1e-10", "%lf", &(lpar.res));
         qcd_readOpenQCDLineOpt("redTol", "1e-12", "%lf", &(lpar.redTol));
         qcd_readOpenQCDLineOpt("deg", "4", "%d", &(lpar.deg));
         qcd_readOpenQCDLineOpt("upper", "10.0", "%lf", &(lpar.upper));
         qcd_readOpenQCDLineOpt("tolReortho", "1e-10", "%lf", &(lpar.tolReortho));
         qcd_readOpenQCDLineOpt("maxRestarts", "30", "%d", &(lpar.maxRestarts));
         qcd_readOpenQCDLineOpt("nLanczosFactor", "6", "%d", &(lpar.nLanczosFactor));
         qcd_readOpenQCDLineOpt("redFreqFactor", "1", "%d", &(lpar.redFreqFactor));
      }
   }
   if( params.hasDer )
   {
      MPI_Bcast(&(params.der_nAPE),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&(params.der_alphaAPE),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   if( params.compute_ev )
   {
      MPI_Bcast(&(params.laplace_nAPE),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&(params.laplace_alphaAPE),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.res),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.redTol),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.deg),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.upper),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.tolReortho),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.maxRestarts),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.nLanczosFactor),1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&(lpar.redFreqFactor),1,MPI_INT,0,MPI_COMM_WORLD);
   }
   if(params.hasDer || params.compute_ev)
      params.needsAPE = 1;
   else
      params.needsAPE = 0;
   params.differentAPE = 1;
   if( params.hasDer && params.compute_ev && (params.der_nAPE == params.laplace_nAPE) && (params.der_alphaAPE) == (params.laplace_alphaAPE) )
      params.differentAPE = 0;
   if (myid==0)
   {
      fclose(fin);
   }

   return;
}

static void print_params(qcd_geometry *geo)
{
   int i,j;
   if(myid==0)
   {
      printversion();
      printf(" Global lattice: %3i x %3i x %3i x %3i\n",geo->L[0],geo->L[1],geo->L[2],geo->L[3]);
      printf(" Local  lattice: %3i x %3i x %3i x %3i\n",geo->lL[0],geo->lL[1],geo->lL[2],geo->lL[3]);
      printf(" Process grid:   %3i x %3i x %3i x %3i\n",geo->L[0]/geo->lL[0],geo->L[1]/geo->lL[1],geo->L[2]/geo->lL[2],geo->L[3]/geo->lL[3]);
      //print dirs & lat-print_params
      printf("Dirs:\n");
      printf("   nbase   = %s\n",params.nbase);
      printf("   output  = %s\n",params.outbase);
      
      printf("Configurations %d:%d:%d\n",params.nc_first, params.nc_step, params.nc_last);
      if(params.compute_ev && params.write_ev)
         printf("Number of laplace eigenvectors: %d, to be computed and written\n",params.Nv);
      else if(params.compute_ev)
         printf("Number of laplace eigenvectors: %d, to be computed\n",params.Nv);
      else
         printf("Number of laplace eigenvectors: %d, to be read from disk\n",params.Nv);

      if(params.compute_ev)
         printf("Gauge fields in Laplace operator will be smeared with: n_APE=%d, alpha_APE=%f\n",params.laplace_nAPE,params.laplace_alphaAPE);
      printf("Gauge fields in derivatives will be smeared with:      n_APE=%d, alpha_APE=%f\n",params.der_nAPE,params.der_alphaAPE);
   
      if(params.compute_ev)
      {
         printf("Lanczos Parameters\n");
         printf("   res =            %e\n", lpar.res);
         printf("   redTol =         %e\n", lpar.redTol);
         printf("   deg =            %d\n", lpar.deg);
         printf("   upper =          %lf\n", lpar.upper);
         printf("   tolReortho =     %e\n", lpar.tolReortho);
         printf("   maxRestarts =    %d\n", lpar.maxRestarts);
         printf("   nLanczosFactor = %d\n", lpar.nLanczosFactor);
         printf("   redFreqFactor =  %d\n\n", lpar.redFreqFactor);
      }
      
      printf("Will compute %i sets of elementals\n",params.nset);
      for(i=0; i<params.nset; i++)
      {
         printf("Set %d:\n",i);
         printf("   p2max=%d\n",params.p2max[i]);
         printf("   nder =%d\n",params.nder[i]);
         if(params.nder[i]>0)
         {
            printf("   mu   =");
            for(j=0; j<params.nder[i]; j++)
               printf("%d ",params.dermu[i][j]);
            printf("\n");
         }
      }
   }
   
   print_p_tab();
}

static int check_endflag()
{
   FILE *fend;
   int iend;

   if (myid==0)
   {
      fend=fopen(params.end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(params.end_file);
         printf("End flag set, run stopped\n\n");
         iend = 1;
      }
      else
         iend = 0;
   }
   MPI_Bcast(&iend,1,MPI_INT,0,MPI_COMM_WORLD);
   return iend;
}


static int get_muidx(int ls, int nder)
{
   //return mu-combination index in the writer table for set number ls
   int i,j,t;
   for(i=0; i < params.wr_nd[nder]; i++)
   {
      t = 1;
      for(j=0; j<nder; j++)
      {
         if( params.wr_muscombs[nder][nder*i+j] != params.dermu[ls][j] )
            t=0;
      }
      if( t )
      {
         return(i);
      }
   }
   qcd_errorLoc(1,1,"get_muidx [mental.c]","mu combination not in writer table");
   return -1;
}


static int  get_pidx(int lp, int nder)
{
   //return p-index in writer table for momentum with index lp in tabs
   int i;
   for(i = 0; i < params.wr_np[nder]; i++)
   {
      if( (params.wr_moms[nder][3*i]   == tabs.p1[lp]) &&
          (params.wr_moms[nder][3*i+1] == tabs.p2[lp]) &&
          (params.wr_moms[nder][3*i+2] == tabs.p3[lp]) )
      {
         return(i);
      }
   }
   qcd_errorLoc(1,1,"get_pidx [mental.c]","Momentum not in writer table");
   return -1;
}


int main(int argc,char* argv[])
{
   char has_o=0;                               // true if optional output-file was specified
   char has_i=0;                               // true if input-file was specified

   double   plaq;                              //average plaquette

   //loop variables
   int i, j, k, l, nc, lt, lp, ls, lx,ly,lz, gx,gy,gz;
   
   char fname[qcd_MAX_STRING_LENGTH];

   // gauge fields to actually allocate
   qcd_gaugeField u1; u1.initialized=0;
   qcd_gaugeField u2; u2.initialized=0;
   qcd_gaugeField u3; u3.initialized=0;
   // pointers to u1,u2,u3, actual locations change depending on smearing
   qcd_gaugeField *u_ptr;
   qcd_gaugeField *uAPEl_ptr;
   qcd_gaugeField *uAPEd_ptr;
   qcd_gaugeField *utmp_ptr;
   qcd_gaugeField *utmp2_ptr;
   // 3d gauge field
   qcd_gaugeField3d u3d;
   

   //laplace eigenvectors
   qcd_spinorComponent3d *v, vp, Dv, vtmp;

   //laplace eigenvalues
   double *lambda;
   
   //geometry parameters
   qcd_geometry geo;
   double theta[4]={0.,0.,0.,0.};              //boundary conditions, this program supports only 0

   //timing
   double time1, time2;
   
   //
   double nrm,phi;
   
   // phase factors
   qcd_complex_16 eipx, eipy, eipz, eipxy, eipxyz;

   // elemental data
   qcd_complex_16 *redelem;
   
   int muidx, pidx, tidx, wrflag, ret;
   
   // lanczos
   qcd_lanczosStatus lstat;
   int lret;
   
   // for the writer
   int wr_mom[3];
   int *wr_mus;
   DistEigvecsHdf5_ComplexDouble *disteigvecs;
   
   size_t Ncp;
   
   
   /////////////////////////////////////////////////////////////////////////////
   // set up MPI                                                              //
   /////////////////////////////////////////////////////////////////////////////

   /* Standalone mental: same MPI thread request as solver_mpiInit. */
   int mpi_provided, mpi_min_provided, mpi_init_rank;
   int mpi_rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,
                                &mpi_provided);
   if (mpi_rc != MPI_SUCCESS) {
      fprintf(stderr, "MPI_Init_thread failed: %d\n", mpi_rc);
      exit(EXIT_FAILURE);
   }
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_init_rank);
   MPI_Allreduce(&mpi_provided, &mpi_min_provided, 1, MPI_INT,
                 MPI_MIN, MPI_COMM_WORLD);
   if (mpi_min_provided < MPI_THREAD_FUNNELED) {
      if (mpi_init_rank == 0)
         fprintf(stderr, "MPI_THREAD_FUNNELED is not supported on all ranks\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      exit(EXIT_FAILURE);
   }
   if (mpi_init_rank == 0) {
      printf("MPI_INIT requested=FUNNELED provided_min=%d status=PASS\n",
             mpi_min_provided);
      fflush(stdout);
   }
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs); // number of processes
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);     // each process gets its ID

   /////////////////////////////////////////////////////////////////////////////
   // parse command line                                                      //
   /////////////////////////////////////////////////////////////////////////////
   int option_index = 0;
   int oc;
   static struct option long_options[] =
   {
      {"input",    1, NULL, 'i'},
      {"output",    1, NULL, 'o'},
      {"help",    0, NULL, 'h'},
      {"version", 0, NULL, 'v'},
      {"template", 0, NULL, 't'},
      {NULL,      0, NULL, 0}
   };

   /* parse options */
   do
   {
      /* getopt_long stores the option index here.   */
      oc = getopt_long (argc, argv, "i:o:hvt",
             long_options, &option_index);

      switch (oc)
      {
         case 'v':
            if(myid==0)
               printversion();
            MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
         case 'h':
            if(myid==0)
            {
               printversion();
               printhelp();
            }
            MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
         case 't':
            if(myid==0)
               printtemplate();
            MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS);
         case 'i':
            has_i = 1;
            memcpy(params.in_file,optarg,strlen(optarg)+1);
            if(myid==0)
               printf("input-file: %s\n",params.in_file);
            break;
         case 'o':
            memcpy(params.out_file,optarg,strlen(optarg)+1);
            has_o = 1;
            if(myid==0)
               printf("output-file: %s\n",params.out_file);
            break;
         case '?':
            if(myid==0)
               printhelp();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            break;
         case -1: break;
      }
   }while(oc != -1 );
   
   if ( !has_i )
   {
      if(myid==0)
      {
         fprintf(stderr,"Error, no input file specified!\n");
         printhelp();
      }
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
   
   /////////////////////////////////////////////////////////////////////////////
   // set up output file (log-file)                                           //
   /////////////////////////////////////////////////////////////////////////////   
   if( has_o )
   {
      // re-direct all output of process 0 to out-file
      if( myid==0 )
      {
         setup_out();
      }
   }

   /////////////////////////////////////////////////////////////////////////////
   // read input file                                                         //
   /////////////////////////////////////////////////////////////////////////////   
   read_infile();
   
   /////////////////////////////////////////////////////////////////////////////
   // set up qcdlib                                                           //
   /////////////////////////////////////////////////////////////////////////////
   if(qcd_initGeometry(&geo,params.L,params.P, theta, myid, numprocs)) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   if(qcd_initEO(&geo)) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   
   /////////////////////////////////////////////////////////////////////////////
   // prepare momenta list, lookup-tables, out-file-tables                    //
   /////////////////////////////////////////////////////////////////////////////
   make_tables(&geo);
   
   /////////////////////////////////////////////////////////////////////////////
   // set up files, write headers                                             //
   /////////////////////////////////////////////////////////////////////////////
   setup_files();
   
   /////////////////////////////////////////////////////////////////////////////
   // print log-file header                                                   //
   /////////////////////////////////////////////////////////////////////////////
   print_params(&geo);
   
   /////////////////////////////////////////////////////////////////////////////
   // allocate memory                                                         //
   /////////////////////////////////////////////////////////////////////////////
   j = qcd_initGaugeField(&u1,&geo);
   u_ptr = &u1;
   if( params.compute_ev )
   {
      if( params.laplace_nAPE != 0 )
      {
         j += qcd_initGaugeField(&u2,&geo);
         uAPEl_ptr = &u2;
      }
      else
         uAPEl_ptr = u_ptr;
   }
   if( params.hasDer )
   {
      if( params.differentAPE )
      {
         if( params.der_nAPE !=0 )
         {
            j += qcd_initGaugeField(&u3,&geo);
            uAPEd_ptr = &u3;
         }
         else
            uAPEd_ptr = u_ptr;
      } else
      {
         uAPEd_ptr = uAPEl_ptr;  
      }
   }
   j += qcd_initGaugeField3d(&u3d,&geo);
   qcd_errorLoc(j,1,"main [mental.c]","Error in gauge field initialization");
   lambda = (double*) malloc(params.Nv * sizeof(double));
   qcd_errorLoc(lambda==NULL,1,"main [mental.c]","Out of memory for eigenvalues");
   
   /////////////////////////////////////////////////////////////////////////////
   // initialize 3d laplace eigenvectors                                      //
   /////////////////////////////////////////////////////////////////////////////
   v = (qcd_spinorComponent3d*) malloc(params.Nv*sizeof(qcd_spinorComponent3d));
   qcd_errorLoc(v==NULL,1,"main [mental.c]","Error in eigenvector initialization");
   j = 0;
   for(i=0; i<params.Nv; i++)
      j += qcd_initSpinorComponent3d(&(v[i]), &geo);
   j += qcd_initSpinorComponent3d(&vp, &geo);
   j += qcd_initSpinorComponent3d(&Dv, &geo);
   j += qcd_initSpinorComponent3d(&vtmp, &geo);
   qcd_errorLoc(j>0,1,"main [mental.c]","Error in eigenvector initialization");
   
   if(myid==0){printf("eigenvectors initialized...\n"); fflush(stdout);}
   
   /////////////////////////////////////////////////////////////////////////////
   // initialize memory for reduced elementals                                //
   /////////////////////////////////////////////////////////////////////////////
   
   redelem = malloc(params.Nv*params.Nv*params.nset*sizeof(qcd_complex_16));
   qcd_errorLoc(redelem==NULL,1,"main [mental.c]","Out of memory for elementals",myid);
   
   //if(myid==0){printf("redelem initialized...\n"); fflush(stdout);}
   
   if( params.write_ev || !params.compute_ev )
   {
      disteigvecs = (DistEigvecsHdf5_ComplexDouble*) malloc(geo.lV3*3*params.Nv*sizeof(DistEigvecsHdf5_ComplexDouble));
      qcd_errorLoc(disteigvecs==NULL,1,"main [mental.c]","Out of memory for disteigs");
      //if(myid==0){printf("disteigvecs initialized...\n"); fflush(stdout);}
   }
   if(myid==0){printf("\nmemory allocated\n"); fflush(stdout);}
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);

   /*
    * Calibration mode:
    * compute the requested Laplace eigenpairs, report lambda_Nv for every
    * time-slice, and skip reduced-elemental construction.
    *
    * Normal behaviour is unchanged unless MENTAL_EVAL_ONLY is non-zero.
    */
   const char *eval_only_env = getenv("MENTAL_EVAL_ONLY");
   const int eval_only =
      (eval_only_env != NULL && atoi(eval_only_env) != 0);

   if(myid==0 && eval_only)
   {
      printf("MENTAL_EVAL_ONLY enabled: eigenvalue calibration only.\n");
      fflush(stdout);
   }
   
   /////////////////////////////////////////////////////////////////////////////
   // main loop                                                               //
   /////////////////////////////////////////////////////////////////////////////
   for (nc=params.nc_first; nc<=params.nc_last; nc+=params.nc_step)
   {  
      if( check_endflag() )
         break;
      
      j=snprintf(fname,qcd_MAX_STRING_LENGTH,"%s/%sn%i",params.cnfg_dir,params.nbase,nc);
      qcd_errorLoc(j>=qcd_MAX_STRING_LENGTH,1,"main [mental.c]",
              "file name for gauge configuration too long");
      if(myid==0) {printf("\nProcessing gauge field %s\n",fname); fflush(stdout);}

      time1 = MPI_Wtime();
      //read in the gauge field in "CERN" format
      qcd_errorLoc(
         qcd_getGaugeField(fname,qcd_GF_OPENQCD,u_ptr), 1,
         "main [mental.c]","Error reading gauge field!\n");
      time2 = MPI_Wtime();
      if(myid==0) {printf("gauge-field loaded in %fsek\n",time2-time1); fflush(stdout);}

      plaq = qcd_calculatePlaquette(u_ptr);
      if(myid==0) {printf("plaquette (torus def) = %e\n",plaq); fflush(stdout);}

      
      
      initDatFiles(&geo, nc);
      
      // if necessary: smear gauge field
      if( params.compute_ev)
      {
         if(params.laplace_nAPE !=0)
         {
            if(params.hasDer && params.differentAPE && (params.der_nAPE !=0))
            {
               // use uAPE_d as temporary space
               utmp_ptr = uAPEd_ptr;
               qcd_copyGaugeField(utmp_ptr, u_ptr);  
            } else
            {
               // work on original u, overwriting it
               utmp_ptr = u_ptr;
            }
            time1=MPI_Wtime();
            for( i=0; i<params.laplace_nAPE; i++ )
            {
               qcd_apeSmear3d(uAPEl_ptr, utmp_ptr, params.laplace_alphaAPE);
               if( i != params.laplace_nAPE-1 )
               {
                  utmp2_ptr = utmp_ptr;
                  utmp_ptr = uAPEl_ptr;
                  uAPEl_ptr = utmp2_ptr;
               }
            }
            time2=MPI_Wtime();
            if(myid==0)
            {
               printf("APE smearing for cov. Laplace done in %e sec\n",time2-time1); fflush(stdout);
            }
         }else
         {
            uAPEl_ptr = u_ptr;
         }
      }
      if( params.hasDer )
      {
         if( params.differentAPE == 0)
         {
            uAPEd_ptr = uAPEl_ptr;
         }else
         if(params.der_nAPE !=0)
         {
            time1=MPI_Wtime();
            // work on original u, overwriting it
            utmp_ptr = u_ptr;
            for( i=0; i<params.der_nAPE; i++ )
            {
               qcd_apeSmear3d(uAPEd_ptr, utmp_ptr, params.der_alphaAPE);
               if( i != params.der_nAPE-1 )
               {
                  utmp2_ptr = utmp_ptr;
                  utmp_ptr = uAPEd_ptr;
                  uAPEd_ptr = utmp2_ptr;
               }
            }
            time2=MPI_Wtime();
            if(myid==0)
            {
               printf("APE smearing for cov. derivatives done in %e sec\n",time2-time1); fflush(stdout);
            }
         } else
         {
            uAPEd_ptr = u_ptr;
         }
      }

      for( lt=0; lt<geo.lL[0]; lt++ )
      {
         // read/compute eigenvectors
         if( params.compute_ev )
         {
            time1=MPI_Wtime();
            qcd_copyGaugeField3dGaugeField(&u3d, uAPEl_ptr, lt+geo.Pos[0]*geo.lL[0]);
            lret=qcd_lanczosCheby(lambda, v, &lstat, &u3d, params.Nv, &lpar);
            if(lret != 0)
            {
               fprintf(stderr,"Error! Lanczos on time-slice %i did not converge!\n",lt+geo.Pos[0]*geo.lL[0]);
               MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            
            time2 = MPI_Wtime()-time1;
            MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, 0, geo.comm3d);
            if(myid==geo.root3d[lt+geo.lL[0]*geo.Pos[0]])
            {
               printf("Eigenvalues/vectors on x0=%i computed in %e sec (max)\n",lt+geo.lL[0]*geo.Pos[0],time1);

               if(eval_only)
               {
                  printf(
                     "EVAL_ONLY x0=%i Nv=%i lambda_Nv=%.17e\n",
                     lt+geo.lL[0]*geo.Pos[0],
                     params.Nv,
                     lambda[params.Nv-1]
                  );
                  fflush(stdout);
               }
            }

            if(eval_only)
               continue;

            if( params.write_ev )
            {
               // write eigen-vectors to files
               ret = snprintf(fname, qcd_MAX_STRING_LENGTH,"%s/%sn%i_evec_t%i.h5",params.evec_dir,params.nbase,nc,lt+geo.Pos[0]*geo.lL[0]);
               qcd_errorLoc(ret>=qcd_MAX_STRING_LENGTH,1,"main [mental.c]",
               "file name for eigenvectors too long");
               // copy eigenvectors to disteigvecs
               Ncp = geo.lV3*3*sizeof(DistEigvecsHdf5_ComplexDouble);
               for(j=0; j<params.Nv; j++)
               {
                  //TODO: write a direct writing routine without copy
                  memcpy(&(disteigvecs[j*geo.lV3*3]),&(v[j].D[0][0].re),Ncp);
               }
               int P3[3];
               P3[0] = geo.L[1]/geo.lL[1];
               P3[1] = geo.L[2]/geo.lL[2];
               P3[2] = geo.L[3]/geo.lL[3];
               ret = DistEigvecsHdf5Writer_CreateFile(&(params.evec_writer),
                        geo.comm3d, fname, params.nbase, nc,
                        &(geo.L[1]), P3, lt+geo.Pos[0]*geo.lL[0], params.Nv,
                        (DistEigvecsHdf5_Metadata *) &(params.wr_mdat),
                        (DistEigvecsHdf5_Programdata *) &(params.wr_pdat));
               qcd_errorLoc(ret != 0,1,"main [mental.c]",
               "DistEigvecsHdf5Writer_CreateFile failed.");
               ret = DistEigvecsHdf5Writer_WriteDistEigvecs(&(params.evec_writer),
                        &(geo.Pos[1]), disteigvecs);
               qcd_errorLoc(ret != 0,1,"main [mental.c]",
               "DistEigvecsHdf5Writer_WriteDistEigvecs failed.");
               ret = DistEigvecsHdf5Writer_CloseFile(&(params.evec_writer));
               qcd_errorLoc(ret != 0,1,"main [mental.c]",
               "DistEigvecsHdf5Writer_CloseFile failed.");
               // write eigen-values to files
               ret = DistEigvalsHdf5Writer_WriteDistEigvals(&(params.eval_writer),
                         lt+geo.Pos[0]*geo.lL[0], myid==geo.root3d[lt+geo.Pos[0]*geo.lL[0]], 
                         lambda);
               qcd_errorLoc(ret != 0,1,"main [mental.c]",
               "DistEigvalsHdf5Writer_WriteDistEigvals failed.");
            }
         }else
         {
            //load eigen-vectors from files
            ret = snprintf(fname, qcd_MAX_STRING_LENGTH,"%s/%sn%i_evec_t%i.h5",params.evec_dir,params.nbase,nc,lt+geo.Pos[0]*geo.lL[0]);
               qcd_errorLoc(ret>=qcd_MAX_STRING_LENGTH,1,"main [mental.c]",
               "file name for eigenvectors too long");
            int P3[3];
            P3[0] = geo.L[1]/geo.lL[1];
            P3[1] = geo.L[2]/geo.lL[2];
            P3[2] = geo.L[3]/geo.lL[3];
            char *tmpstr;
            int tmpnc;
            int tmpt;
            int tmpNv;
            DistEigvecsHdf5_Metadata tmpmdata;
            DistEigvecsHdf5_Programdata tmppdata;
            ret = DistEigvecsHdf5Reader_OpenFile(&(params.evec_reader), geo.comm3d,
                                                 fname, &tmpstr, &tmpnc, &(geo.L[1]),
                                                 P3, &tmpt, &tmpNv, 
                                                 &tmpmdata,
                                                 &tmppdata);
            qcd_errorLoc(ret != 0,1,"main [mental.c]",
              "DistEigvecsHdf5Reader_OpenFile failed.");
            qcd_errorLoc((tmpnc != nc) || (tmpt != lt+geo.Pos[0]*geo.lL[0]) || (tmpNv < params.Nv)
                      || (strcmp(tmpstr,params.nbase)!=0), 1, "main [mental.c]",
                      "Eigenvector file does not match run");
            
            qcd_errorLoc(cmp_meta((DistEigvecsHdf5_Metadata *) &(params.wr_mdat),&(tmpmdata)) != 0,1,"main [mental.c]",
              "meta-data in eigenvector files is inconsistent.");
            
            ret = DistEigvecsHdf5Reader_ReadDistEigvecsRange(&(params.evec_reader), &(geo.Pos[1]),
                                                             0, params.Nv, disteigvecs);
            qcd_errorLoc(ret != 0,1,"main [mental.c]",
              "DistEigvecsHdf5Reader_ReadDistEigvecsRange failed.");
            Ncp = geo.lV3*3*sizeof(DistEigvecsHdf5_ComplexDouble);
            for(j=0; j<params.Nv; j++)
            {
               //TODO: write a direct reading routine without copy
               memcpy(&(v[j].D[0][0].re), &(disteigvecs[j*geo.lV3*3]),Ncp);
            }
            
            ret =DistEigvecsHdf5Reader_CloseFile(&(params.evec_reader));
            qcd_errorLoc(ret != 0,1,"main [mental.c]",
              "DistEigvecsHdf5Reader_CloseFile failed.");
            if(myid==geo.root3d[lt+geo.lL[0]*geo.Pos[0]])
            {
               printf("Eigenvectors on x0=%i loaded from %s\n",lt+geo.lL[0]*geo.Pos[0],fname);
            }
         }
        
         if( params.hasDer && (uAPEd_ptr != uAPEl_ptr) )
         {
            qcd_copyGaugeField3dGaugeField(&u3d, uAPEd_ptr, lt+geo.Pos[0]*geo.lL[0]);
         }
        
         for(lp=0; lp<tabs.np; lp++)
         {
            memset(redelem, 0, params.Nv*params.Nv*params.nset*sizeof(qcd_complex_16));
            for(i=0; i<params.Nv; i++)
            {
               // multiply eigenvector v_i with exp(+ipx)
               vp.BPup2date = 0;
               vp.BMup2date = 0;
               for(lx=0; lx<geo.lL[1]; lx++)
               {
                  gx = lx + geo.lL[1]*geo.Pos[1];
                  phi = gx*tabs.p1[lp]*2.0*M_PI/((double) geo.L[1]);
                  eipx = (qcd_complex_16) { cos(phi), sin(phi)};
                  for(ly=0; ly<geo.lL[2]; ly++)
                  {
                     gy = ly + geo.lL[2]*geo.Pos[2];
                     phi = gy*tabs.p2[lp]*2.0*M_PI/((double) geo.L[2]);
                     eipy = (qcd_complex_16) { cos(phi), sin(phi)};
                     eipxy = qcd_CMUL(eipx,eipy);
                     for(lz=0; lz<geo.lL[3]; lz++)
                     {
                        gz = lz + geo.lL[3]*geo.Pos[3];
                        phi = gz*tabs.p3[lp]*2.0*M_PI/((double) geo.L[3]);
                        eipz = (qcd_complex_16) { cos(phi), sin(phi)};
                        eipxyz=qcd_CMUL(eipxy,eipz);
                        l = qcd_LEXIC0(lx,ly,lz,geo.lL);
                        vp.D[l][0] = qcd_CMUL(v[i].D[l][0],eipxyz);
                        vp.D[l][1] = qcd_CMUL(v[i].D[l][1],eipxyz);
                        vp.D[l][2] = qcd_CMUL(v[i].D[l][2],eipxyz);
                     }
                  }
               }//space loop
            
               // loop over sets, find those that need this momentum
               for(ls=0; ls<params.nset; ls++)
               {
                  if(tabs.pp[lp] <= params.p2max[ls])
                  {
                     //if(myid==0) printf("   set %i\n", ls);

                     //loop over eigenvectors j
                     for(j=0; j<params.Nv; j++)
                     {
                        //apply derivatives
                        // result in Dv
                        switch (params.nder[ls])
                        {
                           case 0:
                              qcd_copySpinorComponent3d(&Dv,&v[j]);                              
                              break;
                           case 1: 
                              qcd_DSymSpinorComponent3d(&Dv,&v[j],&u3d,params.dermu[ls][0]);
                              break;
                           case 2:
                              qcd_DSymSpinorComponent3d(&vtmp,&v[j],&u3d,params.dermu[ls][1]);
                              qcd_DSymSpinorComponent3d(&Dv,&vtmp,&u3d,params.dermu[ls][0]);
                              break;
                           case 3:
                              qcd_DSymSpinorComponent3d(&Dv,&v[j],&u3d,params.dermu[ls][2]);
                              qcd_DSymSpinorComponent3d(&vtmp,&Dv,&u3d,params.dermu[ls][1]);
                              qcd_DSymSpinorComponent3d(&Dv,&vtmp,&u3d,params.dermu[ls][0]);
                              break;
                           case 4:
                              qcd_DSymSpinorComponent3d(&vtmp,&v[j],&u3d,params.dermu[ls][3]);
                              qcd_DSymSpinorComponent3d(&Dv,&vtmp,&u3d,params.dermu[ls][2]);
                              qcd_DSymSpinorComponent3d(&vtmp,&Dv,&u3d,params.dermu[ls][1]);
                              qcd_DSymSpinorComponent3d(&Dv,&vtmp,&u3d,params.dermu[ls][0]);
                              break;
                        }
                        //contract v_i exp(-ipx) D v_j
                        redelem[ls*params.Nv*params.Nv + i*params.Nv + j] = qcd_mulAdjointSpinorComponent3dSpinorComponent3d(&vp,&Dv);
                     }//end eigenvector j loop
                  }
               }//end set-loop
            }//end eigenvector i loop
            // write Nv x Nv sized data chunk for all sets and one p
            for(ls=0; ls<params.nset; ls++)
            {
               if(tabs.pp[lp] <= params.p2max[ls])
               {
                  //muidx = get_muidx(ls, params.nder[ls]);
                  //pidx  = get_pidx(lp, params.nder[ls]);
                  tidx  = lt+geo.Pos[0]*geo.lL[0];
                  wr_mom[0] = tabs.p1[lp];
                  wr_mom[1] = tabs.p2[lp];
                  wr_mom[2] = tabs.p3[lp];
                  wr_mus = params.dermu[ls];
                  if(geo.root3d[lt+geo.Pos[0]*geo.lL[0]] == myid)
                     wrflag = 1;
                  else
                     wrflag = 0;
                  ret = DistRedelemsHdf5Writer_WriteDistRedelems(
                               &(params.writer[params.nder[ls]]), wr_mus, wr_mom, tidx, wrflag,
                               &(redelem[ls*params.Nv*params.Nv]));
                  qcd_errorLoc(ret!=0,1,"main [mental.c]","DistRedelemsHdf5Writer_WriteDistRedelems failed");
               }
            }
         }//end momenta loop
      }//end time-slice sloop
      closeDatFiles();
   }//end main loop (configs)

   ////////////////////////////////////// CLEAN UP AND EXIT ////////////////////
   qcd_destroySpinorComponent3d(&vp);
   qcd_destroySpinorComponent3d(&Dv);
   qcd_destroySpinorComponent3d(&vtmp);
   
   qcd_destroyGaugeField(&u1);
   if(u2.initialized)
      qcd_destroyGaugeField(&u2);
   if(u3.initialized)
      qcd_destroyGaugeField(&u3);
   qcd_destroyGaugeField3d(&u3d);
 //  qcd_destroyVector(&src);
 //  qcd_destroyVector(&sol);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
