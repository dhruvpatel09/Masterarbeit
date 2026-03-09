/* polya.c
 *
 * Computes polyakov loops
 *
 * where m are bare mass parameters 
 *
 * Tomasz Korzec 2022-2023
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <mpi.h>
#include <math.h>
#include <qcd.h>
#include <bdio.h>
#include <solver.h>

#define polya_RELEASE "polya 1.0"

static struct
{
   int append;  /* command line flag */
   char nbase[qcd_MAX_STRING_LENGTH];
   char outbase[qcd_MAX_STRING_LENGTH];
   char dat_file[qcd_MAX_STRING_LENGTH];
   char log_file[qcd_MAX_STRING_LENGTH];
   char end_file[qcd_MAX_STRING_LENGTH];
   char dat_save[qcd_MAX_STRING_LENGTH];
   char log_save[qcd_MAX_STRING_LENGTH];
   char log_dir[qcd_MAX_STRING_LENGTH];
   char cnfg_dir[qcd_MAX_STRING_LENGTH];
   char dat_dir[qcd_MAX_STRING_LENGTH];
   int nc_first;
   int nc_last;
   int nc_step;
} params;

static struct
{
   double *corr;
   char cnfg[qcd_MAX_STRING_LENGTH];
   int nc;
} data;

//MPI related
static int myid,numprocs;


static void error_root(int test,int no,char *name,char *format,...)
{
   va_list args;

   if ((myid==0)&&(test!=0))
   {
      fprintf(stderr,"\nError in %s:\n",name);
      va_start(args,format);
      vfprintf(stderr,format,args);
      va_end(args);
      fprintf(stderr,"\nProgram aborted\n\n");
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD,no);
   }
}


static char* find_val(char *str, char* tag)
{
   /* scan string for the occurence of tag, return pointer to substring
    * between tag and the next '\0'
    */
   int i;
   int len;
   int nm=0;
   char* valstr;

   len=strlen(tag);
   for (i=0; i<512-len; i++)
   {
      if (str[i]==tag[nm])
         nm++;
      else
         nm=0;
      if (nm==len)
      {
         valstr = str+i+1;
         return valstr;
      }
   }
   return NULL;
}

static int find_val_int(char *str, char* tag)
{
   /* scan string for the occurence of tag, return substring
    * between tag and the next '\0'
    */
   int i;
   int len;
   int nm=0;
   char* valstr;
   char* endptr;
   long val;

   len=strlen(tag);
   for (i=0; i<512-len; i++)
   {
      if (str[i]==tag[nm])
         nm++;
      else
         nm=0;
      if (nm==len)
      {
         valstr = str+i+1;
         val = strtol(valstr,&endptr,10);
         if (endptr[0]!='\0' || (int) val != val)
            return -1;
         else
            return (int) val;
      }
   }
   return -1;
}

static int find_val_dbl(char *str, char* tag, double *val)
{
   int i;
   int len;
   int nm=0;
   char* valstr;
   char* endptr;

   len=strlen(tag);
   for (i=0; i<512-len; i++)
   {
      if (str[i]==tag[nm])
         nm++;
      else
         nm=0;
      if (nm==len)
      {
         valstr = str+i+1;
         *val = strtod(valstr,&endptr);
         if (endptr[0]!='\0')
            return 1;
         else
            return 0;
      }
   }
   return 1;
}


static void write_file_head(qcd_geometry *geo)
{
   int i,j,lf;
   char str[512];
   int len;
   BDIO *fbdio;

   bdio_set_dflt_verbose(1);
   fbdio=bdio_open(params.dat_file,"a","Generic Correlator Format 1.0");
   error_root(fbdio==NULL,1,"write_file_head [polya.c]",
                    "Unable to open data file");

   /*write the lattice-info record*/
   error_root(bdio_start_record(BDIO_ASC_GENERIC, 0, fbdio)!=0,1,
              "write_file_head [polya.c]",
              "Could not start a lattice-info record");
   sprintf(str,"CREATOR=%s",polya_RELEASE);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"ENSEMBLE=%s",params.nbase);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"L0=%d",geo->L[0]);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"L1=%d",geo->L[1]);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"L2=%d",geo->L[2]);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"L3=%d",geo->L[3]);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"BC0=periodic");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"BC1=periodic");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"BC2=periodic");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");
   sprintf(str,"BC3=periodic");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
              "write_file_head [polya.c]",
              "Incorrect write count");

   error_root(bdio_start_record(BDIO_ASC_GENERIC, 1, fbdio)!=0,1,
               "write_file_head [polya.c]",
               "Could not start a correlator-info record");
   sprintf(str,"CORR_ID=0");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
               "write_file_head [polya.c]",
               "Incorrect write count");
   sprintf(str,"CORR_NAME=tr[P_mu]");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
               "write_file_head [polya.c]",
               "Incorrect write count");
   sprintf(str,"NDIM=1");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
               "write_file_head [polya.c]",
               "Incorrect write count");
   sprintf(str,"D0=4");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
               "write_file_head [polya.c]",
               "Incorrect write count");
   sprintf(str,"D0_NAME=mu");
   sprintf(str,"DATATYPE=complex");
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
               "write_file_head [polya.c]",
               "Incorrect write count");
   bdio_close(fbdio);
}

static void write_data()
{
   int iw,nw;
   int chunk;
   int i,len;
   double tmp;
   char str[512];
   BDIO *fbdio;

   fbdio=bdio_open(params.dat_file,"a","Generic Correlator Format 1.0");
   error_root(fbdio==NULL,1,"write_data [polya.c]",
              "Unable to open data file");

   /* write configuration description record */
   error_root(bdio_start_record(BDIO_ASC_GENERIC, 4, fbdio)!=0,1,
         "write_data [polya.c]",
         "Could not start a configuration description record");
   sprintf(str,"CNFG_ID=%i",data.nc);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
         "write_dat [polya.c]",
         "Incorrect write count");
   sprintf(str,"CNFG_NAME=%s",data.cnfg);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
         "write_dat [polya.c]",
         "Incorrect write count");
   sprintf(str,"NC=%i",data.nc);
   len = strlen(str)+1;
   error_root(bdio_write(str, len, fbdio)!=len,1,
         "write_dat [polya.c]",
         "Incorrect write count");

   /* write data records */
   error_root(bdio_start_record(BDIO_BIN_F64LE, 5, fbdio)!=0,1,
      "write_data [polya.c]",
      "Could not start a data record");
   tmp = (double) data.nc;
   iw = bdio_write_f64(&tmp,sizeof(double),fbdio);
   tmp = 0.0;
   iw+= bdio_write_f64(&tmp,sizeof(double),fbdio);
   nw = 2*sizeof(double);
   iw+=bdio_write_f64(data.corr, sizeof(double)*8,fbdio);
   nw+=sizeof(double)*8;
   error_root(iw!=nw,1,"write_data [polya.c]",
                       "Incorrect write count");
   bdio_close(fbdio);
}

void read_dirs(void)
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

      qcd_findOpenQCDSection("Configurations");
      qcd_readOpenQCDLine("first","%d",&(params.nc_first));
      qcd_readOpenQCDLine("last","%d",&(params.nc_last));
      qcd_readOpenQCDLine("step","%d",&(params.nc_step));

      error_root((params.nc_last<params.nc_first)||
                 (params.nc_step<1)||
                 (((params.nc_last-params.nc_first)%params.nc_step)!=0),1,
                 "read_dirs [polya.c]","Improper configuration range");
   }

   MPI_Bcast(params.nbase,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.outbase,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(params.log_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.cnfg_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(params.dat_dir,qcd_MAX_STRING_LENGTH,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&(params.nc_first),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.nc_last),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&(params.nc_step),1,MPI_INT,0,MPI_COMM_WORLD);
}

static void setup_files(void)
{
   error_root(qcd_nameSize("%s/%sn%d",params.cnfg_dir,params.nbase,params.nc_last)
              >=qcd_MAX_STRING_LENGTH, 1,"setup_files [polya.c]",
                             "cnfg_dir name is too long");

   error_root(!qcd_isDirWritable(params.dat_dir), 1, "setup_files [polya.c]",
                                              "dat directory is not writable");
   error_root(qcd_nameSize("%s/%s.polya.bdio~",params.dat_dir,params.outbase)
              >=qcd_MAX_STRING_LENGTH,
              1,"setup_files [polya.c]","dat_dir name is too long");

   error_root(!qcd_isDirWritable(params.log_dir),1, "setup_files [polya.c]",
                                             "log directory is not writable");
   error_root(qcd_nameSize("%s/%s.polya.log~",params.log_dir,params.outbase)
              >=qcd_MAX_STRING_LENGTH,
              1,"setup_files [polya.c]","log_dir name is too long");

   sprintf(params.log_file,"%s/%s.polya.log",params.log_dir,params.outbase);
   sprintf(params.end_file,"%s/%s.polya.end",params.log_dir,params.outbase);
   sprintf(params.dat_file,"%s/%s.polya.bdio",params.dat_dir,params.outbase);
   sprintf(params.log_save,"%s~",params.log_file);
   sprintf(params.dat_save,"%s~",params.dat_file);
}


static void read_lat_parms(int *L, int *P)
{
   if( myid==0 )
   {
      qcd_findOpenQCDSection("Lattice");
      printf("lattice section found\n"); fflush(stdout);
      qcd_readOpenQCDLine("L","%d %d %d %d",&(L[0]),&(L[1]),&(L[2]),&(L[3]));
      printf("L read\n"); fflush(stdout);
      qcd_readOpenQCDLine("P","%d %d %d %d",&(P[0]),&(P[1]),&(P[2]),&(P[3]));
      printf("P read\n");fflush(stdout);
   }

   MPI_Bcast(L,4,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(P,4,MPI_INT,0,MPI_COMM_WORLD);
}//end read_lat_parms


static void check_log()
{
   FILE *flog;
   flog=fopen(params.log_file,"r");
   error_root(flog==NULL,1,"check_log [polya.c]",
              "could not open the log file for reading");
   /*TODO parse log file, check whether parameters match with input-file */
   fclose(flog);
   flog = freopen(params.log_file,"a",stdout);
   error_root(flog==NULL,1,"check_log [polya.c]",
              "could not open log file for appending");
   flog = freopen(params.log_file,"a",stderr);
   error_root(flog==NULL,1,"check_log [polya.c]",
              "could not open log file for appending");
   remove("STARTUP_ERROR");
   return;
}

static void check_dat(int *L, int *P)
{
   BDIO *fbdio;
   int i, ir, lf;
   char str[1024];
   char *substr;
   int uinfo=-1;
   double d;
   int id,nc=-1,nc_last=-1;

   if(myid==0)
   {
      fbdio=bdio_open(params.dat_file,"r","Generic Correlator Format 1.0");
      error_root(fbdio==NULL,1,"check_dat [polya.c]",
               "Unable to open data file");

      /* check the lattice description record */
      i=bdio_seek_record(fbdio);
      error_root(i==EOF,1,"check_dat [polya.c]",
                     "No records in old data file found");
      error_root(bdio_get_ruinfo(fbdio)!=0,1,"check_dat [polya.c]",
                     "Data file does not start with a lattice record");
      error_root(bdio_get_rfmt(fbdio)!=BDIO_ASC_GENERIC,1,
                     "check_dat [polya.c]","Wrong record type");
      error_root(bdio_get_rlen(fbdio)>1024,1,"check_dat [dSdm.c]",
                     "Lattice record too long");
      ir=bdio_read(str, bdio_get_rlen(fbdio), fbdio);
      error_root(ir!=bdio_get_rlen(fbdio),1,"check_dat [polya.c]",
                     "Could not read from data file");
      substr=find_val(str,"CREATOR=");
      error_root(strcmp(substr,polya_RELEASE)!=0,1,"check_dat [polya.c]",
                     "Not a valid %s file",polya_RELEASE);
      substr=find_val(str,"ENSEMBLE=");
      error_root(strcmp(substr,params.nbase)!=0,1,"check_dat [polya.c]",
                     "Cannot append: ensemble mismatch");
      i=find_val_int(str,"L0=");
      error_root(i!=L[0],1,"check_dat [polya.c]",
                     "Cannot append: lattice size mismatch");
      i=find_val_int(str,"L1=");
      error_root(i!=L[1],1,"check_dat [polya.c]",
                     "Cannot append: lattice size mismatch");
      i=find_val_int(str,"L2=");
      error_root(i!=L[2],1,"check_dat [polya.c]",
                     "Cannot append: lattice size mismatch");
      i=find_val_int(str,"L3=");
      error_root(i!=L[3],1,"check_dat [polya.c]",
                     "Cannot append: lattice size mismatch");
      substr=find_val(str,"BC0=");
      error_root(strcmp(substr,"periodic")!=0,1,"check_dat [polya.c]",
                     "Cannot append: b.c. mismatch");
      substr=find_val(str,"BC1=");
      error_root(strcmp(substr,"periodic")!=0,1,"check_dat [polya.c]",
                     "Cannot append: b.c. mismatch");
      substr=find_val(str,"BC2=");
      error_root(strcmp(substr,"periodic")!=0,1,"check_dat [polya.c]",
                     "Cannot append: b.c. mismatch");
      substr=find_val(str,"BC3=");
      error_root(strcmp(substr,"periodic")!=0,1,"check_dat [polya.c]",
                     "Cannot append: b.c. mismatch");
      i=bdio_seek_record(fbdio);
      error_root(i!=0,1,"check_dat [polya.c]",
                     "Expected correlator description record not found");
      error_root(bdio_get_ruinfo(fbdio)!=1,1,"check_dat [polya.c]",
                     "Expected correlator description record not found");
      error_root(bdio_get_rfmt(fbdio)!=BDIO_ASC_GENERIC,1,"check_dat [polya.c]",
                     "Correlator description record has wrong format");
      error_root(bdio_get_rlen(fbdio)>1024,1,"check_dat [polya.c]",
                        "Correlator record too long");
      ir=bdio_read(str, bdio_get_rlen(fbdio), fbdio);
      error_root(ir!=bdio_get_rlen(fbdio),1,"check_dat [polya.c]",
                        "Could not read from data file");
      i=find_val_int(str,"CORR_ID=");
      error_root(i!=0,1,"check_dat [polya.c]",
                        "Unexpected CORR_ID");
      substr=find_val(str,"CORR_NAME=");
      error_root(strcmp(substr,"tr[P_mu]")!=0,1,"check_dat [polya.c]",
                        "Unexpected CORR_NAME");
      i=find_val_int(str,"NDIM=");
      error_root(i!=1,1,"check_dat [polya.c]",
                        "Unexpected NDIM");
      i=find_val_int(str,"D0=");
      error_root(i!=4,1,"check_dat [polya.c]",
                        "Unexpected D0");
      substr=find_val(str,"DATATYPE=");
      error_root(strcmp(substr,"complex")!=0,1,"check_dat [polya.c]",
                        "Unexpected DATATYPE");
      //skip till first cnfg-description record
      while((bdio_get_ruinfo(fbdio)!=4) && 
         (fbdio->state==BDIO_R_STATE || fbdio->state==BDIO_H_STATE))
      {
         i=bdio_seek_record(fbdio);
      }

      //read measured correlators, check spacings
      while((bdio_get_ruinfo(fbdio)==4) && 
         (fbdio->state==BDIO_R_STATE || fbdio->state==BDIO_H_STATE))
      {
         //read cnfg description record
         error_root(bdio_get_rlen(fbdio)>512,1,"check_data [polya.c]",
                  "Too long configuration description record encountered");
         i=bdio_read(str,bdio_get_rlen(fbdio),fbdio);
         error_root(bdio_get_rlen(fbdio)!=i,1,"check_data [polya.c]",
                  "Could not read from data file");
         id=find_val_int(str,"CNFG_ID=");
         nc=find_val_int(str,"NC=");
         if(nc_last >0)
         {
            error_root((nc-nc_last)!=params.nc_step,1,"check_data [polya.c]",
                     "Wrong cnfg spacing in data file");
         }
         nc_last=nc;
         i=bdio_seek_record(fbdio);
         error_root(i==EOF,1,"check_data [polya.c]",
                           "Missing data record");
         error_root(bdio_get_ruinfo(fbdio)!=5,1,"check_data [polya.c]",
                           "Not a valid %s file", polya_RELEASE);
         error_root((8+2)*8!=bdio_get_rlen(fbdio),1,
                     "check_data [polya.c]",
                     "Incomplete or oversized data record encountered");
         i=bdio_read_f64(&d, 8, fbdio);
         error_root(i!=8,1,"check_data [polya.c]",
                           "Could not read from data file");
         error_root(id!=(int) d,1,"check_data [polya.c]",
                           "Not a valid %s file", polya_RELEASE);
         i=bdio_read_f64(&d, 8, fbdio);
         error_root(i!=8,1,"check_data [polya.c]",
                           "Could not read from data file");
         error_root(lf!=(int) d,1,"check_data [polya.c]",
                           "Not a valid %s file", polya_RELEASE);
         i=bdio_seek_record(fbdio);
      }
      if(nc>0)
         error_root(params.nc_first!=(nc+params.nc_step),1,"check_data [polya.c]",
                  "Cannot append: nc_first does not continue the data file");

      bdio_close(fbdio);
   }
}


static void read_infile(int argc, char *argv[], int *L, int *P)
{
   int lf;
   // command line option positions/flags
   int ifile;
   int append;

   // file pointers for log and input files
   FILE *flog, *fin;

   if(myid==0)
   {
      flog=freopen("STARTUP_ERROR","a",stdout);
      flog=freopen("STARTUP_ERROR","a",stderr);
 
      ifile=qcd_findOpt(argc,argv,"-i");

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [polya.c]",
                 "Syntax: polya -i <input file> [-a]");

      error_root(qcd_isBigEndian()==-1,1,"read_infile [polya.c]",
                 "Machine has unknown endianness");

      append=qcd_findOpt(argc,argv,"-a");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [polya.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);
   params.append = append;


   read_dirs();
   if(myid==0){ printf("dirs read\n"); fflush(stdout);}
   setup_files();
   if(myid==0){ printf("files set up\n"); fflush(stdout);}
   read_lat_parms(L,P);
   if(myid==0){ printf("lat parms read\n"); fflush(stdout);}


   if (myid==0)
   {
      fclose(fin);
   }


   if (append)
   {
      // read old bdio, read input, check for consistency
      if(myid==0)
      {
         check_log();
         check_dat(L,P);
         printf("\n\n%s continuation run\n\n",polya_RELEASE);
      }
   }else
   {
      if(myid==0)
      {
         flog=NULL;
         flog=fopen(params.log_file,"r");
         error_root(flog!=NULL,1,"read_infile [polya.c]",
              "log file already exists");
         flog = freopen(params.log_file,"a",stdout);
         error_root(flog==NULL,1,"read_infile [polya.c]",
              "could not open log file for writing");
         flog = freopen(params.log_file,"a",stderr);
         error_root(flog==NULL,1,"read_infile [polya.c]",
              "could not open log file for writing");
         remove("STARTUP_ERROR");

         printf("%s log file\n\n",polya_RELEASE);
         printf("Global lattice: %i x %i x %i x %i\n",L[0],L[1],L[2],L[3]);
         printf("Local lattice:  %i x %i x %i x %i\n",L[0]/P[0],L[1]/P[1],L[2]/P[2],L[3]/P[3]);
         printf("Process grid:   %i x %i x %i x %i\n",P[0],P[1],P[2],P[3]);
         printf("\nPeriodic boundary conditions\n\n");

         printf("ensemble: %s\n",params.nbase);
         printf("output:   %s\n",params.outbase);
         printf("\nProcessing cnfgs %i to %i in steps of %i\n"
                          ,params.nc_first,params.nc_last,params.nc_step);
      }
   }
   return;
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


void compute(double *corr, qcd_gaugeField *u)
{
   // main function: computes average polyakov-loops in all directions
   qcd_geometry *geo = u->geo;
   qcd_complex_16 C[3][3],Cp[3][3],CC[3][3],CP[3][3],trC,trCtmp;
   MPI_Status stat;
   int mu,x0,x1,x2,x3,ip,np;
   long l;

   // direction 0
   mu=0;
   np = geo->L[mu]/geo->lL[mu];
   corr[2*mu]   = 0.0;
   corr[2*mu+1] = 0.0;
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      qcd_unit3x3(C);
      for(x0=0; x0<geo->lL[0]; x0++)
      {
	 l = qcd_LEXIC(x0,x1,x2,x3,geo->lL);
         qcd_mul3x3(CC,C,u->D[l][mu]);
	 qcd_copy3x3(C,CC);
      }
      // multiply local results
      for(ip=1; ip<np; ip++)
      {
	 MPI_Sendrecv_replace(CC, 18, MPI_DOUBLE, geo->Pplus[mu], 0, geo->Pminus[mu], 0, MPI_COMM_WORLD, &stat);
	 qcd_mul3x3(CP,CC,C);
	 qcd_copy3x3(CC,CP);
      }
      trC = qcd_trace3x3(CC);
      MPI_Allreduce(&(trC.re), &(trCtmp.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      corr[2*mu]   += trCtmp.re/np;
      corr[2*mu+1] += trCtmp.im/np;
   }
   corr[2*mu]   /= (geo->V/geo->L[mu]);
   corr[2*mu+1] /= (geo->V/geo->L[mu]);

   // direction 1
   mu=1;
   np = geo->L[mu]/geo->lL[mu];
   corr[2*mu]   = 0.0;
   corr[2*mu+1] = 0.0;
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x2=0; x2<geo->lL[2]; x2++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      qcd_unit3x3(C);
      for(x1=0; x1<geo->lL[1]; x1++)
      {
         l = qcd_LEXIC(x0,x1,x2,x3,geo->lL);
         qcd_mul3x3(CC,C,u->D[l][mu]);
         qcd_copy3x3(C,CC);
      }
      // multiply local results
      for(ip=1; ip<np; ip++)
      {
         MPI_Sendrecv_replace(CC, 18, MPI_DOUBLE, geo->Pplus[mu], 0, geo->Pminus[mu], 0, MPI_COMM_WORLD, &stat);
         qcd_mul3x3(CP,CC,C);
         qcd_copy3x3(CC,CP);
      }
      trC = qcd_trace3x3(CC);
      MPI_Allreduce(&(trC.re), &(trCtmp.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      corr[2*mu]   += trCtmp.re/np;
      corr[2*mu+1] += trCtmp.im/np;
   }
   corr[2*mu]   /= (geo->V/geo->L[mu]);
   corr[2*mu+1] /= (geo->V/geo->L[mu]);
   
   // direction 2
   mu=2;
   np = geo->L[mu]/geo->lL[mu];
   corr[2*mu]   = 0.0;
   corr[2*mu+1] = 0.0;
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x3=0; x3<geo->lL[3]; x3++)
   {
      qcd_unit3x3(C);
      for(x2=0; x2<geo->lL[2]; x2++)
      {
         l = qcd_LEXIC(x0,x1,x2,x3,geo->lL);
         qcd_mul3x3(CC,C,u->D[l][mu]);
         qcd_copy3x3(C,CC);
      }
      // multiply local results
      for(ip=1; ip<np; ip++)
      {
         MPI_Sendrecv_replace(CC, 18, MPI_DOUBLE, geo->Pplus[mu], 0, geo->Pminus[mu], 0, MPI_COMM_WORLD, &stat);
         qcd_mul3x3(CP,CC,C);
         qcd_copy3x3(CC,CP);
      }
      trC = qcd_trace3x3(CC);
      MPI_Allreduce(&(trC.re), &(trCtmp.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      corr[2*mu]   += trCtmp.re/np;
      corr[2*mu+1] += trCtmp.im/np;
   }
   corr[2*mu]   /= (geo->V/geo->L[mu]);
   corr[2*mu+1] /= (geo->V/geo->L[mu]);
   
   // direction 3
   mu=3;
   np = geo->L[mu]/geo->lL[mu];
   corr[2*mu]   = 0.0;
   corr[2*mu+1] = 0.0;
   for(x0=0; x0<geo->lL[0]; x0++)
   for(x1=0; x1<geo->lL[1]; x1++)
   for(x2=0; x2<geo->lL[2]; x2++)
   {
      qcd_unit3x3(C);
      for(x3=0; x3<geo->lL[3]; x3++)
      {
         l = qcd_LEXIC(x0,x1,x2,x3,geo->lL);
         qcd_mul3x3(CC,C,u->D[l][mu]);
         qcd_copy3x3(C,CC);
      }
      // multiply local results
      for(ip=1; ip<np; ip++)
      {
         MPI_Sendrecv_replace(CC, 18, MPI_DOUBLE, geo->Pplus[mu], 0, geo->Pminus[mu], 0, MPI_COMM_WORLD, &stat);
         qcd_mul3x3(CP,CC,C);
         qcd_copy3x3(CC,CP);
      }
      trC = qcd_trace3x3(CC);
      MPI_Allreduce(&(trC.re), &(trCtmp.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      corr[2*mu]   += trCtmp.re/np;
      corr[2*mu+1] += trCtmp.im/np;
   }
   corr[2*mu]   /= (geo->V/geo->L[mu]);
   corr[2*mu+1] /= (geo->V/geo->L[mu]);
}



int main(int argc,char* argv[])
{
   double   plaq;                              //average plaquette

   //loop variables
   int j,nc;

   //gauge field
   qcd_gaugeField u;

   //geometry parameters
   qcd_geometry geo;
   int P[4];                                   //number of processes in x0,x1,x2,x3 directions
   int L[4];                                   //global lattice size
   double theta[4]={0.,0.,0.,0.};              //boundary conditions, gauge-fields must be periodic



   //timing
   double time1, time2;

   /////////////////////////////////////////////////////////////////////////////
   // set up MPI                                                              //
   /////////////////////////////////////////////////////////////////////////////

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs); // number of processes
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);     // each process gets its ID

   read_infile(argc, argv, L, P);
   if(qcd_initGeometry(&geo,L,P, theta, myid, numprocs)) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   //if(myid==0) printf(" Local lattice: %i x %i x %i x %i\n",geo.lL[0],geo.lL[1],geo.lL[2],geo.lL[3]);
   if(qcd_initEO(&geo)) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);//needed for cern format
   //if(myid==0){printf(" E/O initialized\n"); fflush(stdout);}


   /////////////////////////////////////////////////////////////////////////////
   // initialize gauge field                                                  //
   /////////////////////////////////////////////////////////////////////////////

   j = qcd_initGaugeField(&u,&geo);
   if(j>0)
   {
      if(myid==0) fprintf(stderr,"Error, not enough memory\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
   
   // initialize the data file
   if((myid==0) && (!params.append))
      write_file_head(&geo);
   //if(myid==0){printf("file header written\n"); fflush(stdout);}

   data.corr = malloc(4*2*sizeof(double));
   if(data.corr==NULL)
   {
      fprintf(stderr,"process %i: out of memory for correlators.\n",myid);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }

   for(nc=params.nc_first; nc<=params.nc_last; nc+=params.nc_step)
   {
      sprintf(data.cnfg,"%s/%sn%i",params.cnfg_dir,params.nbase,nc);
      if(myid==0) printf("\nProcessing gauge field %s\n",data.cnfg);
      data.nc = nc;
      time1 = MPI_Wtime();
      //read in the gauge field in "CERN" format
      if(qcd_getGaugeField(data.cnfg,qcd_GF_OPENQCD,&u))
      {
         fprintf(stderr,"process %i: Error reading gauge field!\n",myid);
         fflush(stdout);
         MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      time2 = MPI_Wtime();
      if(myid==0) {printf("gauge-field loaded in %fsek\n",time2-time1); fflush(stdout);}

      plaq = qcd_calculatePlaquette(&u);
      if(myid==0) {printf("plaquette (torus def) = %e\n",plaq); fflush(stdout);}


      compute(data.corr, &u);
      if(myid==0)
      {
         printf("ReTr P_0=%f\n",data.corr[0]);
	 fflush(stdout);
      }

      if(myid==0)
         write_data();

      if(myid==0)
      {
         error_root(qcd_copyFile(params.log_file, params.log_save),1,
                       "main [polya.c]","Could not copy the log file");
         error_root(qcd_copyFile(params.dat_file, params.dat_save),1,
                    "main [polya.c]","Could not copy the dat file");
      }

      if( check_endflag() )
         break;
   }

   ////////////////////////////////////// CLEAN UP AND EXIT ////////////////////
   free(data.corr);
   qcd_destroyGaugeField(&u);
   qcd_destroyGeometry(&geo);
   MPI_Finalize();
   return(EXIT_SUCCESS);
}//end main
