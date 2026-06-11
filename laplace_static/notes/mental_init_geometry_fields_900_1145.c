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

   MPI_Init(&argc, &argv);
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
