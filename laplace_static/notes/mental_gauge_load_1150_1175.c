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
      
