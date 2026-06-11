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
            }
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
