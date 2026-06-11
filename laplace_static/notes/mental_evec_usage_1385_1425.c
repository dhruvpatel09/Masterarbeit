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
