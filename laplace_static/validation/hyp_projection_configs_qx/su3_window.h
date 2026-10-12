#ifndef PROFILE_SU3_WINDOW_H
#define PROFILE_SU3_WINDOW_H

/* Diagnostics only: no projection or changes to the transported field. */
static void profile_link_defects(const qcd_complex_16 u[3][3],
                                 double *unitarity, double *determinant)
{
    double square=0.0;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) {
        double re=-(i==j),im=0.0;
        for(int k=0;k<3;k++) {
            re+=u[k][i].re*u[k][j].re+u[k][i].im*u[k][j].im;
            im+=u[k][i].re*u[k][j].im-u[k][i].im*u[k][j].re;
        }
        square+=re*re+im*im;
    }
    const int perm[6][3]={{0,1,2},{1,2,0},{2,0,1},
                          {0,2,1},{1,0,2},{2,1,0}};
    double dr=0.0,di=0.0;
    for(int p=0;p<6;p++) {
        double re=1.0,im=0.0;
        for(int row=0;row<3;row++) {
            qcd_complex_16 a=u[row][perm[p][row]];
            double next=re*a.re-im*a.im;
            im=re*a.im+im*a.re;re=next;
        }
        double sign=p<3?1.0:-1.0;dr+=sign*re;di+=sign*im;
    }
    *unitarity=sqrt(square);*determinant=hypot(dr-1.0,di);
}

#ifndef PROFILE_SU3_KERNEL_ONLY
static void profile_window_su3(qcd_gaugeField *u,int cfg,int piter,
                               int t0,int t1,int rank)
{
    qcd_geometry *g=u->geo;
    double window_max[2]={0,0},window_sq[2]={0,0};
    unsigned long long window_count[3]={0,0,0};
    for(int t=t0;t<t1;t++) {
        double maximum[2]={0,0},squares[2]={0,0};
        unsigned long long counts[3]={0,0,0};
        int bad=0,any;
        for(int x=0;x<g->lL[1];x++)
        for(int y=0;y<g->lL[2];y++)
        for(int z=0;z<g->lL[3];z++) {
            long v=qcd_LEXIC(t,x,y,z,g->lL);
            double a,b;profile_link_defects(u->D[v][0],&a,&b);
            if(!isfinite(a)||!isfinite(b))bad=1;
            maximum[0]=fmax(maximum[0],a);maximum[1]=fmax(maximum[1],b);
            squares[0]+=a*a;squares[1]+=b*b;
            counts[0]++;counts[1]+=(a>1e-2);counts[2]+=(b>1e-2);
        }
        MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        if(any) {if(!rank)fprintf(stderr,"FAIL non-finite SU(3) diagnostic\n");MPI_Abort(MPI_COMM_WORLD,1);}
        double global_max[2],global_sq[2];unsigned long long global_count[3];
        MPI_Allreduce(maximum,global_max,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(squares,global_sq,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(counts,global_count,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
        if(global_count[0]!=64ULL*64*64)MPI_Abort(MPI_COMM_WORLD,1);
        if(!rank)printf("SU3_TIME cfg=%d piter=%d t=%d count=%llu unit_max=%.17e unit_rms=%.17e det_max=%.17e det_rms=%.17e unit_above_1em2=%llu det_above_1em2=%llu\n",
            cfg,piter,t,global_count[0],global_max[0],sqrt(global_sq[0]/global_count[0]),
            global_max[1],sqrt(global_sq[1]/global_count[0]),global_count[1],global_count[2]);
        for(int k=0;k<2;k++) {window_max[k]=fmax(window_max[k],global_max[k]);window_sq[k]+=global_sq[k];}
        for(int k=0;k<3;k++)window_count[k]+=global_count[k];
    }
    if(!rank)printf("SU3_WINDOW cfg=%d piter=%d t0=%d t1=%d count=%llu unit_max=%.17e unit_rms=%.17e det_max=%.17e det_rms=%.17e unit_above_1em2=%llu det_above_1em2=%llu\n",
        cfg,piter,t0,t1,window_count[0],window_max[0],sqrt(window_sq[0]/window_count[0]),
        window_max[1],sqrt(window_sq[1]/window_count[0]),window_count[1],window_count[2]);
    /* piter4 is a reported finite-iteration prescription, without a physics
       acceptance threshold. Check the converged piter40 control separately. */
    if(piter==40 && (window_max[0]>1e-12 || window_max[1]>1e-12)) {
        if(!rank)fprintf(stderr,"FAIL piter40 control SU(3) tolerance\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
}
#endif
#endif
