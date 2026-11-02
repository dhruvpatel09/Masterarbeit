static void rank_neighbours(const qcd_geometry *g,int minus[3],int plus[3])
{
    int ranks; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    if(ranks!=64) MPI_Abort(MPI_COMM_WORLD,1);
    int local[4]={g->Pos[0],g->Pos[1],g->Pos[2],g->Pos[3]},positions[64][4];
    MPI_Allgather(local,4,MPI_INT,positions,4,MPI_INT,MPI_COMM_WORLD);
    for(int a=0;a<3;++a) {
        minus[a]=plus[a]=-1;
        for(int r=0;r<ranks;++r) {
            int isminus=1,isplus=1;
            for(int k=0;k<4;++k) {
                int m=local[k],p=local[k];
                if(k==a+1) { m=(m+3)%4; p=(p+1)%4; }
                if(positions[r][k]!=m) isminus=0;
                if(positions[r][k]!=p) isplus=0;
            }
            if(isminus) { if(minus[a]!=-1) MPI_Abort(MPI_COMM_WORLD,1); minus[a]=r; }
            if(isplus) { if(plus[a]!=-1) MPI_Abort(MPI_COMM_WORLD,1); plus[a]=r; }
        }
        if(minus[a]<0 || plus[a]<0) MPI_Abort(MPI_COMM_WORLD,1);
    }
}

static void shift_endpoint(const endpoint_packet *in,endpoint_packet *out,
                           endpoint_packet *send,endpoint_packet *recv,
                           const qcd_geometry *g,int axis,int minus,int plus)
{
    int b=(axis+1)%3,c=(axis+2)%3;
    int nb=g->lL[b+1],nc=g->lL[c+1],na=g->lL[axis+1];
    for(int j=0;j<nb;++j) for(int k=0;k<nc;++k) {
        int p[3]={0}; p[b]=j; p[c]=k;
        send[j*nc+k]=in[qcd_LEXIC0(p[0],p[1],p[2],g->lL)];
    }
    size_t bytes=(size_t)nb*nc*sizeof(*send);
    if(bytes>INT_MAX) MPI_Abort(MPI_COMM_WORLD,1);
    int rc=MPI_Sendrecv(send,(int)bytes,MPI_BYTE,minus,170+axis,
                        recv,(int)bytes,MPI_BYTE,plus,170+axis,
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    if(rc) MPI_Abort(MPI_COMM_WORLD,1);
    for(int x=0;x<g->lL[1];++x)
    for(int y=0;y<g->lL[2];++y)
    for(int z=0;z<g->lL[3];++z) {
        int p[3]={x,y,z}; long v=qcd_LEXIC0(x,y,z,g->lL);
        if(p[axis]+1==na) out[v]=recv[p[b]*nc+p[c]];
        else { ++p[axis]; out[v]=in[qcd_LEXIC0(p[0],p[1],p[2],g->lL)]; }
    }
}

static void verify_shift(const endpoint_packet *base,const endpoint_packet *shifted,
                         long volume,int axis,int r)
{
    int bad=0,any;
    for(long v=0;v<volume;++v) for(int k=0;k<3;++k)
        if(shifted[v].xyz[k]!=(base[v].xyz[k]+(k==axis?r:0))%NS) bad=1;
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(any) MPI_Abort(MPI_COMM_WORLD,1);
}

static void kahan_add(double value,double *sum,double *correction)
{
    double y=value-*correction, t=*sum+y;
    *correction=(t-*sum)-y;
    *sum=t;
}
