/* Compare original and observation-only solver objects on a known spectrum. */
#include <qcd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char **argv)
{
    if(argc!=5)return 1;
    qcd_geometry g;memset(&g,0,sizeof(g));g.lV3=16;g.lL[0]=1;
    qcd_gaugeField3d gauge={&g};qcd_spinorComponent3d v[3];double eigenvalues[3]={NAN,NAN,NAN};
    for(int i=0;i<3;i++)if(qcd_initSpinorComponent3d(&v[i],&g))return 1;
    qcd_lanczosParams p;qcd_defaultLanczosParams(&p);p.res=1e-9;p.deg=10;p.nLanczosFactor=atoi(argv[4]);
    p.upper=strtod(argv[3],NULL);p.maxRestarts=atoi(argv[2]);
    qcd_lanczosStatus status;memset(&status,0,sizeof(status));
    int code=qcd_lanczosCheby(eigenvalues,v,&status,&gauge,3,&p);
    FILE *f=fopen(argv[1],"wb");if(!f)return 1;
    fwrite(&code,sizeof(code),1,f);fwrite(&status,sizeof(status),1,f);
    fwrite(eigenvalues,sizeof(*eigenvalues),3,f);
    for(int i=0;i<3;i++)fwrite(v[i].D,sizeof(*v[i].D),g.lV3,f);
    fclose(f);
    printf("FIXTURE code=%d restarts=%d residual=%.17e lambda0=%.17e lambda2=%.17e\n",code,status.restarts,status.maxRes,eigenvalues[0],eigenvalues[2]);
    for(int i=0;i<3;i++)qcd_destroySpinorComponent3d(&v[i]);
    return 0;
}
