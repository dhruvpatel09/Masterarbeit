#ifndef QX_AUDITED_VALUES_H
#define QX_AUDITED_VALUES_H
#include <stdio.h>
#include <string.h>
#include <math.h>
/* The audit prints 17 significant decimal digits, preserving each double. */
static int read_audited_values(const char *path,double values[192][100])
{
    FILE *in=fopen(path,"r");char line[128];int bad=!in;
    if(!bad && (!fgets(line,sizeof(line),in) || strcmp(line,"t,i,lambda\n")))bad=1;
    for(int t=0;t<192 && !bad;++t)for(int i=0;i<100;++i) {
        int a,b;double x;
        if(fscanf(in,"%d,%d,%lf",&a,&b,&x)!=3 || a!=t || b!=i || !isfinite(x) || x<=0 || x>=0.15 ||
           (i && x+1e-12<values[t][i-1])) {bad=1;break;}
        values[t][i]=x;
    }
    if(in) {char extra;if(!bad && fscanf(in," %c",&extra)!=EOF)bad=1;fclose(in);}
    return bad;
}
#endif
