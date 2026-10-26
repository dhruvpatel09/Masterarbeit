#ifndef QX_RECOVERY_SELECTION_H
#define QX_RECOVERY_SELECTION_H
static const int selected_times[]={0,1,32,64,95,96,97,98,99,100,101,128,160,191};
static int qx_selected(int t)
{
    for(unsigned k=0;k<sizeof(selected_times)/sizeof(selected_times[0]);++k)
        if(t==selected_times[k])return 1;
    return 0;
}
#endif
