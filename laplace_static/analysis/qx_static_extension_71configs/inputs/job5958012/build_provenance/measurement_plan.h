/* Generated from measurement_plan.py. */
static const int mode_cutoffs[]={40,60,80,100};
static int selected_radius(int r) {
    static const int radii[]={1,2,4,8,12,16,20,23,24,25,28,32};
    for(unsigned int i=0;i<sizeof(radii)/sizeof(radii[0]);i++)if(r==radii[i])return 1;
    return 0;
}
