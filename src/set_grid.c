#include "pdisk.h"

void set_grid(void) {
    MALLOC_SAFE(( rc =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((rmin =  (double *)malloc(sizeof(double)*(NR+1))));
    MALLOC_SAFE((lam =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((mdot = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((dr  =  (double *)malloc(sizeof(double)*NR)));
    dlr = log(params.ro/params.ri)/((double)NR);
    MALLOC_SAFE((lrc = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((lrmin = (double *)malloc(sizeof(double)*(NR+1))));

    MALLOC_SAFE(( taumin =  (double *)malloc(sizeof(double)*(NR+1))));
    MALLOC_SAFE(( tauc =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE(( weights =  (double *)malloc(sizeof(double)*NR)));

    int set_flag = TRUE;
    int i;

    for(i=0;i<NR+1;i++) {
        lrmin[i] = log(params.ri) + i*dlr;
        rmin[i] = exp(lrmin[i]);
        taumin[i] = 0;
    }
    for (i=0;i<NR;i++) {
        rc[i] = .5*(rmin[i] + rmin[i+1]);
        lrc[i] = log(rc[i]);
        dr[i] = rmin[i+1]-rmin[i];
        lam[i] = 0;
        tauc[i] = 0;
        weights[i] = 1;
    }
    weights[0] = 0.5; weights[NR-1] = 0.5;
   
    printf("Grid spacing is dlr = %.3e\n",dlr);    

    if (params.forced_torque) {
        read_torque_file(&fld,params.torque_file);
    }
    set_mdot(params.planet_torque);
    return;

}

void free_grid(void) {
    free(rc); free(rmin); free(lam); free(dr);
    free(lrc); free(lrmin);
    free(mdot); free(taumin);
    free(tauc);
    free(weights);
     return;
}

