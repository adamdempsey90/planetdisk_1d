#include "pdisk.h"

void set_grid(void) {
    MALLOC_SAFE(( rc =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((rmin =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((lam =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((mdot = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((dr  =  (double *)malloc(sizeof(double)*NR)));
    dlr = log(params.ro/params.ri)/((double)(NR-1));
    MALLOC_SAFE((lrc = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE((lrmin = (double *)malloc(sizeof(double)*NR)));

    MALLOC_SAFE(( taumin =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE(( tauc =  (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE(( weights =  (double *)malloc(sizeof(double)*NR)));

    int set_flag = TRUE;
    int i;

    for (i=0;i<NR;i++) {
        lrc[i] = log(params.ri) + i*dlr;
        rc[i] = exp(lrc[i]);
        dr[i] = rc[i]*dlr;
        lam[i] = 0;
        taumin[i] = 0;
        tauc[i] = 0;
        weights[i] = 1;
    }
    weights[0] = 0.5; weights[NR-1] = 0.5;
    lrmin[0] = (log(params.ri) - dlr + lrc[0])*.5;
    rmin[0] = exp(lrmin[0]);//.5*(rc[0] + exp(log(params.ri)-dlr));
    for(i=1;i<NR;i++) {
        lrmin[i] = .5*(lrc[i]+lrc[i-1]);
        rmin[i] = exp(lrmin[i]);
    }
   
    printf("Grid spacing is dlr = %.3e\n",dlr);    

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

