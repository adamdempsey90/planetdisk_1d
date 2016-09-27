#include "pdisk.h"
void allocate_field( Field *tmpfld) {
     int i;
     MALLOC_SAFE( (tmpfld->times = (double *)malloc(sizeof(double)*params.nt)));
     MALLOC_SAFE( (tmpfld->avals =  (double *)malloc(sizeof(double)*params.nt)));
     MALLOC_SAFE( (tmpfld->vs = (double *)malloc(sizeof(double)*params.nt)));
     MALLOC_SAFE( (tmpfld->sol = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->sol_mdot =  (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->torque = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->lami =  (double *)malloc(sizeof(double)*NR)));
     MALLOC_SAFE( (tmpfld->mdoti = (double *)malloc(sizeof(double)*NR)));

     MALLOC_SAFE( (tmpfld->mdot_ss =  (double *)malloc(sizeof(double)*params.nt)));
     MALLOC_SAFE( (tmpfld->vs_ss = (double *)malloc(sizeof(double)*params.nt)));
     MALLOC_SAFE( (tmpfld->sol_ss = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->ivals_ss = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->kvals_ss = (double *)malloc(sizeof(double)*NR*params.nt)));

     MALLOC_SAFE( (tmpfld->dTr = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->dep_func = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->mdL = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->mdR = (double *)malloc(sizeof(double)*NR*params.nt)));
    MALLOC_SAFE( (tmpfld->lamp = (double *)malloc(sizeof(double)*NR*params.nt)));
    MALLOC_SAFE( (tmpfld->lam0 = (double *)malloc(sizeof(double)*NR*params.nt)));
     MALLOC_SAFE( (tmpfld->efficiency = (double *)malloc(sizeof(double)*params.nt)));

     MALLOC_SAFE( (tmpfld->nu_grid =  (double *)malloc(sizeof(double)*NR)));
     MALLOC_SAFE( (tmpfld->grid_torque =  (double *)malloc(sizeof(double)*(NR+1))));
     MALLOC_SAFE( (tmpfld->grid_torquec =  (double *)malloc(sizeof(double)*(NR))));
    for(i=0;i<NR;i++) {
        tmpfld->lami[i] = 0;
        tmpfld->mdoti[i] = 0;
        tmpfld->nu_grid[i] = 0;
        tmpfld->grid_torque[i] = 0;
        tmpfld->grid_torquec[i] = 0;
    }

    for(i=0;i<params.nt;i++) {
        tmpfld->times[i] = 0;
        tmpfld->avals[i] = 0;
        tmpfld->vs[i] = 0;
        tmpfld->vs_ss[i] = 0;
        tmpfld->mdot_ss[i] = 0;
        tmpfld->efficiency[i] = 0;

    }

    for(i=0;i<NR*params.nt;i++) {
        tmpfld->sol[i]= 0;
        tmpfld->sol_mdot[i]=0;
        tmpfld->lamp[i] = 0;
        tmpfld->lam0[i] = 0;
        tmpfld->sol_ss[i] = 0;
        tmpfld->ivals_ss[i] = 0;
        tmpfld->kvals_ss[i] = 0;
        tmpfld->dTr[i] = 0;
        tmpfld->dep_func[i] = 0;

    }
    return;
}

void free_field( Field *tmpfld) {

    free(tmpfld->times);
    free(tmpfld->avals);
    free(tmpfld->vs);
    free(tmpfld->sol);
    free(tmpfld->sol_mdot);
    free(tmpfld->torque);
    free(tmpfld->lami);
    free(tmpfld->mdoti);
    free(tmpfld->mdot_ss);
    free(tmpfld->vs_ss);
    free(tmpfld->sol_ss);
    free(tmpfld->lamp);
    free(tmpfld->lam0);
    free(tmpfld->efficiency);
    free(tmpfld->nu_grid);
    free(tmpfld->ivals_ss);
    free(tmpfld->kvals_ss);
    free(tmpfld->dTr);
    free(tmpfld->dep_func);
    free(tmpfld->grid_torque);
    free(tmpfld->grid_torquec);

    return;
}

