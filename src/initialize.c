#include "pdisk.h"

void init_lam_from_file(void) {
    FILE *f = fopen("lambda_init.dat","r");
    if (f == NULL) {
        printf("Can't find the data file lambda_init.dat!");
    }
    double x,temp;
    int i=0;
    while (fscanf(f,"%lg\t%lg\n",&x,&temp) != EOF) {
          if (fabs(x -rc[i]) < 1e-4) {
                 lam[i] = temp;
          }
          else {
             printf("Grid is not the same! %.12e\t%.12lg\t.12%lg\n", fabs(x-rc[i]),rc[i],x);
          }
        i++;
    }
 
    fclose(f);

    set_mdot(params.planet_torque);
    return;

}


void init_lam(void) {
    int i;
    printf("Init lam\n");
    for(i=0;i<NR;i++) {
        lam[i] = params.bc_mdot*2*rc[i]/(3*nu(rc[i]));
    }
    printf("Done\n");
    return;
}
/*
void init_lam(void) {
    int i;
    double mdot0, sig0, x;

    if (params.flux_bc) {
        for(i=0;i<NR;i++) {
     //       lam[i] = params.bc_mdot*2*rc[i]*(1 - sqrt(params.ri/rc[i]))/(3*nu(rc[i]));
            lam[i] = params.bc_mdot*2*rc[i]/(3*nu(rc[i]));
        }
    }
    else {
         mdot0 = 1.5*params.nu0 * (params.bc_lam[1]*pow(params.ro,params.gamma - .5) - params.bc_lam[0]*pow(params.ri,params.gamma-.5))/(sqrt(params.ro)-sqrt(params.ri));
     sig0 = 2*mdot0/(3*params.nu0);

     for(i=0;i<NR;i++) {
         x = rc[i]/params.ri;
         lam[i] = params.bc_lam[0]*pow(x,.5-params.gamma) + sig0*pow(rc[i],1-params.gamma)*(1-sqrt(1/x));
     }
    }
#ifdef INITIAL_NOISE
    double locs[7] = {.1,.6, 2, 5, 9, 20 , 30};
    int nlocs = 7; 
    int j;
    for(i=0;i<NR;i++) {
        for(j=0;j<nlocs;j++) {
            lam[i] += params.bc_lam[1] * exp( -(rc[i]-locs[j])*(rc[i]-locs[j])/.2);
        }
    }
#endif
    set_mdot(FALSE);
    return;
}
*/
