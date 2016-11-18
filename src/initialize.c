#include "pdisk.h"
//#include <gsl/gsl_sf_bessel.h>
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

void init_ring_test(void) {

    double t = 100.;
    double nu0 = nu(1.0);
    double r0 = 1.0;
    double tau = 12*nu0*t/(r0*r0);
    double m = 1;
    double norm = m/(M_PI*r0*r0) * (1/tau);
    double fac,u;
    int i;

    printf("initializing to viscously spreading ring\n");
    for(i=0;i<NR;i++) {
        u = rc[i]/r0;
        fac = 2*M_PI*rc[i];
        lam[i] = 0;
     //   lam[i] = norm * pow(u,-.25) * exp(-(1 + u*u)/tau) * gsl_sf_bessel_Inu(.25,2*u/tau);
    }

    return;
}
void init_constant_mdot_sig(void) {
    printf("Initializing profile to connstant mdot %.e\n", params.bc_val[2]);
    int i;
    double rg =.5*(rmin[0] + exp(log(rmin[0]) -dlr));
    for(i=0;i<NR;i++) {
        lam[i]  = 2*M_PI*rc[i]*params.bc_val[0]*nu(rg)/nu(rc[i]);
    }
    
    return;
}
void init_constant_mdot_lam(void) {
    printf("Initializing profile to connstant mdot %.e\n", params.bc_val[2]);
    int i;
    double rg =.5*(rmin[0] + exp(log(rmin[0]) -dlr));
    for(i=0;i<NR;i++) {
        lam[i]  = params.bc_val[0]*(nu(rg)/rg)/(nu(rc[i])/rc[i]);
    }
    
    return;
}

void init_constant_mdot(void) {
    printf("Initializing profile to connstant mdot %.e\n", params.bc_val[2]);
    int i;
    double zero_fac = ((params.bc_val[0] == 0)&&(params.bc_val[2]==0)) ? 1.0 : 0.0;
    for(i=0;i<NR;i++) {
        lam[i] = 2*M_PI*rc[i] * ( params.bc_val[2]/(3*M_PI*nu(rc[i])) + zero_fac);
    }

  

}
void init_grad_prof(void) {

    int i;
    for(i=0;i<NR;i++) {
        lam[i] = 1.0;
    }
    printf("Init to 1.\n");
    return;

}

void init_sigma_prof(void) {
    printf("Initializing profile to linear between %lg and %lg\n", params.bc_val[0],params.bc_val[2]);
    int i;
    double rp = exp( log(rmin[NR]) + dlr);
    double rm = exp( log(rmin[0]) - dlr);
    for(i=0;i<NR;i++) {
        lam[i] = 2*M_PI*((rm*params.bc_val[0] + (rp*params.bc_val[2]-rm*params.bc_val[0])*(rc[i] - rm)/(rp-rm)));
    }
    return;

}
void init_zero_torque(void) {
    int i;
    double rp = exp( log(rmin[NR]) + dlr);
    double rm = exp( log(rmin[0]) - dlr);
    printf("Initializing zero torque\n");
    params.bc_val[0] = params.bc_val[2]/(3*M_PI*params.alpha*params.h*params.h*pow(rm,params.gamma));
    for(i=0;i<NR;i++) {
        lam[i] = 2*M_PI*rc[i]*params.bc_val[0] * pow(rc[i]/rm,-params.gamma);
    }
    return;

}

void init_linear_prof(void) {
    printf("Initializing profile to linear between %lg and %lg\n", params.bc_val[0],params.bc_val[2]);
    int i;
    double rp = exp( log(rmin[NR]) + dlr);
    double rm = exp( log(rmin[0]) - dlr);
    for(i=0;i<NR;i++) {
        lam[i] = (params.bc_val[0] + (params.bc_val[2]-params.bc_val[0])*(rc[i] - rm)/(rp-rm));
    }


    double mfinal = 1.5*( params.bc_val[2]*nu(rp)/sqrt(rp) - params.bc_val[0]*nu(rm)/sqrt(rm))/(sqrt(rp)-sqrt(rm));
    printf("Final Mdot should be %lg\n",mfinal);
    return;
}

void init_lam(void) {
    int i;

    //init_ring_test();
    //
    if ((params.bc_type[0] == BCMDOTIN) && (params.bc_type[1] == BCMDOTOUT)) {
        init_constant_mdot();
    }
    else if ((params.bc_type[0] == BCLAMIN) && (params.bc_type[1] == BCLAMOUT)) {
        init_linear_prof();
    }
    else if ((params.bc_type[0] == BCSIGIN) && (params.bc_type[1] == BCSIGOUT)) {
        init_sigma_prof();
    }
    else if (params.bc_type[1] == BCGRADOUT) {
        init_grad_prof();
    }
    else if ( (params.bc_type[0] == BCLAMIN) && (params.bc_type[1] == BCMDOTOUT)) {
        init_constant_mdot_lam();
    }
    else if ( (params.bc_type[0] == BCSIGIN) && (params.bc_type[1] == BCMDOTOUT)) {
        init_constant_mdot_sig();
    }
    else if ( (params.bc_type[0] == BCZTIN) && (params.bc_type[1] == BCMDOTOUT)) {
        init_zero_torque();
    }
    printf("Init lam\n");
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
