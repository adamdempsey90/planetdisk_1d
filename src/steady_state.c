#include "pdisk.h"

void allocate_steady_state_field(SteadyStateField *tmpfld) {
    int i;
    MALLOC_SAFE( ( tmpfld->lam = (double *)malloc(sizeof(double)*NR)) );
    MALLOC_SAFE( ( tmpfld->lam0 = (double *)malloc(sizeof(double)*NR)) );
    MALLOC_SAFE( ( tmpfld->lamp = (double *)malloc(sizeof(double)*NR)) );
    MALLOC_SAFE(( tmpfld->ivals = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE(( tmpfld->kvals = (double *)malloc(sizeof(double)*NR)));


    for(i=0;i<NR;i++) {
        tmpfld->lam[i] = 0;
        tmpfld->lamp[i] = 0;
        tmpfld->lam0[i] = 0;
        tmpfld->ivals[i] = 0;
        tmpfld->kvals[i] = 0;
    }
    return;
}

void free_steady_state_field(SteadyStateField *tmpfld) {
    free(tmpfld->lam);
    free(tmpfld->lam0);
    free(tmpfld->lamp);
    free(tmpfld->ivals); 
    free(tmpfld->kvals);

    return;
}
void steadystate_config(SteadyStateField *tmpfld, double a) {
    int i;
    double res,resp;

    tmpfld->a = a;


    tmpfld->ivals[0] = 0;
    if (params.nonlocal_torque) {
        resp = -2*dlr*dTr_nl(rc[0],0,a,TRUE)*rc[0]*sqrt(rc[0])/(3*nu(rc[0]));

    }
    else {
        resp = -2*dlr*dTr(rc[0],a)*sqrt(rc[0])/(3*nu(rc[0]));
    }
    for(i=1;i<NR;i++) {
        if (params.nonlocal_torque) {
            res = -2*dlr*dTr_nl(rc[i],i,a,TRUE)*rc[i]*sqrt(rc[i])/(3*nu(rc[i]));

        }
        else {
            res = -2*dlr*dTr(rc[i],a)*sqrt(rc[i])/(3*nu(rc[i]));
        }
        tmpfld->ivals[i] = tmpfld->ivals[i-1] + .5*(res+resp);
        resp = res;
    }

#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        tmpfld->ivals[i] = exp(tmpfld->ivals[i]); // * pow(rc[i]/params.ri,params.gamma-.5);
    }


    tmpfld->kvals[0] = 0;
    resp = dlr * (tmpfld->ivals[0])*sqrt(rc[0])/2;
    for(i=1;i<NR;i++) {
//        res = dr[i]*tmpfld->ivals[i]/(3*nu(rc[i]));
        res = dlr * (tmpfld->ivals[i])*sqrt(rc[i])/2;
        tmpfld->kvals[i] = tmpfld->kvals[i-1] + .5*(res+resp);
        resp = res;
    }

    if (params.flux_bc) {
        tmpfld->mdot = params.bc_mdot;
    }
    else {
        tmpfld->mdot = (params.bc_lam[1]*tmpfld->ivals[NR-1] - params.bc_lam[0])/tmpfld->kvals[NR-1];
    }
    for(i=0;i<NR;i++) {
        tmpfld->lam0[i] = params.bc_mdot*2*rc[i]/(3*nu(rc[i])); 
    }
    if (params.flux_bc) {
        
#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
        for(i=0;i<NR;i++) {
         tmpfld->lam[i] = (1 + tmpfld->kvals[i]/sqrt(params.ri))/(tmpfld->ivals[i]);
         tmpfld->lam[i] *= tmpfld->lam0[i]*sqrt(params.ri/rc[i]);
/*
         tmpfld->lam[i] =  (tmpfld->mdot * tmpfld->kvals[i])/tmpfld->ivals[i];
         tmpfld->lam[i] += tmpfld->mdot * sqrt(params.ri)/tmpfld->ivals[i];
         tmpfld->lam[i] *= 2*sqrt(rc[i])/(3*nu(rc[i]));
*/
         }
    }
    else {
#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
        for(i=0;i<NR;i++) {
         tmpfld->lam[i] = (params.bc_lam[0] + tmpfld->mdot * tmpfld->kvals[i])/tmpfld->ivals[i];

         }
    }
    if (params.flux_bc) {
        tmpfld->mdot0 = params.bc_mdot;
    }
    else {
        tmpfld->mdot0 = 1.5*params.alpha*params.h*params.h;
        tmpfld->mdot0  *= (params.bc_lam[1]*pow(params.ro,params.gamma-.5)-params.bc_lam[0]*pow(params.ri,params.gamma-.5))/(sqrt(params.ro)-sqrt(params.ri));
    }
/*
    tmpfld->vs = (tmpfld->lam[NR-1])*nu(params.ro)*1.5/sqrt(params.ro)
                     -(tmpfld->mdot)*(sqrt(params.ro)-sqrt(params.ri));
    tmpfld->vs *=  -2*sqrt(planet.a)/(planet.mp*params.mth);
*/
    /*
    tmpfld->vs = -1.5*params.alpha*params.h*params.h*2*sqrt(a)*(params.bc_lam[1]-params.bc_lam[0]) *(1 - (tmpfld->mdot)/(tmpfld->mdot0));

    tmpfld->vs /= (planet.mp*params.mth);
    */
#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        
            //        tmpfld->lam0[i] += params.bc_lam[0] /pow(rc[i]/params.ri,params.gamma-.5);
     
        if (tmpfld->lam0[i] == 0) {
            tmpfld->lamp[i] = 0;
        }
        else {
            tmpfld->lamp[i] = (tmpfld->lam[i] - tmpfld->lam0[i])/(tmpfld->lam0[i]);
        }

    }
    tmpfld->vs = tmpfld->lamp[NR-1]*sqrt(rc[NR-1]);
    tmpfld->vs *= -params.bc_mdot * 2*sqrt(a)/(planet.mp*params.mth);
    return;

}

void steadystate_config_nl(SteadyStateField *tmpfld, double a) {
    int i;
    double res;

    tmpfld->a = a;

    double k = pow(planet.mp*params.mth,2)/(params.alpha*pow(params.h,5));
    k /= (3*M_PI);
    double f = 0;
    double fac = 0;
    double fp = planet.eps;
    double tauval = 0;
    double rval;
    double tsh = 1.89 + .53/planet.mp;
    set_tau(a);

    if (params.flux_bc) {
        tmpfld->mdot = params.bc_mdot;
        tmpfld->mdot0 = params.bc_mdot;
    }
#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        
        tmpfld->lam0[i] = (tmpfld->mdot0) * 2*rc[i] *( 1- sqrt(params.ri/rc[i]))/(3*nu(rc[i]));
//        tmpfld->lam0[i] += params.bc_lam[0] /pow(rc[i]/params.ri,params.gamma-.5);
       

    }

    for(i=0;i<NR;i++) {
        tauval = tauc[i];
        rval = rc[i];
        if (tauval <= tsh) {
            f = fp;
        }
        else {
            f = fp*sqrt(tsh/tauval);
        }

        fac = (sqrt(a)-sqrt(params.ri))/(sqrt(rval)-sqrt(params.ri));


        tmpfld->lamp[i] = -(k*f*fac/(1+k*fp));
        tmpfld->lam[i] = tmpfld->lam0[i]*(1+tmpfld->lamp[i]);
    
    }

    tmpfld->vs = tmpfld->lamp[NR-1];
    tmpfld->vs *= -2*sqrt(planet.a) * params.bc_mdot * (sqrt(rc[NR-1])-sqrt(params.ri));
    tmpfld->vs /= (planet.mp*params.mth);

    
    return;

}

