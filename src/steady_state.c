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
    tmpfld->a = planet.a;
    tmpfld->vs = 0;
    tmpfld->mdot0 = params.bc_mdot;
    tmpfld->mdot = params.bc_mdot;
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
    resp = -2*dr[0]*dTr_ex(rc[0],a)*sqrt(rc[0])/(3*nu(rc[0]));
    for(i=1;i<NR;i++) {
        res = -2*dr[i]*dTr_ex(rc[i],a)*sqrt(rc[i])/(3*nu(rc[i]));
        tmpfld->ivals[i] = tmpfld->ivals[i-1] + .5*(res+resp);
        resp = res;
    }

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        tmpfld->ivals[i] = exp(tmpfld->ivals[i]); // * pow(rc[i]/params.ri,params.gamma-.5);
    }


    tmpfld->kvals[0] = 0;
    resp = dr[0] * (tmpfld->ivals[0])/(2*sqrt(rc[0]));
    for(i=1;i<NR;i++) {
//        res = dr[i]*tmpfld->ivals[i]/(3*nu(rc[i]));
        res = dr[i] * (tmpfld->ivals[i])/(2*sqrt(rc[i]));
        tmpfld->kvals[i] = tmpfld->kvals[i-1] + .5*(res+resp);
        resp = res;
    }

    tmpfld->mdot = params.bc_mdot;
    tmpfld->mdot0 = params.bc_mdot;
    for(i=0;i<NR;i++) {
        tmpfld->lam0[i] = params.bc_mdot*2*rc[i]/(3*nu(rc[i])); 
    }
        
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
         tmpfld->lam[i] = (1 + tmpfld->kvals[i]/sqrt(params.ri))/(tmpfld->ivals[i]);
         tmpfld->lam[i] *= tmpfld->lam0[i]*sqrt(params.ri/rc[i]);
    }
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
     
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

void steadystate_config_nl(SteadyStateField *tmpfld,double a) {
    int i,j;
    double w, xd;
    double c1, c2, x;
    double mdot_fac, lam_fac;
    double F_fac;
    double res;
    double dtr;
    double Fa;

    double TL, TR;
    mdot_fac = params.bc_mdot;
    xd = planet.xd;
    w = params.h*a;


    tmpfld->mdot = mdot_fac;
    tmpfld->mdot0 = mdot_fac;

    c1=0;
    c2=1.0;
    for(i=0;i<NR;i++) {
        tmpfld->lam0[i] = (2./3)*mdot_fac*rc[i]/nu(rc[i]);
    }

    for(i=0;rc[i]<=a;i++) {
        x = (rc[i] - a)/w;
        x += (xd + planet.c);
        x = x/sqrt(2);
        F_fac = 1.5*nu(rc[i])/sqrt(rc[i]);
        dtr = dTr_ex(rc[i],a)/F_fac;
        c1 += dr[i]*dtr* sqrt(rc[i]);
        c2 -= dr[i]*dtr*erfc(-x)/2.;
    }
    Fa = c1/(c2*sqrt(a));
    for(i=0;rc[i]<=a;i++) {
        x = (rc[i] - a)/w;
        x += (xd + planet.c);
        x = x/sqrt(2);
        tmpfld->lamp[i] = .5*erfc(-x)*c1/( c2*sqrt(rc[i])) ;
        tmpfld->lam[i] = (tmpfld->lamp[i] + 1)*tmpfld->lam0[i];
    }


    c1 = 0;
    c2 = 1.0;

    for(i=NR-1;rc[i]>a;i--) {
        x = (rc[i] - a)/w;
        x -= (xd + planet.c);
        x /= sqrt(2);
        F_fac = 1.5*nu(rc[i])/sqrt(rc[i]);
        dtr = dTr_ex(rc[i],a)/F_fac;
        c1 += dr[i]*dtr* (sqrt(a)*Fa + sqrt(rc[i]));
        c2 -= dr[i]*dtr*erfc(-x)/2.;
    }

    for(i=NR-1;rc[i]>a;i--) {
        x = (rc[i] - a)/w;
        x -= (xd + planet.c);
        x /= sqrt(2);
        tmpfld->lamp[i] = Fa*sqrt(a/rc[i]) + .5*erfc(-x)*c1/(c2*sqrt(rc[i])) ;
        tmpfld->lam[i] = (tmpfld->lamp[i] + 1)*tmpfld->lam0[i];
    }


    tmpfld->vs = tmpfld->lamp[NR-1]*sqrt(rc[NR-1]);
    tmpfld->vs *= -params.bc_mdot * 2*sqrt(a)/(planet.mp*params.mth);

    for(i=0;i<NR;i++) {
        tmpfld->ivals[i] = 0;
        tmpfld->kvals[i] = 0;
    }

    return;
}

