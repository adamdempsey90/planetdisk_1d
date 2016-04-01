#include "pdisk.h"
double dTr(double x,double a) {
    double left_fac, right_fac; 
    double norm, xi,res;

    if (planet.gaussian == TRUE) {
        xi = (x-a)/planet.dep;
        left_fac = (xi-planet.beta)/planet.delta;
        right_fac = (xi+planet.beta)/planet.delta;
        left_fac = exp(-left_fac*left_fac);
        right_fac = exp(-right_fac*right_fac);

        norm = planet.T0/(planet.delta * sqrt(M_PI));

        res = -norm*( (planet.G1+1)*right_fac - left_fac);
    }
   
    else {
        xi = (x-a) / scaleH(a);

 //       norm = planet.eps * a*M_PI*(planet.mp*params.mth)*(planet.mp*params.mth); 
        norm = planet.eps * (planet.mp*params.mth)*(planet.mp*params.mth)/2.;
        if (planet.symmetric_torque) {
            right_fac = norm*pow(a/fmax(scaleH(x),fabs(x-a)),4);
        }
        else {
            right_fac = norm*pow(x/fmax(scaleH(x),fabs(x-a)),4);
        }
        left_fac = -norm*pow(x/fmax(scaleH(x),fabs(x-a)),4);    
        
        left_fac *= (1-smoothing(xi, -planet.c, planet.delta));
        right_fac *= smoothing(xi,-planet.c, planet.delta)*smoothing(xi,planet.c,planet.delta);
        
        res = left_fac*(1-planet.onesided) + right_fac;


    
    }

    return res;
}

double smoothing(double x, double x0, double w) {
    return 0.5*(1 + tanh( (x-x0)/w));
}


double calc_total_torque(double a, double *y) {
    int i;
    double res = 0;

#ifdef OPENMP
#pragma omp parallel for private(i) reduction(+:res) 
#endif
    for(i=0;i<NR;i++) {
        if (params.nonlocal_torque) {
            res += weights[i]*dTr_nl(rc[i],i,a,TRUE)*y[i];  // Factor of 1/2 from torque normalization

        }
        else {
            res += weights[i]*dTr(rc[i],a)*y[i];  // Factor of 1/2 from torque normalization
        }
    }
//    res *= -2*dlr*sqrt(a)/(planet.mp*params.mth);

    return res*dlr;

}

