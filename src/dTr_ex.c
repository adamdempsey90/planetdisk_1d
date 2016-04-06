#include "pdisk.h"

double dTr_ex(double x,double a) {
    double left_fac, right_fac; 
    double norm, xi,res;

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


    

    return res;
}

double smoothing(double x, double x0, double w) {
    return 0.5*(1 + tanh( (x-x0)/w));
}


