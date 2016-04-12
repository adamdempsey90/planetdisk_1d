#include "pdisk.h"



double dTr_linear(double x,double a) {
    double left_fac, right_fac; 
    double norm, xi,res;

    xi = (x-a) / scaleH(a);

 //       norm = planet.eps * a*M_PI*(planet.mp*params.mth)*(planet.mp*params.mth); 
    norm = planet.eps * (planet.mp*params.mth)*(planet.mp*params.mth)/2.;
    if (planet.symmetric_torque) {
        right_fac = norm*pow(a/fmax(scaleH(a),fabs(x-a)),4);
    }
    else {
        right_fac = norm*pow(x/fmax(scaleH(a),fabs(x-a)),4);
    }
    left_fac = -norm*pow(x/fmax(scaleH(a),fabs(x-a)),4);    

/*
#ifndef TRANSLATE
    left_fac *= (1-smoothing(xi, -(planet.c+planet.xd), planet.delta));
    right_fac *= smoothing(xi,-(planet.c+planet.xd), planet.delta)*smoothing(xi,(planet.c+planet.xd),planet.delta);
#else 
*/
    left_fac *= (1-smoothing(xi, -planet.c, planet.delta));
    right_fac *= smoothing(xi,-planet.c, planet.delta)*smoothing(xi,planet.c,planet.delta);
    res = left_fac*(1-planet.onesided) + right_fac;


    

    return res/x;
}



double dTr_ex(double x, double a) {

#ifdef TRANSLATE
    double xd = planet.xd*params.h*a;
    if (x > a) {
        xd = fmax(a,x-xd);
    }
    else {
        xd = fmin(a,x+xd);
    }
    return dTr_linear(xd,a);
#else
    
    return dTr_linear(x,a);

#endif


}
