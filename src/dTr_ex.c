#include "pdisk.h"


double fits[3][4] ={ { 0.0012, -.14, 0.0017, -0.115},
                    {0.08,-0.2,0.11,-0.25},
                    {0.7,0.06,0.8,0.12}};


double dTr_kanagawa(double x, double a) {
    
    double norm = 2*M_PI*.4*planet.eps*planet.q*planet.q;
    
    if (x < a) norm *= -1;

    double fac = 2*M_PI*x;

    double rfac = (x>=a) ? 1 : x/a;
    double dist = fabs(x/a-1);
    

/*
    if ((x>=rmin[400])&&(x<=rmin[401])) {
        return 1.0;
    }
    
    return 0.0;
    */
    if (dist > 1.3*params.h) {
        return norm*pow(dist/rfac,-4.)/fac;
    }



    
    return 0;

}

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


double dTr_fit(double x, double a) {
    
    double eps,mu,s;

    double alpha = params.alpha;
    double hp = params.h * a;
    double dist;
    
    dist = (x-a)/hp;


    if (dist >= 2./3) {

        eps = fits[0][2] * pow(alpha,fits[0][3]);
        mu = fits[1][2] * pow(alpha,fits[1][3]);
        s = fits[2][2] * pow(alpha,fits[2][3]);

    }
    else if (dist <= -2./3) {
        eps = -fits[0][0] * pow(alpha,fits[0][1]);
        mu = fits[1][0] * pow(alpha,fits[1][1]);
        s = fits[2][0] * pow(alpha,fits[2][1]);
        
    }
    else {
        return 0;
    }
    
    dist = log(fabs(dist)) - mu;
    return eps*exp(-dist*dist/(2*s*s))/(sqrt(2*M_PI)*s);
}

double dTr_ex(double x, double a) {

    return dTr_kanagawa(x,a);
/*
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
#ifdef TORQUEFIT
    return dTr_fit(x,a);
#else
    return dTr_linear(x,a);
#endif
#endif
*/

}
