#include "pdisk.h"
double dTr_nl(double x, int i, double rp, int centered) {

    if (params.shock_dep) {

        double tsh = 1.89 + .53/planet.mp;
        double f0 = planet.eps; //.45;
        double T0 = (planet.mp*params.mth)*(planet.mp*params.mth)/params.mth;
        double norm = sqrt(tsh) * f0 * T0 * .5/M_PI;
        norm /= rp;
        
        double tauval;
        if (centered) {
            tauval = tauc[i];
        }
        else {
            tauval= taumin[i];
        }


        if (tauval <= tsh) {
            return 0;
        }
        if (x < rp) {
            return -norm*tau_integrand(x/rp)*pow(tauval,-1.5);//sqrt(x);
        }
        else {
            return norm*tau_integrand(x/rp)*pow(tauval,-1.5);//sqrt(x);
        } 
    }
    else {
///        double fac = dTr(x,rp)/x;
        double nlfac=0;
        double xd = planet.xd;
        double hp = params.h*rp;
        double left_region = (xd-1)*hp/2;
        double right_region = (xd+1)*hp/2;
        if ( (fabs(x-rp) > left_region) && (fabs(x-rp) < right_region)) {
       //     printf("%lg inside %lg and %lg\n",(x-rp)/hp,left_region,right_region);
            nlfac = 1./hp;
        }
        return nlfac;

    }
}
double dTr_kana(double x, double a) {
    
    double norm = planet.eps*.798*pow(planet.mp*params.mth,2)*pow(a,4);
    double dist = x-a;

    if (fabs(dist) < params.h*a*1.3) {
        return 0;
    }
    else {
        if (dist<0) norm *= -1;
        return norm/(x*pow(fabs(x-a),4));
    }


}

