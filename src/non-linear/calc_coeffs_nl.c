#include "pdisk.h"
double calc_coeffs_nl(double xm,double rp,int i) {

/*
    double tsh = 1.89 + .53/planet.mp;
    double f0 = .45;
    double T0 = (planet.mp*params.mth)*(planet.mp*params.mth)*pow(params.h,-3.);
    double norm = -sqrt(tsh) * f0 * T0 * .5/M_PI;
    if (taumin[i] <= tsh) {
        return 0;
    }
    if (xm < rp) {
        return -norm*tau_integrand(xm/rp)*pow(taumin[i],-1.5); ///sqrt(xm);
    }
    else {
        return norm*tau_integrand(xm/rp)*pow(taumin[i],-1.5); ///sqrt(xm);
    } 
*/
    return -dTr_nl(xm,i,rp,FALSE)*sqrt(xm)*2;

}

