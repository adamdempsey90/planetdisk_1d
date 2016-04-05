#include "pdisk.h"
double dTr_nl(double x, double a) {
        double tsh = 1.89 + .53/planet.mp;
        double f0 = planet.eps; //.45;
        double T0 = (planet.mp*params.mth)*(planet.mp*params.mth)/params.mth;
        double norm = sqrt(tsh) * f0 * T0 * .5/M_PI;
        norm /= a;
        
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
        if (x < a) {
            return -norm*tau_integrand(x/a)*pow(tauval,-1.5);//sqrt(x);
        }
        else {
            return norm*tau_integrand(x/a)*pow(tauval,-1.5);//sqrt(x);
        }
}
