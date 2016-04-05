#include "pdisk.h"
double cs(double x) {
    return params.h*pow(x,(2*params.gamma - 3.)/4.);

}
double nu(double x) {
      double res = params.nu0 * pow(x,params.gamma);

      if (params.hs_visc) {
          double rh = pow(planet.mp*params.mth/3.,1./3) * planet.a;
          if (fabs(x-planet.a) <= rh) {
            res += 1.5*pow(rh/planet.a,3) * sqrt(planet.a)/(2*M_PI);
          }

      }

      return res;
}
double scaleH(double x) {
    return params.h * x *  pow(x, (params.gamma -.5)/2);
}
