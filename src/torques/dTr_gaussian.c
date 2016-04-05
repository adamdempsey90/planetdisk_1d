#include "pdisk.h"

double dTr(double x,double a) {
    double left_fac, right_fac; 
    double norm, xi,res;

     xi = (x-a)/planet.dep;
     left_fac = (xi-planet.beta)/planet.delta;
     right_fac = (xi+planet.beta)/planet.delta;
     left_fac = exp(-left_fac*left_fac);
     right_fac = exp(-right_fac*right_fac);

     norm = planet.T0/(planet.delta * sqrt(M_PI));

     res = -norm*( (planet.G1+1)*right_fac - left_fac);
    return res;
}


