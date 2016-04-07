#include "pdisk.h"
void calc_coeffs(double x, double dx,double *a, double *b,double rp,int planet_torque) {
    // dx = rc[i]-rc[i-1] or rc[i+1] - rc[i] or dr[i]
    *a = 3*nu(x) * (params.gamma - .5)/(x);
    *b = 3*nu(x)/(dx);

    if (planet_torque == TRUE) {
        *a -= dTr(x,rp)/( sqrt(x));
    }
//    *a /= 2;
    return;
}
