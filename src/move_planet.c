#include "pdisk.h"

double calc_drift_speed(double a,double *y) {
    double T = calc_total_torque(a,y);
    return -2*sqrt(a)*T/(planet.mp*params.mth);
}


void move_planet(double dt, double *y, double *vs, double *a) {
    double T = calc_total_torque(*a,y);
    double q = planet.mp*params.mth;

    double L0 = sqrt(*a) * q;
    double L1;
    L1 = L0 - dt*T;
    
    *vs = -2 * (*a) * T / L0;
    *a = pow(L1/q,2);
/*
    *vs = calc_drift_speed(*a,y);
    planet.a += dt*(*vs);
    *a += dt*(*vs);
*/
    return;
}

