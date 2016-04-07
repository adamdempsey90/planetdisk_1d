#include "pdisk.h"
int in_region(double x, double a, double leftr, double rightr) {

        double dist = fabs(x-a);
        if ((dist >= leftr) && (dist <= rightr)) {
            return TRUE;
        }
        else {
            return FALSE;
        }
}

double dep_func(double x, double a, double xd, double w) {
    double dist = x-a;
/*
    if (dist < 0) dist += xd*w;
    else        dist -= xd*w;
    dist *= dist;
    return exp(- dist * 2*M_PI*M_PI/(w*w))/(w);
*/
    double leftr = (xd-.5)*w;
    double rightr = (xd+.5)*w;
        
    double norm = 1./w;
    double res = smoothing(x,a+leftr,params.h/2)*(1-smoothing(x,a+rightr,params.h/2));
    res += smoothing(x,a-rightr,params.h/2)*(1-smoothing(x,a-leftr,params.h/2));
    return norm*res;




}

