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


void set_uw(double *u, double *w, double a, int n) {
/* Set the u w vectors for non-local deposition.
 * There are two of each, u[i],u[i+n],w[i],w[i+n], 
 * for the inner and outer disks.
 */
    int i;

    double hp = params.h*a;
    double xd = planet.xd;
    double rm, rp;
    double facm, facp;

#ifdef OPENMP
#pragma omp parallel for private(i,rm,rp,facm,facp) 
#endif
    for(i=0;i<n;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        
        if (rc[i] < a) {
        /* Inner disk */
            w[i] = dr[i]*dTr_ex(rc[i],a); // lower integral weights
            w[i+n] = 0; // upper integral weights
            
            facm = dep_func(rm,a,planet.xd,hp);
            facp = dep_func(rp,a,planet.xd,hp);
            u[i] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;
            u[i+n] = 0;

        }
        else {
        /* Outer disk */
            w[i] = 0; // lower integral weights
            w[i+n] = dr[i]*dTr_ex(rc[i],a);; // upper integral weights

            facm = dep_func(rm,a,planet.xd,hp);
            facp = dep_func(rp,a,planet.xd,hp);
            u[i+n] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;

            u[i] = 0;

        }


    }
    return;

}
