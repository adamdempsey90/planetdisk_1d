#include "pdisk.h"


void set_uw(double *u, double *w, double a, int n) {
/* Set the u w vectors for non-local deposition.
 * There are two of each, u[i],u[i+n],w[i],w[i+n], 
 * for the inner and outer disks.
 */
    int i;

    double hp = params.h*a;
    double xd = planet.xd;
    double leftr = (xd-.5)*hp;
    double rightr = (xd+.5)*hp;
    double dist;

#ifdef OPENMP
#pragma omp parallel for private(i,dist) 
#endif
    for(i=0;i<n;i++) {
        dist = fabs(rc[i]-a);
        if (rc[i] < a) {
        /* Inner disk */
            w[i] = dr[i]*dTr_ex(rc[i],a); // lower integral weights
            w[i+n] = 0; // upper integral weights

            if ( dist>leftr && dist<rightr) {
                u[i] = 1/hp;
            }
            else {
                u[i] = 0;
            }
            u[i+n] = 0;

        }
        else {
        /* Outer disk */
            w[i] = 0; // lower integral weights
            w[i+n] = dr[i]*dTr_ex(rc[i],a);; // upper integral weights

            if ( dist>leftr && dist<rightr) {
                u[i+n] = 1/hp;
            }
            else {
                u[i] = 0;
            }
            u[i] = 0;

        }


    }
    return;

}




