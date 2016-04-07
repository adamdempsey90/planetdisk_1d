#include "pdisk.h"

void set_boundary(void) {
    int i;
    double em,ec,ep;
    double rb = rmin[0];
    double r0 = rc[0];


    i=0;
    ec = 1.5*nu(rc[i])/sqrt(rc[i]);
    ep = 1.5*nu(rc[i+1])/sqrt(rc[i+1]);
    matrix.md[i] = -ec/(sqrt(rc[i+1]) - sqrt(rc[i]))  - 1.5*nu(rb)/rb;
    matrix.ud[i] = ep/(sqrt(rc[i+1]) - sqrt(rc[i]));

    i=NR-1;
    ec = 1.5*nu(rc[i])/sqrt(rc[i]);
    em = 1.5*nu(rc[i-1])/sqrt(rc[i-1]);
    matrix.md[i] = -ec/(sqrt(rc[i]) - sqrt(rc[i-1]));
    matrix.ld[i-1] = em/(sqrt(rc[i]) - sqrt(rc[i-1]));

    matrix.fm[NR-1] = params.bc_mdot;

    matrix.u[0] = 0;
    matrix.u[NR] = 0;

   
    return;

}
