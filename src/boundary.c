#include "pdisk.h"

void set_boundary(void) {
    
    double rb = rmin[0];
    double r0 = rc[0];
    double bb = 3*nu(rb);
    double ab = bb*(params.gamma - .5)/rb;

    matrix.md[0] += ab*(r0*nu(r0))/(rb*nu(rb)) + (3/rb)*(rb*nu(rb)-r0*nu(r0))/(r0-rb);
   // matrix.ud[0] = 0;
   // matrix.ld[NR-2] = 0;
   // matrix.md[NR-1] = 0;
    matrix.fm[NR-1] = params.bc_mdot;
    matrix.u[0] = 0;
    matrix.u[NR] = 0;
/*
   matrix.md[0] = 0;
   matrix.ud[0] = 0;
   matrix.md[NR-1] = 0;
   matrix.ld[NR-2] = 0;
   matrix.fm[NR-1] = params.bc_mdot;
  */ 
    return;

}
