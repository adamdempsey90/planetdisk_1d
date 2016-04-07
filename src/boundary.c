#include "pdisk.h"

void set_boundary(void) {
    
    double rb =rmin[0];
    double r0 = rc[0];
    double bb = 3*nu(rb);
    double ab = bb*(params.gamma - .5)/rb;

    mat.md[0] = ab*(r0*nu(r0))/(rb*nu(rb)) + (3/rb)*(rb*nu(rb)-r0*nu(r0))/(r0-rb);
    mat.ud[0] = 0;
    mat.ld[NR-2] = 0;
    mat.md[NR-1] = 0;
    mat.u[0] = 0;
    mat.u[NR] = 0;
    mat.fm[NR-1] = params.bc_mdot;

    return;

}
