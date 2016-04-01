#include "pdisk.h"
double tau_integrand(double x) {
    double norm = 2 * pow(params.h,-2.5)*pow(2.,-1.25);
    double p1 = params.gamma;
    double p2 = -.5*(params.gamma - 1.5);
    
    return norm* pow(x,(5*p2+p1)*.5 - 11./4)*pow(fabs(pow(x,1.5)-1),1.5);
}

void set_tau(double a) {
    /* Set the value of tau
     * This assumes that matrix.icol has already been set to the 
     * correct index!
     */
    int i;
    int icol = matrix.icol;

    /* taumin[icol] is tau at the inner zone edge that holds the planet */
    taumin[icol] = 0;
    tauc[icol] = 0;
    /* Outer disk */
    taumin[icol+1] = .5*dlr*tau_integrand(rmin[icol+1]/a)*rmin[icol+1]/a;
    tauc[icol+1] = .5*dlr*tau_integrand(rc[icol+1]/a)*rc[icol+1]/a;
    for(i=icol+2;i<NR-1;i++) {
        taumin[i] = taumin[i-1] + dlr*tau_integrand(rmin[i]/a)*rmin[i]/a;
        tauc[i] = tauc[i-1] + dlr*tau_integrand(rc[i]/a)*rc[i]/a;

    }
    taumin[NR-1] = taumin[NR-2] + .5*dlr*tau_integrand(rmin[NR-1]/a)*rmin[NR-1]/a;
    tauc[NR-1] = tauc[NR-2] + .5*dlr*tau_integrand(rc[NR-1]/a)*rc[NR-1]/a;
    /* Inner Disk */
    taumin[icol-1] = .5*dlr*tau_integrand(rmin[icol-1]/a)*rmin[icol-1]/a;
    tauc[icol-1] = .5*dlr*tau_integrand(rc[icol-1]/a)*rc[icol-1]/a;
    for(i=icol-2;i>0;i--) {
        taumin[i] = taumin[i+1] + .5*dlr*tau_integrand(rmin[i]/a)*rmin[i]/a;
        tauc[i] = tauc[i+1] + .5*dlr*tau_integrand(rc[i]/a)*rc[i]/a;
    }
    taumin[0] = taumin[1] + .5*dlr*tau_integrand(rmin[1]/a)*rmin[1]/a;
    tauc[0] = tauc[1] + .5*dlr*tau_integrand(rc[1]/a)*rc[1]/a;
    return;
}

