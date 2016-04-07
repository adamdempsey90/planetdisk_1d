#include "pdisk.h"


void set_mdot(int planet_torque) {
    int i;
    double lamm;
#ifndef NONLOCAL
   for(i=1;i<NR;i++) { 
        
        mdot[i] = 1.5*( lam[i]*nu(rc[i])/sqrt(rc[i]) - lam[i-1]*nu(rc[i-1])/sqrt(rc[i-1]))/(sqrt(rc[i])-sqrt(rc[i-1]));

        if (planet_torque) {
            lamm = lam[i-1]*dr[i]/(dr[i-1]+dr[i]) + lam[i]*dr[i-1]/(dr[i-1]+dr[i]);
            mdot[i] -= 2*sqrt(rmin[i])*dTr_ex(rmin[i],planet.a)*lamm;
        }
        

    }
    mdot[0] = lam[0]*1.5*nu(rmin[0])/rmin[0];
#endif
    return;
}

