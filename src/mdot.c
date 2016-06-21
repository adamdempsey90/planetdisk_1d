#include "pdisk.h"


void set_mdot(int planet_torque) {
    int i;
    double lamm;
#ifdef NONLOCAL
    if (planet_torque && !params.forced_torque) {
        set_torque_nl(planet.a,lam,mdot);
    }
#else
    zero_array(mdot,NR);
#endif
    
    for(i=1;i<NR;i++) { 

        mdot[i] = -mdot[i]*2*sqrt(rmin[i]) + 1.5*( lam[i]*nu(rc[i])/sqrt(rc[i]) - lam[i-1]*nu(rc[i-1])/sqrt(rc[i-1]))/(sqrt(rc[i])-sqrt(rc[i-1]));

        if (planet_torque && !params.forced_torque) {
#ifndef NONLOCAL        
            lamm = lam[i-1]*dr[i]/(dr[i-1]+dr[i]) + lam[i]*dr[i-1]/(dr[i-1]+dr[i]);
            mdot[i] -= 2*sqrt(rmin[i])*dTr_ex(rmin[i],planet.a)*lamm;
#endif
        }
        

    }
    mdot[0] = -mdot[0]*2*sqrt(rc[i]) + lam[0]*1.5*nu(rmin[0])/rmin[0];
    if (planet_torque && params.forced_torque) {
        for(i=0;i<NR;i++) {
            mdot[i] -= 2*sqrt(rmin[i]) * 2*M_PI*rmin[i] * fld.grid_torque[i];
        }
    }
    return;
}

