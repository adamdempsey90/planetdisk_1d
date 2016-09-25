#include "pdisk.h"


void set_mdot(int planet_torque) {
    int i;
    double lamm;
    if (planet_torque && planet.nonlocal_torque) {
        set_torque_nl(planet.a,lam,mdot,TRUE);
    }
    else {
        zero_array(mdot,NR);
    }
    
    double rm,num,am,bm,drm;
    for(i=1;i<NR;i++) { 


        rm = rmin[i];
        num = nu(rm);
        am = 3*num;
        bm = 3*(num/rm)*(params.gamma-.5);

        drm = rc[i]-rc[i-1];
        
        //mdot[i] = -mdot[i]*2*sqrt(rmin[i]);
        mdot[i] = -2*sqrt(rmin[i])*mdot[i] +  (bm*(rc[i]-rm)-am)*lam[i-1]/drm + (am + bm*(rm-rc[i-1]))*lam[i]/drm; 

        //mdot[i] = -mdot[i]*2*sqrt(rmin[i]) + 1.5*( lam[i]*nu(rc[i])/sqrt(rc[i]) - lam[i-1]*nu(rc[i-1])/sqrt(rc[i-1]))/(sqrt(rc[i])-sqrt(rc[i-1]));
/*
        if (planet_torque && !planet.nonlocal_torque) {
            lamm = lam[i-1]*dr[i]/(dr[i-1]+dr[i]) + lam[i]*dr[i-1]/(dr[i-1]+dr[i]);
            mdot[i] -= 2*sqrt(rmin[i])*dTr_ex(rmin[i],planet.a)*lamm;
        }
*/      

    }

    mdot[0] = get_inner_bc_mdot(lam[0]);
    
    //mdot[0] = -mdot[0]*2*sqrt(rc[i]) + lam[0]*1.5*nu(rmin[0])/rmin[0];
    /*
    if (planet_torque && params.forced_torque) {
        for(i=0;i<NR;i++) {
            mdot[i] -= 2*sqrt(rmin[i]) * 2*M_PI*rmin[i] * fld.grid_torque[i];
        }
    }
    */
    return;
}

