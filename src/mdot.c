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

/*
void set_mdot(int planet_torque) {
    int i;
    double ca, cb;


    for(i=0;i<NR;i++) {
        ca = 3*nu(rc[i])/rc[i];
        cb = ca*(params.gamma - .5);

        if (planet_torque) {
                cb -= 2*dTr_ex(rc[i],planet.a)/(sqrt(rc[i]));
        }  
        
        mdot[i] = cb*lam[i]; 
        //printf("%d\t%.2e\n",i,mdot[i]/params.bc_mdot);
        if (i == 0){ 
            if (params.flux_bc) {
                mdot[i] = params.bc_mdot;
                mdot[i] = mdot[1];
            }
            else {
                mdot[i] += ca*(lam[i+1]-params.bc_lam[0])/(dlr);
            }
        }
        else {
            if (i==NR-1) {
                if (params.flux_bc) {
                    mdot[i] = params.bc_mdot;
                }
                else {
                    mdot[i] += ca*(params.bc_lam[1] - lam[i-1])/(2*dlr);
                }
            }
            else {
                mdot[i] += ca*(lam[i+1]-lam[i-1])/(2*dlr);
            }
        }
    }



    return;
}
*/
