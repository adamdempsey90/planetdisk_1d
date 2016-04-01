#include "pdisk.h"
void explicit_step(double dt, double *aplanet, double *y) {
/* TVD RK3 
 */
    int i;
    double a1, a2;
    double ap,am,bp,bm,rm,rp;
    double *y1 = (double *)malloc(sizeof(double)*NR);
    double *y2 = (double *)malloc(sizeof(double)*NR);
    
    a1 = *aplanet; a2 = *aplanet;

#pragma omp parallel for private(i)
    for(i=0;i<NR;i++) {
        y1[i] = 0;
        y2[i] = 0;
    }

#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(rc,rmin,y,y1,dt)
    for(i=0;i<NR;i++) {
        if (i==0) {
            am=0; bm=0;
        }
        else {
            rm = rmin[i];
            calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,*aplanet,params.planet_torque);
        }
        if (i!= NR-1) {
            rp = rmin[i+1];
            calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,*aplanet,params.planet_torque);
        }
        else {
            ap = 0; bp =0;
        }
        y1[i] = (ap-am-bm-bp)*y[i];
        if (i!=0) {
            y1[i] += (bm-am)*y[i-1];
        }
        else {
            y1[i] -= params.bc_mdot;
        }
        if (i!=NR-1) {
            y1[i] += (ap+bp)*y[i+1];
        }
        else {
            y1[i] += params.bc_mdot;
        }
        y1[i] = y[i] + dt*y1[i]/dr[i];
    }



    if (params.move_planet) {
        a1 = *aplanet + dt*calc_drift_speed(*aplanet,y);
        
    }

#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(rc,rmin,y,y2,dt)
    for(i=0;i<NR;i++) {
        if (i==0) {
            am=0; bm=0;
        }
        else {
            rm = rmin[i];
            calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,*aplanet,params.planet_torque);
        }
        if (i!= NR-1) {
            rp = rmin[i+1];
            calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,*aplanet,params.planet_torque);
        }
        else {
            ap = 0; bp =0;
        }

        y2[i] = (ap-am-bm-bp)*y1[i];
        if (i!=0) {
            y2[i] += (bm-am)*y1[i-1];
        }
        else {
            y2[i] -= params.bc_mdot;
        }
        if (i!=NR-1) {
            y2[i] += (ap+bp)*y1[i+1];
        }
        else {
            y2[i] += params.bc_mdot;
        }
        y2[i] = .75*y[i] + .25*y1[i] + .25*dt*y2[i]/dr[i];
    }
    if (params.move_planet) {
        a2 = .75*(*aplanet) + .25*a1 + .25*dt*calc_drift_speed(a1,y1);
    }

#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(rc,rmin,y,y1,y2,dt)
    for(i=0;i<NR;i++) {
        if (i==0) {
            am=0; bm=0;
        }
        else {
            rm = rmin[i];
            calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,*aplanet,params.planet_torque);
        }
        if (i!= NR-1) {
            rp = rmin[i+1];
            calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,*aplanet,params.planet_torque);
        }
        else {
            ap = 0; bp =0;
        }
    
        y[i] = (1./3)*y[i] + (2./3)*y2[i]; 
        y[i] += (2./3)*dt*((ap-am-bm-bp)*y2[i])/dr[i];
        
        if (i != 0) {
            y[i] += (2./3)*dt*(bm-am)*y2[i-1]/dr[i];
        }
        else {
            y[i] -= (2./3)*dt*params.bc_mdot/dr[i];
        }
        if (i != NR-1) {
            y[i] += (2./3)*dt*(ap+bp)*y2[i+1]/dr[i];
        }
        else {
            y[i] += (2./3)*dt*params.bc_mdot/dr[i];
        }
    }
    if (params.move_planet) {
        *aplanet = (1./3)*(*aplanet) + (2./3)*a2 + (2./3)*dt*calc_drift_speed(a2,y2);
    }
    

    free(y1); free(y2);
    return;
}

