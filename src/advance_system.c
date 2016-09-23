#include "pdisk.h"

double set_dt(void) {
    int i;
    double dt1,dt2,dt3;
    double dt_min=params.tvisc;
    for(i=0;i<NR;i++) {
        dt1 = dlr/(rc[i]*cs(rc[i]));
        dt2 = dlr*dlr/nu(rc[i]);
        if (dt1 < dt_min) {
            dt_min = dt1;
        }
        if (dt2 < dt_min) {
            dt_min = dt2;
        }
    }

    return dt_min*params.cfl;
}


void advance_system(double dt, double *t, double tend) { 
    double newa, newv;
    double dt2 = tend-(*t);

/*
    if (params.explicit_stepper) {
        dt = set_dt();
    }
*/
    if (dt > dt2) {
        dt = dt2;
    }

    newa = planet.a;
    newv = planet.vs;

    
    if ((params.move_planet) && (*t >= params.release_time)) {
        move_planet(dt,lam,&newv,&newa);
    }

    crank_nicholson_step(dt,planet.a,lam);

    if ((params.move_planet) && (*t >= params.release_time)) {
        planet.a = newa;
        planet.vs = calc_drift_speed(planet.a,lam);
        set_planet_deposition();
    }

    *t += dt;

    return;
}


