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

    if (params.explicit_stepper) {
        dt = set_dt();
    }

    if (dt > dt2) {
        dt = dt2;
    }

    newa = planet.a;
    newv = planet.vs;

    
    if ((params.move_planet) && (*t >= params.release_time)) {
        move_planet(dt,lam,&newv,&newa);
    }

    if (params.nonlocal_torque) {
        crank_nicholson_step_nl(dt,planet.a,lam);
    }
    else {
        crank_nicholson_step(dt,planet.a,lam);
    }

    if ((params.move_planet) && (*t >= params.release_time)) {
        planet.a = newa;
        planet.vs = calc_drift_speed(planet.a,lam);
    }
 /*
    if (params.explicit_stepper) {
        explicit_step(dt,&planet.a,lam);
        planet.vs = calc_drift_speed(planet.a,lam);
    }
    else {
    if ((params.move_planet) && (*t >= params.release_time)) {

        if (params.move_planet_implicit) {
     //       move_planet_implicit(dt,lam,&planet.vs, &planet.a);
      //      predictor_corrector(dt,lam, &planet.vs,&planet.a); 
            multi_step(dt,lam,&planet.vs,&planet.a);
        }
        else {
            if (params.nonlocal_torque) {
                crank_nicholson_step_nl(dt,planet.a,lam);
             }
            else {
                crank_nicholson_step(dt,planet.a,lam);
              }
            move_planet(dt, lam, &planet.vs, &planet.a); 
        }
    }
    else {
        if (params.nonlocal_torque) {
            crank_nicholson_step_nl(dt,planet.a,lam);
         }
        else {
            crank_nicholson_step(dt,planet.a,lam);
          }
        planet.vs = calc_drift_speed(planet.a,lam);
    }
    }
*/
    *t += dt;

    return;
}


