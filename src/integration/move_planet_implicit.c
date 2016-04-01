#include "pdisk.h"

void move_planet_implicit(double dt, double *y, double *vs, double *a) {
    double tol = 1e-4;
    double a_old = *a;
    double args[2] = {dt,a_old};
    
    double *lam_old;
     MALLOC_SAFE(( lam_old =(double *)malloc(sizeof(double)*NR)));
    memcpy(&lam_old[0],&y[0],sizeof(double)*NR);
    *a = secant_method(&planet_zero_function_euler,a_old, .9*a_old, lam_old,tol,args); 
    memcpy(&y[0],&lam_old[0],sizeof(double)*NR);
    *vs = calc_drift_speed(planet.a,y);
    free(lam_old);
    /*
    double vs = calc_drift_speed(*a, );
    double rhs = planet.a + .5*dt * vs;
    double args[2] = {dt,rhs};
    double a_old = planet.a;
    double tol = 1e-4;
    planet.a = secant_method*(&planet_zero_function,a_old, .5*a_old,tol,args)
    */
    return;

}

double secant_method(double (*function)(double,double *,double[]),double x1, double x2,double *y,double tol, double args[]) {
    int i;

    
    double f1,f2, temp;

    f1 = (*function)(x1,y,args);
    for (i=0;i<MAXITERATIONS;i++) {
        f2 = (*function)(x2,y,args);
        if (fabs(f1-f2) <= tol) {
            printf("f1-f2 < tol:%lg\t%lg\t%lg\t%lg\n",x1,x2,f1,f2);
            break;
        }
        temp = x1;
        x1 -= f1*(x1-x2)/(f1-f2);
        if (fabs(x1-x2) <= tol) {
            printf("Converged to %lg in %d iterations\n",x1,i);
            return x1;
        }
       // printf("%lg\t%lg\t%lg\t%lg\n",x1,x2,f1,f2);
        x2 = temp;
        f1 =f2;
    }
    printf("Max iterations exceeded!\n");

    return x1;
}

double planet_zero_function_euler(double a, double *y,double args[]) {
    double  dt = args[0];
    double a_old = args[1];
    if (params.nonlocal_torque) {
        crank_nicholson_step_nl(dt,a,y);
    }
    else {
        crank_nicholson_step(dt,a,y);
    }
    double vs = calc_drift_speed(a,y);

    return a - a_old - dt*vs;
}

void predictor_corrector(double dt, double *y, double *vs, double *a) {

/*
    double vs1 = calc_drift_speed(*a, y);
    
    double a1 = *a + dt*vs1;

    crank_nicholson_step(dt,a1,y);
    
    *vs = calc_drift_speed(a1,y);
    *a += .5*dt*(vs1 + *vs);
*/

    predict_step(dt,y,vs,a);
    correct_step(dt,y,vs,a);
    predict_step(dt,y,vs,a);
    correct_step(dt,y,vs,a);


    return;
}

void predict_step(double dt, double *y,double *vs, double *a) {
    *vs = calc_drift_speed(*a,y);
    *a += dt * (*vs);
    return;
}

void correct_step(double dt, double *y, double *vs, double *a) {
    if (params.nonlocal_torque) {
        crank_nicholson_step_nl(dt,*a,y);
    }
    else {
        crank_nicholson_step(dt,*a,y);
    }
    double vs1 = calc_drift_speed(*a,y);
    *a += .5*dt*(vs1 - (*vs));
    *vs = vs1;
    return;
}

void multi_step(double dt, double *y, double *vs, double *a) {
    double vs1 = calc_drift_speed(*a,y);
    *a += .5*dt*vs1;

    if (params.nonlocal_torque) {
        crank_nicholson_step_nl(dt,*a,y);
    }
    else {
        crank_nicholson_step(dt,*a,y);
    }

    *vs = calc_drift_speed(*a,y);

    *a += .5*dt*(*vs);

    return;
}



