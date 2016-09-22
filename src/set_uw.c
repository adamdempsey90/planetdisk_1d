#include "pdisk.h"


int in_region(double x, double a, double leftr, double rightr) {

        double dist = fabs(x-a);
        if ((dist >= leftr) && (dist <= rightr)) {
            return TRUE;
        }
        else {
            return FALSE;
        }
}

double dep_func(double x, double a, double xd, double w) {
    double hp = scaleH(a);

    double sigma = w;
    sigma *= sigma;

    //double dist = (x-a)/hp;
    double dist = fabs(x-a);
    double res = 1./w;
    double leftr, rightr;
    if (xd == 0) {
        return dTr_ex(x,a);
    }
    else {
       /* 
        if (dist < 0) dist += (xd);
        else        dist -= (xd);
        dist *= dist;
        return exp(- dist/(2.*sigma)) / (sqrt(2*M_PI)*sigma);
       */ 
        if (x < a) {
            leftr = a - xd - w/2;
            rightr = a - xd + w/2;

            if ( (x > leftr) && (x < rightr) ) {
                return res;
            }
            else {
                return 0;
            }

        }
        else {
            
            leftr = a + xd - w/2;
            rightr = a + xd + w/2;

            if ( (x > leftr) && (x < rightr) ) {
                return res;
            }
            else {
                return 0;
            }
        }
        
    }
/*
    double leftr = (xd-.5)*w;
    double rightr = (xd+.5)*w;
    double res;   
    double norm = 1./w;

    if ( (fabs(x-a) > leftr) && (fabs(x-a) <rightr)) {
        res = norm;
    }
    else {
        res = 0;
    }
    return res;
*/
/*
    double res = smoothing(x,a+leftr,params.h/2)*(1-smoothing(x,a+rightr,params.h/2));
    res += smoothing(x,a-rightr,params.h/2)*(1-smoothing(x,a-leftr,params.h/2));
    return norm*res;
*/


}


void set_uw(double *u, double *w, double a, int n) {
/* Set the u w vectors for non-local deposition.
 * There are two of each, u[i],u[i+n],w[i],w[i+n], 
 * for the inner and outer disks.
 */
    int i;

    double hp = params.h*a;
    double gdepth = 1./(1 + .04*planet.K);
    double kp = pow(planet.K * params.h*params.h,.25); 
    double dr1 = (gdepth/4. + .08)*kp*a;
    double dr2 = .33*kp*a;
    double xd = .5*(dr1 +dr2);
    double wd = (dr2-dr1);
    double rm, rp;
    double facm, facp;


#ifdef OPENMP
#pragma omp parallel for private(i,rm,rp,facm,facp) 
#endif
    for(i=0;i<n;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        
        if (rc[i] < a) {
        /* Inner disk */
            w[i] = dr[i]*dTr_ex(rc[i],a); // lower integral weights
            w[i+n] = 0; // upper integral weights
            
            facm = dep_func(rm,a,xd,wd);
            facp = dep_func(rp,a,xd,wd);
            u[i] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;
            u[i+n] = 0;

        }
        else {
        /* Outer disk */
            w[i] = 0; // lower integral weights
            w[i+n] = dr[i]*dTr_ex(rc[i],a);; // upper integral weights

            facm = dep_func(rm,a,xd,wd);
            facp = dep_func(rp,a,xd,wd);
            u[i+n] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;

            u[i] = 0;

        }


    }
    return;

}
void set_torque_nl(double a, double *y, double *res) {
    int i;
    double TL, TR;
    double hp, xd;

    hp = params.h*a;
    xd = planet.xd;

    TL = calc_inner_torque(a,y);
    TR = calc_outer_torque(a,y);

#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        if (rc[i] <= a) {
            res[i] = TL*dep_func(rc[i],a,xd,hp);
        }
        else {
            res[i] = TR*dep_func(rc[i],a,xd,hp);
        }

    }

    return;
}
