#include "pdisk.h"


void set_planet_deposition(void) {

    if (planet.scaling_dep) {
        double gdepth = 1./(1 + .04*planet.K);
        double kp = pow(planet.K * params.h*params.h,.25); 
        double dr1 = (gdepth/4. + .08)*kp*planet.a;
        double dr2 = .33*kp*planet.a;
        planet.xd =  .5*(dr1 +dr2);
        planet.wd = (dr2-dr1);
        
        printf("Gap depth should be %.2e\n",gdepth); 
    }

    return;
}

int in_region(double x, double a, double leftr, double rightr) {

        double dist = fabs(x-a);
        if ((dist >= leftr) && (dist <= rightr)) {
            return TRUE;
        }
        else {
            return FALSE;
        }
}


void set_dep_box(double a,double xd, double wd, double *box_region) {
    int i,j;
    int have_set[4] = {FALSE, FALSE, FALSE, FALSE};
    box_region[0] = a - xd - wd/2;
    box_region[1] = a - xd + wd/2;
    box_region[2] = a + xd - wd/2;
    box_region[3] = a + xd + wd/2;

    printf("\n%lg\t%lg\t%lg\n",a,xd,wd);
    for(i=0;i<4;i++) printf("%lg\t",box_region[i]);
    printf("\n");

    for(i=0;i<NR;i++) {
        for(j=0;j<4;j++) {
            if (!have_set[j]) {
                if (box_region[j] >= rmin[i]) {
                    box_region[j] = rmin[i];
                    have_set[j] = TRUE;
                }
            }
        }
    }

    return;
}

double dep_func(double x, double a, double xd, double w) {
    if (xd == 0) {
        return dTr_ex(x,a);
    }

    double distL,distR, res;

    double wd;
    if (planet.scaling_dep) {
        wd = w/(2*4);
    }
    else {
        wd = w;
    }

    wd *= wd;


   // if (x <=a) {
   //     xd *= -1;
    //}
    distL = (x-(a+xd));
    distR = (x-(a-xd));
    distL *= distL;
    distR *= distR;
    return exp(-distL/(2*wd))/sqrt(2*M_PI*wd) + exp( - distR/(2*wd))/sqrt(2*M_PI*wd); 
/*
    if (x<a) {

        if ((xd -.5*w < -fabs(x-a)) && (-fabs(x-1)<xd+.5*w)){
            return 1./w;
        }
        else {
            return 0;
        }
    }
    else {

        if ((xd -.5*w < fabs(x-a)) && (fabs(x-1)<xd+.5*w)){
            return 1./w;
        }
        else {
            return 0;
        }
    }
*/

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




void set_uw(double *u, double *w, double a, int n) {
/* Set the u w vectors for non-local deposition.
 * There are two of each, u[i],u[i+n],w[i],w[i+n], 
 * for the inner and outer disks.
 */
    int i;

    double rm, rp;
    double facm, facp;

#ifdef _OPENMP
#pragma omp parallel for private(i,rm,rp,facm,facp) 
#endif
    for(i=0;i<n;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        
        if (rc[i] < a) {
        /* Inner disk */
            if (planet.symmetric_torque) {
                w[i] = 0;
            }
            else {
                w[i] = dr[i]*dTr_ex(rc[i],a); // lower integral weights
            }
            w[i+n] = 0; // upper integral weights
            
            facm = dep_func(rm,a,planet.xd,planet.wd);
            facp = dep_func(rp,a,planet.xd,planet.wd);
            u[i] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;
            u[i+n] = 0;

        }
        else {
        /* Outer disk */
            w[i+n] = dr[i]*dTr_ex(rc[i],a); // upper integral weights
            if (planet.symmetric_torque) {
                w[i] = -w[i+n]; // lower integral weights
            }
            else {
                w[i] = 0;
            }

            facm = dep_func(rm,a,planet.xd,planet.wd);
            facp = dep_func(rp,a,planet.xd,planet.wd);
            u[i+n] = -2*sqrt(rp)*facp + 2*sqrt(rm)*facm;

            u[i] = 0;

        }


    }
    return;

}
void set_torque_nl(double a, double *y, double *res, int edge) {
    int i;
    double TL, TR;

    //TL = calc_inner_torque(a,y);
    //TR = calc_outer_torque(a,y);

    TL = 0;
    TR = 0;
//#ifdef _OPENMP
//#pragma omp parallel for reduction(+:TL,TR) private(i)
//#endif
    for(i=0; i < NR; i++) {
        if (rc[i] < a) {
            TL += dr[i] * dTr_ex(rc[i],a) * y[i];
        }
        else {
            TR += dr[i] * dTr_ex(rc[i],a) * y[i];
        }
    }
    if (planet.symmetric_torque) {
        TL = - TR;
    }
    //printf("total torque TL = %lg, TR = %lg, dT = %lg \n",TL,TR,TL+TR);

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        if (rc[i] < a) {
            if (edge) {
                res[i] = TL*dep_func(rmin[i],a,planet.xd,planet.wd);
            }
            else {
                fld.dep_func[i] = dep_func(rc[i],a,planet.xd,planet.wd);
                res[i] = TL* fld.dep_func[i];
            }
        }
        else {
            if (edge) {
                res[i] = TR*dep_func(rmin[i],a,planet.xd,planet.wd);
            }
            else {
                fld.dep_func[i] = dep_func(rc[i],a,planet.xd,planet.wd);
                res[i] = TR*fld.dep_func[i];
            }
        }

    }

    return;
}

void calculate_linear_torque(double aplanet, double *y, double *TL, double *TR) {
    int i;    
    *TL = 0;
    *TR = 0;

    for(i=0; i < NR; i++) {
        if (rc[i] < aplanet) {
            *TL += dr[i] * dTr_ex(rc[i],aplanet) * y[i];
        }
        else {
            *TR += dr[i] * dTr_ex(rc[i],aplanet) * y[i];
        }
    }
    return;
}


void set_torque_linear(double aplanet, double *y) {
    int i;
    double TL,TR;
    double rm, rp;
    double facm,facp;
    calculate_linear_torque(aplanet, y,&TL,&TR);

#ifdef _OPENMP
#pragma omp parallel for private(i,rm,rp,facm,facp) 
#endif
    for(i=0;i<NR;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        
        if (rc[i] < aplanet) { 
            facm = dep_func(rm,aplanet,planet.xd,planet.wd);
            facp = dep_func(rp,aplanet,planet.xd,planet.wd);
            matrix.fm[i] = TL*(-2*sqrt(rp)*facp + 2*sqrt(rm)*facm);

        }
        else {
        /* Outer disk */
            facm = dep_func(rm,aplanet,planet.xd,planet.wd);
            facp = dep_func(rp,aplanet,planet.xd,planet.wd);
            matrix.fm[i] = TR*(-2*sqrt(rp)*facp + 2*sqrt(rm)*facm);

        }


    }

    return;
}
