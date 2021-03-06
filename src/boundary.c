#include "pdisk.h"


/*
void set_boundary(void) {
    int i;
    double em,ec,ep;
    double rb = rmin[0];
    double r0 = rc[0];


    i=0;
    ec = 1.5*nu(rc[i])/sqrt(rc[i]);
    ep = 1.5*nu(rc[i+1])/sqrt(rc[i+1]);
    matrix.md[i] = -ec/(sqrt(rc[i+1]) - sqrt(rc[i]))  - 1.5*nu(rb)/rb;
    matrix.ud[i] = ep/(sqrt(rc[i+1]) - sqrt(rc[i]));

    i=NR-1;
    ec = 1.5*nu(rc[i])/sqrt(rc[i]);
    em = 1.5*nu(rc[i-1])/sqrt(rc[i-1]);
    matrix.md[i] = -ec/(sqrt(rc[i]) - sqrt(rc[i-1]));
    matrix.ld[i-1] = em/(sqrt(rc[i]) - sqrt(rc[i-1]));

    matrix.fm[NR-1] = params.bc_mdot;

    matrix.u[0] = 0;
    matrix.u[NR] = 0;

   
    return;

}
*/
double mixed_inner_boundary_value(double c1, double c2, double y0) {
/* Inner boundary
 * B.C of y' + c1*y  = c2
 */
    double rm = rmin[0];
    double r0 = rc[0];


    double num = nu(rm);
    double am = 3*num;
    double bm = 3*(num/rm)*(params.gamma - .5);

    double fac1 = bm/am;
    double fac2 = ( 1- fac1*(r0 - rm))/(1 - c1*(r0-rm));


    
    return  y0*am*(fac1 - c1*fac2) + am * c2*fac2;


}
double mixed_outer_boundary_value(double c1, double c2,double y0) {
/* Outer boundary
 * B.C of y' + c1*y  = c2
 */
    double rp = rmin[NR];
    double rn = rc[NR-1];


    double nup = nu(rp);
    double ap = 3*nup;
    double bp = 3*(nup/rp)*(params.gamma - .5);

    double fac1 = bp/ap;
    double fac2 = ( 1 - fac1*(rn - rp))/(1 - c1*(rn-rp));

    return y0 * ap*(fac1 - c1*fac2) + c2*ap*fac2; 


}

void mixed_inner_boundary(double c1, double c2) {
/* Inner boundary
 * B.C of y' + c1*y  = c2
 */
    double rm = rmin[0];
    double rp = rmin[1];
    double r0 = rc[0];

    double drp = rc[1]-rc[0];

    double num = nu(rm);
    double nup = nu(rp);
    double am = 3*num;
    double bm = 3*(num/rm)*(params.gamma - .5);
    double ap = 3*nup;
    double bp = 3*(nup/rp)*(params.gamma - .5);

    double fac1 = bm/am;
    double fac2 = ( 1- fac1*(r0 - rm))/(1 - c1*(r0-rm));


    
    matrix.md[0] = (bp*(rc[1]-rp) - ap)/drp - am*(fac1 - c1*fac2);
    matrix.ud[0] = (ap + bp*(rp-rc[0]))/drp;
    matrix.fm[0] = -am * c2*fac2;
    matrix.u[0] = 0;

    return;

}
void mixed_outer_boundary(double c1, double c2) {
/* Outer boundary
 * B.C of y' + c1*y  = c2
 */
    double rm = rmin[NR-1];
    double rp = rmin[NR];
    double rn = rc[NR-1];

    double drm = rc[NR-1]-rc[NR-2];

    double num = nu(rm);
    double nup = nu(rp);
    double am = 3*num;
    double bm = 3*(num/rm)*(params.gamma - .5);
    double ap = 3*nup;
    double bp = 3*(nup/rp)*(params.gamma - .5);

    double fac1 = bp/ap;
    double fac2 = ( 1 - fac1*(rn - rp))/(1 - c1*(rn-rp));

    matrix.md[NR-1] = ap*(fac1 - c1*fac2) - (bm*(rm - rc[NR-2])  + am)/drm; 
    matrix.ld[NR-2] = -(bm * (rc[NR-1] - rm) - am)/drm;
    matrix.fm[NR-1] = c2*ap*fac2; 
    matrix.u[NR-1] = 0;

    return;

}

void fixed_inner_boundary_mdot(double val) {
    double rm = rmin[0];
    double num = nu(rm);
    double am = 3*num;
    double bm = 3*(num/rm)*(params.gamma - .5);
    mixed_inner_boundary(bm/am, val/am);
    return;
}
void fixed_outer_boundary_mdot(double val) {
    double rp = rmin[NR];
    double nup = nu(rp);
    double ap = 3*nup;
    double bp = 3*(nup/rp)*(params.gamma - .5);
    mixed_outer_boundary(bp/ap, val/ap);
    return;
}



void fixed_inner_boundary_lam(double val) {
    double rm = .5*(rmin[0] + exp( log(rmin[0]) - dlr));
    double dr1 = rm - rmin[0];
    mixed_inner_boundary(1./dr1, val/dr1);
    return;
}
void fixed_outer_boundary_lam(double val) {
    double rp = .5*(rmin[NR] + exp( log(rmin[NR]) + dlr));
    double dr1 = rp - rmin[NR];
    mixed_outer_boundary(1./dr1, val/dr1);
    return;
}
void fixed_inner_boundary_sigma(double val) {
    double rm = .5*(rmin[0] + exp( log(rmin[0]) - dlr));
    fixed_inner_boundary_lam(2*M_PI*rm*val);
    return;
}
void fixed_outer_boundary_sigma(double val) {
    double rp = .5*(rmin[NR] + exp( log(rmin[NR]) + dlr));
    fixed_outer_boundary_lam(2*M_PI*rp*val);
    return;
}
void fixed_inner_boundary_grad(double val) {
    mixed_inner_boundary(0, val);
    return;
}
void fixed_outer_boundary_grad(double val) {
    mixed_outer_boundary(0,val);
    return;
}
void fixed_inner_boundary_zero_torque(double val) {
    double rm = .5*(rmin[0] + exp( log(rmin[0]) - dlr));
    val = params.bc_val[2]/(3*M_PI*params.alpha*params.h*params.h*pow(rm,params.gamma));
    fixed_inner_boundary_sigma(val);
    return;
}

double fixed_inner_boundary_mdot_value(double val,double y0) {
    double rm = rmin[0];
    double num = nu(rm);
    double am = 3*num;
    double bm = 3*(num/rm)*(params.gamma - .5);
    return mixed_inner_boundary_value(bm/am, val/am,y0);
}
double fixed_outer_boundary_mdot_value(double val,double y0) {
    double rp = rmin[NR];
    double nup = nu(rp);
    double ap = 3*nup;
    double bp = 3*(nup/rp)*(params.gamma - .5);
    return mixed_outer_boundary_value(bp/ap, val/ap,y0);
}

double fixed_inner_boundary_lam_value(double val,double y0) {
    double rm = .5*(rmin[0] + exp( log(rmin[0]) - dlr));
    double dr1 = rm - rmin[0];
    return mixed_inner_boundary_value(1./dr1, val/dr1,y0);
}
double fixed_outer_boundary_lam_value(double val,double y0) {
    double rp = .5*(rmin[NR] + exp( log(rmin[NR]) + dlr));
    double dr1 = rp - rmin[NR];
    return mixed_outer_boundary_value(1./dr1, val/dr1,y0);
}
double fixed_inner_boundary_grad_value(double val, double y0) {
    return mixed_inner_boundary_value(0, val,y0);
}
double fixed_outer_boundary_grad_value(double val, double y0) {
    return mixed_outer_boundary_value(0,val,y0);
}
    

double get_inner_bc_mdot(double y0) {
    switch(params.bc_type[0]) {

        case BCMDOTIN:
            return fixed_inner_boundary_mdot_value(params.bc_val[0],y0);
        case BCLAMIN:
            return fixed_inner_boundary_lam_value(params.bc_val[0],y0);
        case BCGRADIN:
            return fixed_inner_boundary_grad_value(params.bc_val[0],y0);
        default:
            // Mixed BC
            return mixed_inner_boundary_value(params.bc_val[0],params.bc_val[1],y0);
    }
}
double get_outer_bc_mdot(double y0) {

    switch(params.bc_type[1]) {

        case BCMDOTOUT:
            return fixed_outer_boundary_mdot_value(params.bc_val[2],y0);
        case BCLAMOUT:
            return fixed_outer_boundary_lam_value(params.bc_val[2],y0);
        case BCGRADOUT:
            return fixed_outer_boundary_grad_value(params.bc_val[2],y0);
        default:
            // Mixed BC
            return mixed_outer_boundary_value(params.bc_val[2],params.bc_val[3],y0);
    }

}

void set_boundary(void) {
    /* Set the boundary conditions
     */
        
    switch(params.bc_type[0]) {

        case BCMDOTIN:
            fixed_inner_boundary_mdot(params.bc_val[0]);
            break;
        case BCLAMIN:
            fixed_inner_boundary_lam(params.bc_val[0]);
            break;
        case BCGRADIN:
            fixed_inner_boundary_grad(params.bc_val[0]);
            break;
        case BCSIGIN:
            fixed_inner_boundary_sigma(params.bc_val[0]);
        case BCZTIN:
            fixed_inner_boundary_zero_torque(params.bc_val[0]);
            break;
        default:
            // Mixed BC
            mixed_inner_boundary(params.bc_val[0],params.bc_val[1]);
    }

    switch(params.bc_type[1]) {

        case BCMDOTOUT:
            fixed_outer_boundary_mdot(params.bc_val[2]);
            break;
        case BCLAMOUT:
            fixed_outer_boundary_lam(params.bc_val[2]);
            break;
        case BCGRADOUT:
            fixed_outer_boundary_grad(params.bc_val[2]);
            break;
        case BCSIGOUT:
            fixed_outer_boundary_sigma(params.bc_val[2]);
            break;
        default:
            // Mixed BC
            mixed_outer_boundary(params.bc_val[2],params.bc_val[3]);
    }




    return;
}
