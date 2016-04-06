#include "pdisk.h"
void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;
    double ap,am, bm,bp, rm, rp;

    for(i=0;i<NR;i++) {
        matrix.fm[i] = 0;
        matrix.md[i] = 0;
        matrix.u[i] = 0;
        matrix.u[i+NR] = 0;
        matrix.w[i] = 0;
        matrix.w[i+NR] = 0;

        if (i < NR-1) {
            matrix.ud[i] = 0;
            matrix.ld[i] = 0;
        }
    }



#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(matrix,rc,rmin)
    for(i=1;i<NR-1;i++) {
        rm = rmin[i];
        rp = rmin[i+1];

        am = 3*nu(rm) * (params.gamma - .5)/(rm);
        bm = 3*nu(rm)/dr[i];
        ap = 3*nu(rp) * (params.gamma - .5)/(rp);
        bp = 3*nu(rp)/dr[i];
        //calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,aplanet,params.planet_torque);
        //calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,aplanet,params.planet_torque);
#ifndef NONLOCAL
        if (params.planet_torque) {
            am -= 3*sqrt(rm)*dTr_ex(rm,aplanet);
            ap -= 3*sqrt(rp)*dTr_ex(rp,aplanet);
        }
#endif
        matrix.md[i] = (ap-am - bm - bp)*dt/2.;
        matrix.ld[i-1] = (-am + bm)*dt/2.;
        matrix.ud[i] = (ap + bp)*dt/2.;
    }
#ifdef NONLOCAL
    set_uw(matrix.u,matrix.w,aplanet);
#endif
    

    if (params.flux_bc) {
        am = 3*nu(rmin[1]) * (params.gamma - .5)/(rmin[1]);
        bm = 3*nu(rmin[1])/dr[i];
        ap = 3*nu(rmin[NR-1]) * (params.gamma - .5)/(rmin[NR-1]);
        bp = 3*nu(rmin[NR-1])/dr[i];
#ifndef NONLOCAL
        if (params.planet_torque) {
            am -= 3*sqrt(rmin[1])*dTr_ex(rmin[1],aplanet);
            ap -= 3*sqrt(rmin[NR-1])*dTr_ex(rmin[NR-1],aplanet);
        }
#endif
 //       matrix.md[0] = (ap-bp)*dt/2.;
 //       matrix.ud[0] = (ap+bp)*dt/2.;
        matrix.md[0] = 0; matrix.ud[0] = 0;
        matrix.md[NR-1] = (-am-bm)*dt/2.;
        matrix.ld[NR-2] = (-am+bm)*dt/2.;
 //       matrix.fm[0] = -params.bc_mdot*dt;
        matrix.fm[0] = 0;
        matrix.fm[NR-1] = params.bc_mdot*dt;
    }

    matvec(matrix.ld,matrix.md,matrix.ud,y,matrix.fm,NR);

#pragma omp parallel for private(i) shared(y,dr,matrix)
    for(i=0;i<NR;i++) {
        matrix.fm[i] += dr[i]*y[i];
        matrix.md[i] = dr[i] - matrix.md[i];
        if (i<NR-1) {
            matrix.ld[i] *= -1;
            matrix.ud[i] *= -1;
        }
    }

    if (params.flux_bc) {
        matrix.fm[0] = 0; //-params.bc_mdot;
//        matrix.fm[NR-1] = params.bc_mdot;
        matrix.md[0] = 1.0;
//        matrix.md[NR-1] = 0;
//        matrix.ld[NR-2] = 0;
        matrix.ud[0] = -sqrt(rc[0]/rc[1]);
    }
    else {
        matrix.md[NR-1] = 1;
        matrix.ld[NR-2] = 0;
        matrix.fm[NR-1] = params.bc_lam[1];

        matrix.md[0] = 1;
        matrix.ud[0] = 0;
        matrix.fm[0] = params.bc_lam[0];
    }
/*    
    else {
    //    matrix.md[0] = 0; //1;
        matrix.ud[0] = 0;
    //    matrix.fm[0] = 0; //params.bc_mdot*2*rc[0]/(3*nu(rc[0]));


        if (params.bc_lam[0] == 0) {
            matrix.md[0] = 1;
            matrix.ud[0] = 0;
            matrix.fm[0] = params.bc_lam[0];
        }

    }
*/

#ifdef NONLOCAL
    trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,matrix.u,matrix.w,y,NR);
#else
    trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
#endif

    return;
}

