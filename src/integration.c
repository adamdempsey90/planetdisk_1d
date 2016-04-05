#include "pdisk.h"
void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;
    double ap,am, bm,bp, rm, rp;

    for(i=0;i<NR-1;i++) {
        matrix.fm[i] = 0;
        matrix.md[i] = 0;
        matrix.ud[i] = 0;
        matrix.ld[i] = 0;
    }
    matrix.md[NR-1] = 0;
    matrix.fm[NR-1] = 0;



#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(matrix,rc,rmin)
    for(i=1;i<NR-1;i++) {
        rm = rmin[i];
        rp = rmin[i+1];

        calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,aplanet,params.planet_torque);
        calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,aplanet,params.planet_torque);

        matrix.md[i] = (ap-am - bm - bp)*dt/2.;
        matrix.ld[i-1] = (-am + bm)*dt/2.;
        matrix.ud[i] = (ap + bp)*dt/2.;
      }
    if (params.flux_bc) {
        calc_coeffs(rmin[NR-1],rc[NR-1]-rc[NR-2],&am,&bm,aplanet,params.planet_torque);
        calc_coeffs(rmin[1],rc[1]-rc[0],&ap,&bp,aplanet,params.planet_torque);
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
    trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);


    return;
}

