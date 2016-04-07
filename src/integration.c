#include "pdisk.h"
void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;
    double ap,am, bm,bp, rm, rp;
    double dp_tot, dm_tot, dp,dm;
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

#ifdef OPENMP
#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(matrix,rc,rmin)
#endif
    for(i=0;i<NR;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        dp_tot = dr[i] + dr[i+1];
        dm_tot = dr[i] + dr[i-1];
        dp = dr[i+1];
        dm = dr[i];


        am = 3*nu(rm) * (params.gamma - .5)/(rm);
        bm = 6*nu(rm);
        ap = 3*nu(rp) * (params.gamma - .5)/(rp);
        bp = 6*nu(rp);

#ifndef NONLOCAL
        if (params.planet_torque) {
            am -= 2*sqrt(rm)*dTr_ex(rm,aplanet);
            ap -= 2*sqrt(rp)*dTr_ex(rp,aplanet);
        }
#endif
        matrix.md[i] = (dp*ap-bp)/dp_tot - (bm+dm*am)/dp_tot;
        if (i>0) {
            matrix.ld[i-1] = (bm - dm*am)/dm_tot;
        }
        if (i<NR) {
            matrix.ud[i] = (bp + dm*ap)/dp_tot;
        }
    }
#ifdef NONLOCAL
    if (params.planet_torque) {
        set_uw(matrix.u,matrix.w,aplanet);
    }
#endif
    set_boundary();

    // Coefficient Matrix is all set

#ifdef NONLOCAL    
    matvec_full(matrix.ld,matrix.md,matrix.ud,matrix.u,matrix.w,y,matrix.fm,dt,dt/2.,NR,2);
#else
    matvec(matrix.ld,matrix.md,matrix.ud,y,matrix.fm,dt,dt/2.,NR);
#endif


#ifdef OPENMP
#pragma omp parallel for private(i) shared(y,dr,matrix)
#endif
    for(i=0;i<NR;i++) {
        matrix.fm[i] += dr[i]*y[i];
        matrix.md[i] = dr[i] - matrix.md[i]*dt/2.;
        if (i<NR-1) {
            matrix.ld[i] *= -1;
            matrix.ud[i] *= -1;
        }
    }

 

#ifdef NONLOCAL
    trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,matrix.u,matrix.w,y,NR);
#else
    trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
#endif

    return;
}

