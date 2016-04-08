#include "pdisk.h"
void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;
    double em,ec,ep;
    double lm,lc,lp;
    double rm, rp;
    double dp_tot, dm_tot, dp,dm,dc;

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
#pragma omp parallel for private(i,rm,rp,lm,lc,lp,em,ec,ep,dm,dc,dp,dm_tot,dp_tot) shared(matrix,rc,rmin)
#endif
    for(i=1;i<NR-1;i++) {
        rm = rmin[i];
        rp = rmin[i+1];

        lm = sqrt(rc[i-1]);
        lc = sqrt(rc[i]);
        lp = sqrt(rc[i+1]);

        em = 1.5*nu(rc[i-1])/sqrt(rc[i-1]);
        ec = 1.5*nu(rc[i])/sqrt(rc[i]);
        ep = 1.5*nu(rc[i+1])/sqrt(rc[i+1]);

        dp_tot = dr[i] + dr[i+1];
        dm_tot = dr[i] + dr[i-1];
        dp = dr[i+1];
        dc = dr[i];
        dm = dr[i-1];


        matrix.md[i] = -ec*(1./(lc-lm) + 1./(lp-lc));
        matrix.ld[i-1] = em/(lc-lm);
        matrix.ud[i] = ep/(lp-lc);


#ifndef NONLOCAL
        if (params.planet_torque) {

            em = 2*sqrt(rm)*dTr_ex(rm,aplanet);
            ep = 2*sqrt(rp)*dTr_ex(rp,aplanet);
            matrix.ld[i-1] += em*dc/dm_tot;
            matrix.ud[i] -= ep*dc/dp_tot;
            matrix.md[i] += em*dm/dm_tot - ep*dp/dp_tot;
        }
#endif
        
    }
    
    
#ifdef NONLOCAL
    if (params.planet_torque) {
        set_uw(matrix.u,matrix.w,aplanet,NR);
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
            matrix.ld[i] *= -dt/2.;
            matrix.ud[i] *= -dt/2.;
        }
    }

 

#ifdef NONLOCAL
    trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,matrix.u,matrix.w,NR);
#else
    trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
#endif

    return;
}

