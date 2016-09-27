#include "pdisk.h"



void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;
   // double em,ec,ep;
   // double lm,lc,lp;
    double rm, rp;
    double am,ap,num,nup,bm,bp,drm,drp;
 //   double dp_tot, dm_tot, dp,dm,dc;

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

#ifdef _OPENMP
//#pragma omp parallel for private(i,rm,rp,lm,lc,lp,em,ec,ep,dm,dc,dp,dm_tot,dp_tot) shared(matrix,rc,rmin)
//#pragma omp parallel for private(i,rm,rp,num,nup,am,bm,ap,bp,drm,drp)
#endif
    for(i=1;i<NR-1;i++) {
        rm = rmin[i];
        rp = rmin[i+1];

        num = nu(rm);
        nup = nu(rp);

        am = 3*num;
        bm = 3*(num/rm)*(params.gamma - .5);

        ap = 3*nup;
        bp = 3*(nup/rp)*(params.gamma-.5);

        drm = rc[i]-rc[i-1];
        drp = rc[i+1]-rc[i];

        matrix.md[i] = ( bp * ( rc[i+1] - rp) - ap) / drp - ( am + bm * ( rm - rc[i-1] ) ) / drm;
        matrix.ud[i] = ( ap + bp * ( rp - rc[i] ) ) / drp;
        matrix.ld[i-1] = - ( bm *(rc[i] - rm ) - am ) / drm; 


/*
        if (params.planet_torque) {
            if (params.forced_torque) {
                 matrix.fm[i] += 2*sqrt(rm) * 2*M_PI*rm * fld.grid_torque[i] - 2*sqrt(rp) * 2 *M_PI*rp * fld.grid_torque[i+1];
            }
            else {
            }
        }
*/ 
    }
    
    
    if (params.planet_torque && planet.nonlocal_torque ) {
        set_uw(matrix.u,matrix.w,aplanet,NR);
    }



    set_boundary();
    // Coefficient Matrix is all set

    if (params.planet_torque && planet.nonlocal_torque) {
        matvec_full(matrix.ld,matrix.md,matrix.ud,matrix.u,matrix.w,y,matrix.fm,dt,dt/2.,NR,2);
    }
    else {
        matvec(matrix.ld,matrix.md,matrix.ud,y,matrix.fm,dt,dt/2.,NR);
    }


#ifdef _OPENMP
#pragma omp parallel for private(i) shared(y,dr,matrix)
#endif
    for(i=0;i<NR-1;i++) {
        matrix.fm[i] += dr[i]*y[i];
        matrix.md[i] = dr[i] - matrix.md[i]*dt/2.;
        matrix.ld[i] *= -dt/2.;
        matrix.ud[i] *= -dt/2.;
        matrix.u[i] *= -dt/2;
        matrix.u[i+NR] *= -dt/2;
    }
    matrix.fm[NR-1] += dr[NR-1]*y[NR-1];
    matrix.md[NR-1] = dr[NR-1] - matrix.md[NR-1]*dt/2.;
    matrix.u[NR-1] *= -dt/2.;
    matrix.u[NR-1+NR] *= -dt/2.;

 

    if (params.planet_torque && planet.nonlocal_torque) {
        trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,matrix.u,matrix.w,NR);
    }
    else {
        trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
    }
    return;
}

/*
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
*/

void explicit_step_func(double aplanet, double *y,double *mdL, double *mdR) {
    int i;
    
    double TL = 0;
    double TR = 0;
    double facm,facp,am,ap,bm,bp,nup,num,drm,drp;
    double fac,rm,rp,mfacL,mfacR,lfacL,ufacR;
    double yL,yM,yR;

    double rgL = .5*(rmin[0] + exp(log(rmin[0]) - dlr));
    double rgR = .5*(rmin[NR] + exp(log(rmin[NR]) + dlr));
    double bc_fac_L = 2*M_PI*rgL;
    double bc_fac_R = 2*M_PI*rgR;

    for(i=0;i<NR;i++) {

        fac = dr[i]*dTr_ex(rc[i],aplanet)*y[i];
        if (rc[i] < aplanet) TL += fac;
        else    TR += fac;
    }

    for(i=0;i<NR;i++) {
        rm = rmin[i];
        rp = rmin[i+1];
        facm = dep_func(rm,aplanet,planet.xd,planet.wd);
        facp = dep_func(rp,aplanet,planet.xd,planet.wd);
        facm *= (rm < aplanet) ? TL : TR; 
        facp *= (rp > aplanet) ? TR : TL; 

        num = nu(rm);
        nup = nu(rp);

        am = 3*num;
        bm = 3*(num/rm)*(params.gamma - .5);

        ap = 3*nup;
        bp = 3*(nup/rp)*(params.gamma-.5);

        drm = rc[i]-rc[i-1];
        drp = rc[i+1]-rc[i];

        mfacR = ( bp * ( rc[i+1] - rp) - ap) / drp;
        mfacL = ( am + bm * ( rm - rc[i-1] ) ) / drm;
        ufacR  = ( ap + bp * ( rp - rc[i] ) ) / drp;
        lfacL  = ( bm *(rc[i] - rm ) - am ) / drm; 
        yL = (i==0) ? params.bc_val[0]*bc_fac_L: y[i-1] ;

        yM = y[i];
        yR = (i==NR-1) ? params.bc_val[2]*bc_fac_R: y[i+1];
            
        mdL[i] = (lfacL*yL + mfacL*yM) - 2*sqrt(rp)*facp;
        mdR[i] = (mfacR*yM + ufacR*yR) - 2*sqrt(rm)*facm;
        //dy[i] = dt*(lfac*yL + mfac*yM + ufac*yR) ;
        //dy[i] += dt*( -2*sqrt(rm)*facm + 2*sqrt(rp)*facp);


    }

    return;
}
