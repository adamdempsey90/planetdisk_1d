#include "pdisk.h"


void check_floor(double *y) {
    int i;

    if (params.use_floor) {
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
        for(i=0;i<NR;i++) {
            y[i] = fmax(y[i],2*M_PI*rc[i]*params.mass_floor);
        }
    }
    return;
}


void set_coeffs_matrix(void) {
    int i;
    double rm, rp;
    double am,ap,num,nup,bm,bp,drm,drp;

#ifdef _OPENMP
#pragma omp parallel for private(rm,rp,num,nup,am,bm,ap,bp,drm,drp,i)
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


    }
    return;
}

void explicit_step(double dt, double c1, double c2, double aplanet, double *y, double *res) {
    /* Explicit step, y^{n+1} = c1*y^n + c2*dt*A*y^n
     */
    int i;
    double TL,TR;
    if (params.density_dep) {
        TL = calc_inner_torque(planet.a,y);
        TR = calc_outer_torque(planet.a,y);
    }
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
/* Set the diffusion coefficient matrix */
    set_coeffs_matrix(); 

    
    if (params.planet_torque && planet.nonlocal_torque) {
        if (params.density_dep) {
            set_density_dep(matrix.ld, matrix.md,matrix.ud,y,TL, TR, aplanet);
        }  
        else {
            set_uw(matrix.u,matrix.w,aplanet,NR);
        }
    }


    set_boundary();
    // Coefficient Matrix is all set

    if (params.planet_torque && planet.nonlocal_torque && !params.density_dep) {
        matvec_full(matrix.ld,matrix.md,matrix.ud,matrix.u,matrix.w,y,matrix.fm,dt*c2,dt*c2,NR,2);
    }
    else {
        matvec(matrix.ld,matrix.md,matrix.ud,y,matrix.fm,dt*c2,dt*c2,NR);
    }

    /* Add final result */
#ifdef _OPENP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        res[i] = y[i]*c1 + matrix.fm[i]/dr[i];;
    }

    return;
}
void steady_state_step(double aplanet, double *y) {
    int i;

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
/* Set the diffusion coefficient matrix */
    set_coeffs_matrix(); 

    
    if (params.planet_torque && planet.nonlocal_torque ) {
        set_uw(matrix.u,matrix.w,aplanet,NR);
    }



    set_boundary();
    for(i=0;i<NR;i++) {
        matrix.fm[i] *= -1;
    }
    // Coefficient Matrix is all set
    

    if (params.planet_torque && planet.nonlocal_torque) {
        trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,matrix.u,matrix.w,NR);
    }
    else {
        trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
    }
    check_floor(y);

    return;
}

void crank_nicholson_step(double dt, double aplanet, double *y) {
    int i;

    double TL, TR;
    if (params.density_dep) {
        TL = calc_inner_torque(planet.a,y);
        TR = calc_outer_torque(planet.a,y);
    }
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
/* Set the diffusion coefficient matrix */
        
    set_coeffs_matrix(); 
    
    
    if (params.planet_torque && planet.nonlocal_torque) {
        if (params.density_dep) {
            set_density_dep(matrix.ld, matrix.md,matrix.ud,y,TL, TR, aplanet);
        }  
        else {
            set_uw(matrix.u,matrix.w,aplanet,NR);
        }
    }



    set_boundary();
    // Coefficient Matrix is all set

    if (params.planet_torque && planet.nonlocal_torque && !params.density_dep) {
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

 

    if (params.planet_torque && planet.nonlocal_torque && !params.density_dep) {
        trisolve_sm2(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,matrix.u,matrix.w,NR);
    }
    else {
        trisolve(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,NR);
    }
    check_floor(y);
    return;
}




void tvd_step(double dt, double aplanet, double *y) {
/* TVD step
 * y' = y + dt*f(y)
 * y_new = .5*y + .5*(y' + dt*f(y'))
 */

    int i;
    double *res = (double *)malloc(sizeof(double)*NR);
    double *res2 = (double *)malloc(sizeof(double)*NR);
    explicit_step(dt, 1.0, 1.0,aplanet, y, res);

    check_floor(res);
    explicit_step(dt,1.0,1.0,aplanet,res,res2);
    check_floor(res2);

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<NR;i++) {
        y[i] = .5*y[i] + .5*res2[i];;
    }

    check_floor(y);


    free(res);
    free(res2);

    return;
}
