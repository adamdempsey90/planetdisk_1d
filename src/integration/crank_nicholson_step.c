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

void crank_nicholson_step_nl(double dt, double aplanet, double *y) {
    int i;
    int need_set_flag = TRUE;
    double ap,am, bm,bp, rm, rp;
    
    for(i=0;i<NR-1;i++) {
        matrix.fm[i] = 0;
        matrix.md[i] = 0;
        matrix.ud[i] = 0;
        matrix.ld[i] = 0;
        matrix.u[i] = 0;
        if (!params.shock_dep) {
            matrix.u[i+NR] = 0;
            matrix.u[i+2*NR] = 0;
        }
        matrix.w[i] = 0;
        if (rmin[i] >= aplanet) {
            if (need_set_flag) { 
                matrix.icol = i;
                if (params.shock_dep) {
                    matrix.w[matrix.icol] = 1;
                }
                need_set_flag = FALSE;
            }
        }
    }



    if (need_set_flag) {
        printf("Couldnt find planet in domain, a=%lg\n Setting planet to last grid point at r=%lg\n",aplanet,rc[NR-1]);
        matrix.icol = NR-1;
        if (params.shock_dep) {
             matrix.w[matrix.icol] = NR-1;
         }

    }
    matrix.md[NR-1] = 0;
    matrix.fm[NR-1] = 0;
    
    FILE *fout;

    if (params.shock_dep) {
        set_tau(aplanet);
    }
    else {
        set_weights(matrix.w,matrix.u,aplanet,matrix.icol,NR);
       /*
        fout=fopen("matrix_test.dat","w");
        for(i=0;i<NR;i++) {
            fprintf(fout,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",rc[i],y[i],matrix.w[i],matrix.w[i+NR],matrix.u[i+NR],matrix.u[i+NR*2],dep_func(rmin[i],aplanet,planet.xd,params.h*aplanet));
        }
      fclose(fout);
      */
    }

#pragma omp parallel for private(i,rm,rp,am,ap,bm,bp) shared(matrix,rc,rmin)
    for(i=1;i<NR-1;i++) {
        rm = rmin[i];
        rp = rmin[i+1];

        calc_coeffs(rm,rc[i]-rc[i-1],&am,&bm,aplanet,FALSE);
        calc_coeffs(rp,rc[i+1]-rc[i],&ap,&bp,aplanet,FALSE);

        matrix.md[i] = (ap-am - bm - bp)*dt/2.;
        matrix.ld[i-1] = (-am + bm)*dt/2.;
        matrix.ud[i] = (ap + bp)*dt/2.;

        if (params.planet_torque && params.shock_dep) {
            am  = calc_coeffs_nl(rm,aplanet,i);
            ap = calc_coeffs_nl(rp,aplanet,i+1);
            matrix.u[i] = (ap-am)*dt/2.;
        }

      }
    if (params.flux_bc) {
        calc_coeffs(rmin[NR-1],rc[NR-1]-rc[NR-2],&am,&bm,aplanet,FALSE);
        calc_coeffs(rmin[1],rc[1]-rc[0],&ap,&bp,aplanet,FALSE);
        matrix.md[0] = (ap-bp)*dt/2.;
        matrix.ud[0] = (ap+bp)*dt/2.;
        matrix.md[NR-1] = (-am-bm)*dt/2.;
        matrix.ld[NR-2] = (-am+bm)*dt/2.;
        matrix.fm[0] = -params.bc_mdot*dt;
        matrix.fm[NR-1] = params.bc_mdot*dt;

        if (params.planet_torque && params.shock_dep) {
            am=calc_coeffs_nl(rmin[NR-1],aplanet,NR-1);
            ap=calc_coeffs_nl(rmin[1],aplanet,1);
            matrix.u[0] = ap*dt/2.;
            matrix.u[NR-1] = -am*dt/2.;
        }
    }
    
    if (params.planet_torque && !params.shock_dep) {
        for(i=0;i<matrix.nsol*NR;i++) {
            matrix.u[i] *= dt/2.;
        }
    }
    
    matvec(matrix.ld,matrix.md,matrix.ud,y,matrix.fm,NR);
    
    double nlfac=0;
    double nlfac2 = 0;
    if (params.shock_dep) {
        nlfac = y[matrix.icol];
    }
    else {
#pragma omp parallel for private(i) reduction(+:nlfac,nlfac2)
        for(i=0;i<NR;i++) {
            nlfac += matrix.w[i]*y[i];
            nlfac2 += matrix.w[i+NR]*y[i];
        }
    }

//    printf("Total torques\nOuter = %.2e\nInner = %.2e\n",nlfac2,nlfac);

#pragma omp parallel for private(i) shared(nlfac,nlfac2,y,dr,matrix)
    for(i=0;i<NR;i++) {
        matrix.fm[i] += dr[i]*y[i];
        if (!params.shock_dep) {
            matrix.fm[i] += matrix.u[i+NR*2]*nlfac2 + matrix.u[i+NR]*nlfac;
            matrix.u[i+NR] *= -1;
            matrix.u[i+2*NR] *= -1;
        }
        else {
            matrix.fm[i] += matrix.u[i]*nlfac;
            matrix.u[i] *= -1;
        }
        matrix.md[i] = dr[i] - matrix.md[i];
        if (i<NR-1) {
            matrix.ld[i] *= -1;
            matrix.ud[i] *= -1;
        }
    }
/*
    if (params.flux_bc) {
        matrix.fm[0] = -params.bc_mdot;
        matrix.fm[NR-1] = params.bc_mdot;
        matrix.md[0] = 0;
        matrix.md[NR-1] = 0;
        matrix.ld[NR-2] = 0;
        matrix.ud[0] = 0;
    }
    else {
*/
    if (!params.flux_bc) {
        matrix.md[NR-1] = 1;
        matrix.ld[NR-2] = 0;

        if (params.shock_dep) {
            matrix.u[NR-1] = 0;
            matrix.u[0] = 0;
        }
        else {
            matrix.u[NR-1+NR*2] = 0;
            matrix.u[0 + NR] = 0;
        }
        matrix.fm[NR-1] = params.bc_lam[1];

        matrix.md[0] = 1;
        matrix.ud[0] = 0;
        matrix.fm[0] = params.bc_lam[0];
    }
    else {

        if (params.bc_lam[0] == 0) {
            matrix.md[0] = 1;
            matrix.ud[0] = 0;
            if (params.shock_dep) {
                matrix.u[0] = 0;
            }
            else {
                matrix.u[0+NR] = 0;
            }
            matrix.fm[0] = params.bc_lam[0];
        }
    }


    if (params.shock_dep) {
       trisolve_sm(matrix.ld,matrix.md,matrix.ud,matrix.fm,y,matrix.u,matrix.w,matrix.icol,NR);
    }
    else {
        for (i=0;i<NR;i++) {
            matrix.u[i] = matrix.fm[i];
        }
       trisolve_smn(matrix.ld,matrix.md,matrix.ud,matrix.u,y,matrix.w,matrix.nsol,NR);

    }

    return;
}

