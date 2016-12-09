#include "pdisk.h"
#include <liblinear.h>

Params linear_params;
Disk *disk;
Grid *grid;
double *work_arr;

void init_linearwaves(void) {
    int num_modes = params.mend-params.mstart+1;
    disk = (Disk *)malloc(sizeof(Disk));
    grid = (Grid *)malloc(sizeof(Grid));
    work_arr = (double *)malloc(sizeof(double)*NR);

    linear_params.np= 1 ;
    linear_params.rank= 0 ; 
    linear_params.n= params.nlinear ;
    linear_params.nrhs= 3;
    linear_params.nm= num_modes;
    linear_params.nphi= 2000 ;
    linear_params.mstart= params.mstart;
    linear_params.mend= params.mend ;
    linear_params.h= params.h;
    linear_params.mu= -.5;
    linear_params.delta= params.delta;
    linear_params.nuindx= params.gamma;
    linear_params.eta= 2./3;
    linear_params.alpha= params.alpha;
    linear_params.omf= planet.omp;
    linear_params.f= 0;
    linear_params.sig0= 1.0;
    linear_params.dlr= log(params.ro_lin/params.ri_lin)/(double)params.nlinear;
    linear_params.rmin= params.ri_lin;
    linear_params.rmax= params.ro_lin;
    linear_params.tol= 1e-8;
    linear_params.ieps= 0.0;
    linear_params.iso= TRUE;
    linear_params.pcorrect= FALSE;
    linear_params.zero_inner_bc= FALSE;
    linear_params.zero_outer_bc= FALSE;
    linear_params.simple_visc= FALSE;
    linear_params.indirect= FALSE;
    linear_params.fromfile= FALSE;
    linear_params.eps= planet.soft;
    linear_params.eps2= planet.soft*planet.soft;
    linear_params.a= planet.a;

    init_grid(params.mstart,params.mend,grid,linear_params,disk);
    int i;
    for(i=0;i<NR;i++) work_arr[i] = log( lam[i]/(2*M_PI*rc[i]));
    interpolate_onto_grid(lrc,work_arr,NR,grid->lr,disk->sigma,disk->dlsdlr,disk->d2lsdlr,linear_params.n);

    for(i=0;i<linear_params.n;i++) {
        disk->sigma[i] = exp(disk->sigma[i]);
        disk->d2lsdlr[i] = (disk->d2lsdlr[i]<0) ? fmax(disk->d2lsdlr[i],-1./(linear_params.h*linear_params.h)) : fmin(disk->d2lsdlr[i],1./(linear_params.h*linear_params.h));
    }

    FILE *f = fopen("test.dat","w");

    for(i=0;i<linear_params.n;i++) {
        fprintf(f,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",grid->r[i],disk->sigma[i],disk->dlsdlr[i],disk->d2lsdlr[i], grid->dppot[i],grid->drpot[i]);
    }
    fclose(f);


    return;
}
void zero_linearwaves(void) {
    int i;
    for(i=0;i<linear_params.n*linear_params.nm;i++) {
        grid->lamdep[i] = 0;
        grid->lamex[i] = 0;
        grid->fw[i] = 0;
        grid->drfw[i] = 0;
        grid->u[i] = 0;
        grid->v[i] = 0;
       // grid->s[i] = 0;
    }
    /*
    for(i=0;i<linear_params.nm;i++) {
        grid->TL[i] = 0;
        grid->TR[i] = 0;
    }
    */
    return;
}

void output_linear_torques(char *fname) {
    FILE *f = fopen(fname,"w");

    printf("Outputting to %s\n",fname);
    double n = (double)linear_params.n;
    double m = (double)grid->nm;
    double mstart = (double)grid->mvals[0];
    double mend = (double)grid->mvals[grid->nm-1];
    fwrite(&n,sizeof(double),1,f);
    fwrite(&m,sizeof(double),1,f);
    fwrite(&mstart,sizeof(double),1,f);
    fwrite(&mend,sizeof(double),1,f);
    fwrite(grid->TL,sizeof(double),grid->nm,f);
    fwrite(grid->TR,sizeof(double),grid->nm,f);
    fwrite(grid->r,sizeof(double),linear_params.n,f);
    fwrite(grid->lamex,sizeof(double),linear_params.n*grid->nm,f);
    fwrite(grid->lamdep,sizeof(double),linear_params.n*grid->nm,f);
    fwrite(grid->drfw,sizeof(double),linear_params.n*grid->nm,f);
    fwrite(grid->fw,sizeof(double),linear_params.n*grid->nm,f);
    fwrite(grid->dppot,sizeof(double),linear_params.n*grid->nm,f);
    fwrite(grid->u,sizeof(double complex),linear_params.n*grid->nm,f);
    fwrite(grid->v,sizeof(double complex),linear_params.n*grid->nm,f);
    fwrite(grid->s,sizeof(double complex),linear_params.n*grid->nm,f);
    fwrite(disk->sigma,sizeof(double),linear_params.n,f);

    fclose(f);


}

void calculate_linearwaves(double *dbar, double *lamex, double *lamdep, double *TL, double *TR) {
    int i;
    int num_modes = linear_params.mend - linear_params.mstart + 1;
    for(i=0;i<NR;i++) work_arr[i] = log( dbar[i]/(2*M_PI*rc[i]));
    interpolate_onto_grid(lrc,work_arr,NR,grid->lr,disk->sigma,disk->dlsdlr,disk->d2lsdlr,linear_params.n);
    for(i=0;i<linear_params.n;i++) {
        disk->sigma[i] = exp(disk->sigma[i]);
        disk->d2lsdlr[i] = (disk->d2lsdlr[i]<0) ? fmax(disk->d2lsdlr[i],-1./(linear_params.h*linear_params.h)) : fmin(disk->d2lsdlr[i],1./(linear_params.h*linear_params.h));
    }
    *TL = 0;
    *TR = 0;
    //zero_linearwaves();
    for(i=0;i<num_modes;i++) {
        linearwaves(i,grid,linear_params,disk,TRUE);
        *TL -= grid->TL[i];
        *TR += grid->TR[i];
    }
    //output_linear_torques("output_test.dat");
/*
    FILE *f = fopen("output_test.dat","w");
    fwrite(grid->r,sizeof(double),linear_params.n,f);
    fwrite(disk->sigma,sizeof(double),linear_params.n,f);
    fwrite(grid->TL,sizeof(double),num_modes,f);
    fwrite(grid->TR,sizeof(double),num_modes,f);

    fclose(f);
    int j;
    f = fopen("test_u_r.dat","w");
    for(j=0;j<num_modes;j++) {
        for(i=0;i<linear_params.n;i++) {
            fprintf(f,"%lg\t",creal(grid->u[i + linear_params.n*j]));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    f = fopen("test_u_i.dat","w");
    for(j=0;j<num_modes;j++) {
        for(i=0;i<linear_params.n;i++) {
            fprintf(f,"%lg\t",cimag(grid->u[i + linear_params.n*j]));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    f = fopen("test_lamex.dat","w");
    for(j=0;j<num_modes;j++) {
        for(i=0;i<linear_params.n;i++) {
            fprintf(f,"%lg\t",(grid->lamex[i + linear_params.n*j]));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    f = fopen("test_drfw.dat","w");
    for(j=0;j<num_modes;j++) {
        for(i=0;i<linear_params.n;i++) {
            fprintf(f,"%lg\t",(grid->drfw[i + linear_params.n*j]));
        }
        fprintf(f,"\n");
    }
    fclose(f);
    */
    return;
}
