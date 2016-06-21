#include "pdisk.h"

int locate(double *xx, int nx, double x) {
/* Bisection lookup routine. 
 * Given the ordered abcissas xx[0:nx-1], find the value of index j
 * such that x lies between xx[j] and xx[j+1].
 * Return *j=-1 or *j=nx if x lies outside the range of xx. 
*/
    
    int ju, jm, jl;
    int ascnd;

    jl = -1;
    ju = nx;

    ascnd = (xx[nx-1] >= xx[0]);

    while (ju-jl > 1) {
        jm = (ju + jl) >> 1;
        if (x >= xx[jm] == ascnd) {
            jl = jm;
        }
        else {
            ju = jm;
        }
    }
    if (x == xx[0]) {
        return 0;
    }
    else {
        if (x == xx[nx-1]) {
            return  nx-2;
        }
        else {
            return  jl;
        }
    }
}

void linear_interpolation(double *x,double *y,double *xd, double *yd,int nx,int nd) {
/* Linear interpolation of the data xd,yd with nd data points onto the nx x-values 
 * stored in x. Result stored in y. Both x and xd should be in order.
*/
    int  i,j;
  
    for(i=0;i<nx;i++) {
        j = locate(xd,nd,x[i]);
        if (j==-1) {
            printf("Grid point %lg below interpolation range %lg\n",x[i],xd[0]);
            y[i] = 0;
        }
        else {
            if (j==nd) {
                printf("Grid point %lg above interpolation range %lg\n",x[i],xd[nd-1]);
                y[i] = 0;
            }
            else {
                y[i] = yd[j] + (yd[j+1] - yd[j])*(x[i]-xd[j])/(xd[j+1]-xd[j]);
            }
        }

    }
    return;
}

void cubic_spline_interpolation(double *x,double *y,double *xd, double *yd,int nx,int nd) {
/* Cubic spline interpolation of the data xd,yd with nd data points onto the nx x-values 
 * stored in x. Result stored in y. Both x and xd should be in order.
*/
    int  i,j,ind;
    double *cs_m = (double *)malloc(sizeof(double)*nd);
    double *cs_ld = (double *)malloc(sizeof(double)*(nd-1));
    double *cs_ud = (double *)malloc(sizeof(double)*(nd-1));
    double *cs_rhs = (double *)malloc(sizeof(double)*nd);
    double *cs_k = (double *)malloc(sizeof(double)*nd);

    for(i=0;i<nd;i++) {
        if (i < nd-1) {
            cs_ld[i] = 0;
            cs_ud[i] = 0;
        }
        cs_m[i] = 0;
        cs_k[i] = 0;
        cs_rhs[i] = 0;
    }
    double dxL, dxR,dyL,dyR;
    dxR = xd[1]-xd[0];
    cs_m[0] = 2./dxR;
    cs_ud[0] = 1./dxR;
    cs_rhs[0] = 3*(yd[1]-yd[0])/(dxR*dxR);
    for(i=1;i<nd-1;i++) {
        dxL = xd[i] - xd[i-1];
        dxR = xd[i+1]-xd[i];
        dyL = yd[i]-yd[i-1];
        dyR = yd[i+1]-yd[i];
        cs_m[i] = 2./dxL + 2./dxR;
        cs_ld[i-1] = 1./dxL;
        cs_ud[i] = 1./dxR;
        cs_rhs[i] = 3*dyL/(dxL*dxL) + 3*dyR/(dxR*dxR);
    }
    dxL = xd[nd-1] - xd[nd-2];
    cs_m[nd-1] = 2./dxL;
    cs_ld[nd-2] = 1./dxL;
    cs_rhs[nd-1] = 3*(yd[nd-1]-yd[nd-2])/(dxL*dxL);

    trisolve(cs_ld,cs_m,cs_ud,cs_rhs,cs_k,nd);

         
    double t,a,b;
    for(i=0;i<nx;i++) {
        j = locate(xd,nd,x[i]);
        if (j==-1) {
            printf("Grid point %lg below interpolation range %lg\n",x[i],xd[0]);
            y[i] = 0;
        }
        else {
            if (j==nd) {
                printf("Grid point %lg above interpolation range %lg\n",x[i],xd[nd-1]);
                y[i] = 0;
            }
            else {
                t = (x[i] - xd[j])/(xd[j+1]-xd[j]);
                a = cs_k[j]*(xd[j+1]-xd[j]) - (yd[j+1]-yd[j]);
                b = -cs_k[j+1]*(xd[j+1]-xd[j]) + (yd[j+1]-yd[j]);
                y[i] = (1-t)*yd[j-1] + t*yd[j]  + t*(1-t)*(a*(1-t)+b*t);
            }
        }

    }
    free(cs_k);
    free(cs_m);
    free(cs_ld);
    free(cs_ud);
    free(cs_rhs);
    return;
}
void read_torque_file(Field *tmpfld,char *trq_file_name) {
    int i;
    int nd;
    double temp;

    printf("Reading torque from file %s\n",trq_file_name);
    FILE *f = fopen(trq_file_name,"r");

    fread(&temp,sizeof(double),1,f); // ny
    nd = (int)temp;
    fread(&temp,sizeof(double),1,f); // alpha
    fread(&temp,sizeof(double),1,f); // mdot
    fread(&temp,sizeof(double),1,f); // h
    fread(&temp,sizeof(double),1,f); // ymin
    fread(&temp,sizeof(double),1,f); // ymax

    double *r_torque = (double *)malloc(sizeof(double)*nd);
    double *dtr_torque = (double *)malloc(sizeof(double)*nd);

    fread(r_torque,sizeof(double),nd,f);
    fread(dtr_torque,sizeof(double),nd,f);
    fclose(f);

    linear_interpolation(rmin,tmpfld->grid_torque,r_torque,dtr_torque,NR,nd);
//    linear_interpolation(rc,tmpfld->grid_torquec,r_torque,dtr_torque,NR,nd);
    free(r_torque);
    free(dtr_torque);
    printf("Done\n");
    return;
}
