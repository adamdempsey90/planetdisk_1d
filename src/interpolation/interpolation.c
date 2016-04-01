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

