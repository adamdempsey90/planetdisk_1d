#include "pdisk.h"
void matvec(double *ld, double *md, double *ud, double *y, double *res, double a, double b,int n) {
    /* res = a*res + b*mat.y, 
     * where mat is a tridiagonal matrix with 
     * md = main diagonal
     * ld = lower diagonal
     * ud = upper diagonal
     */
    int i;

    res[0] = a*res[0] + b*(md[0]*y[0] + ud[0]*y[1]);



#ifdef OPENMP
#pragma omp parallel for private(i) 
#endif
    for(i=1;i<n-1;i++) {
        res[i] = a*res[i] + b*( ld[i-1]*y[i-1] + md[i]*y[i] + ud[i]*y[i+1]);

    }
    res[n-1] = a*res[n-1] + b*( ld[n-2]*y[n-2] + md[n-1] * y[n-1]);

    return;
}

void matvec_full(double *ld, double *md, double *ud,double *u, double *w, double *y, double *res, double a, double b,int n,int nw) {
    /* res = a*res + b*(mat + u1 w1^T + u2 w2^T + ...).y, 
     * where mat is a tridiagonal matrix with 
     * md = main diagonal
     * ld = lower diagonal
     * ud = upper diagonal
     * u and w are vectors which add to mat via an outer product
     */
    int i,j,k;
    double temp;


    res[0] = a*res[0] + b* (md[0]*y[0] + ud[0]*y[1]);

    temp=0;
#ifdef OPENMP
#pragma omp parallel for private(j,k) collapse(2) reduction(+:temp)  
#endif
    for(j=0;j<nw;j++) {
        for(k=0;k<n;k++) {
           temp += b*u[0+n*j]*w[k + n*j]*y[k]; 
        }
    }
    res[0] += temp;

#ifdef OPENMP
#pragma omp parallel for private(i,j,k)  
#endif
    for(i=1;i<n-1;i++) {
        res[i] = a*res[i] + b* (ld[i-1]*y[i-1] + md[i]*y[i] + ud[i]*y[i+1]);
        for(j=0;j<nw;j++) {
            for(k=0;k<n;k++) {
                res[i] += b*u[i+n*j]*w[k + n*j]*y[k]; 
            }
        }
    }
    res[n-1] = a*res[n-1] + b*( ld[n-2]*y[n-2] + md[n-1] * y[n-1]);

    temp=0;
#ifdef OPENMP
#pragma omp parallel for private(j,k) collapse(2) reduction(+:temp)  
#endif
    for(j=0;j<nw;j++) {
        for(k=0;k<n;k++) {
            temp += b*u[n-1+n*j]*w[k + n*j]*y[k]; 
        }
    }
    res[n-1] += temp;
    return;
}
double dotprod(double *v1, double *v2, int n) {
    int i;

    double res=0;

#ifdef OPENMP
#pragma omp parallel for private(i) reduction(+:res)
#endif
    for(i=0;i<n;i++) {
        res += v1[i]*v2[i];
    }
    return res;


}
void outer(double *v1, double *v2, double *res, int n) {
    int i,j;

#ifdef OPENMP
#pragma omp parallel for private(i,j) collapse(2)
#endif
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            res[j+i*n] = v1[i]*v2[j];
        }
    }

    return;
}
void zero_array(double *v, int n) {
    int i;
    for(i=0;i<n;i++) {
        v[i] = 0;
    }
    return;
}
