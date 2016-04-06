#include "pdisk.h"
void matvec(double *ld, double *md, double *ud, double *a, double *b, int n) {
    int i;

    b[0] += md[0]*a[0] + ud[0]*a[1];


    for(i=1;i<n-1;i++) {
        b[i] += ld[i-1]*a[i-1] + md[i]*a[i] + ud[i]*a[i+1];

    }
    b[n-1] += ld[n-2]*a[n-2] + md[n-1] * a[n-1];

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
