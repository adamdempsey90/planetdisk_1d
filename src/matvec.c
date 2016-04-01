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

