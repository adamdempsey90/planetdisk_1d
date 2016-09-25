#include "pdisk.h"
void trisolve(double *ld, double *md, double *ud, double *d,double *sol,int n) {
    int i;
    double *cp, *bp, *dp;
    cp =  (double *)malloc(sizeof(double)*(n-1));
    bp =  (double *)malloc(sizeof(double)*n);
    dp =  (double *)malloc(sizeof(double)*n);

    for(i=0;i<n-1;i++) {
        cp[i] = 0;
        bp[i] = 1;
        dp[i] = 0;
    }
    bp[n-1] = 0;
    dp[n-1] = 0;

    cp[0] = ud[0]/md[0];

    for(i=1;i<n-1;i++) {
        cp[i] = ud[i]/(md[i]- ld[i-1]*cp[i-1]);
    }
    dp[0] = d[0]/md[0];
    for(i=1;i<n;i++) {
        dp[i] = (d[i] - ld[i-1]*dp[i-1])/(md[i]-ld[i-1]*cp[i-1]);

    }

    sol[n-1] = dp[n-1];

    for(i=n-2;i>=0;i--) {
        sol[i] = dp[i] - cp[i]*sol[i+1];
    }

    free(cp); free(bp); free(dp);
    return;
}

void trisolve_sm(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n) {
    /* Thomas algorithm using Shermin-Morrison formula for 1 extra matrix A + u w^T
     */

    int i;
    double *x1;
    x1 = (double *)malloc(sizeof(double)*n);

    double fac1,fac2;

    trisolve(ld,md,ud,d,sol,n);
    trisolve(ld,md,ud,u,x1,n);

    fac1=0;fac2=0;
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:fac1,fac2) 
#endif
    for(i=0;i<n;i++) {
        fac1 += w[i]*sol[i];
        fac2 += w[i]*x1[i];
    }
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(fac1,fac2) 
#endif
    for(i=0;i<n;i++) {
        sol[i] += -fac1/(1+fac2)*x1[i];
    }


    free(x1);
    return;
}
void trisolve_sm2(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n) {
    /* Thomas algorithm using Shermin-Morrison formula for 2 extra matrices A + u1 w1^T + u2 w2^T
     */
    int i;
    double *x1;
    x1 = (double *)malloc(sizeof(double)*n);
    
    double fac1,fac2;

    trisolve_sm(ld,md,ud,d,sol,u,w,n);
    trisolve_sm(ld,md,ud,&u[n],x1,u,w,n);


    fac1=0;
    fac2=0;
#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:fac1,fac2) 
#endif
    for(i=0;i<n;i++) {
        fac1 += w[i+n]*sol[i];
        fac2 += w[i+n]*x1[i];
    }
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(fac1,fac2)
#endif
    for(i=0;i<n;i++) {
        
        sol[i] += -fac1/(1+fac2)*x1[i];

    }


    free(x1);
    return;
}
