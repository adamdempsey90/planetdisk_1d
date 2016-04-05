#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void trisolve_sm(double *ld, double *md, double *ud, double *d,double *sol,double *u, double *w, int n) {
    /* Thomas algorithm using Shermin-Morrison formula for an extra matrix A + u w^T, col, at index icol
     */

    int i;
    double *cp, *bp, *dp;
    double *cp2, *bp2, *dp2;
    double *sol2;
     cp =  (double *)malloc(sizeof(double)*(n-1));
    bp =  (double *)malloc(sizeof(double)*n);
    dp =  (double *)malloc(sizeof(double)*n);

    cp2 =  (double *)malloc(sizeof(double)*(n-1));
    bp2 =  (double *)malloc(sizeof(double)*n);
    dp2 =  (double *)malloc(sizeof(double)*n);
     sol2 =  (double *)malloc(sizeof(double)*n);
  
     
    for(i=0;i<n-1;i++) {
        cp[i] = 0;
        bp[i] = 1;
        dp[i] = 0;
        cp2[i] = 0;
        bp2[i] = 1;
        dp2[i] = 0;
    }
    bp[n-1] = 0;
    dp[n-1] = 0;

    cp[0] = ud[0]/md[0];

    bp2[n-1] = 0;
    dp2[n-1] = 0;

    cp2[0] = ud[0]/md[0];
    for(i=1;i<n-1;i++) {
        cp[i] = ud[i]/(md[i]- ld[i-1]*cp[i-1]);
        cp2[i] = ud[i]/(md[i]- ld[i-1]*cp2[i-1]);
    }
    dp[0] = d[0]/md[0];
    dp2[0] = u[0]/md[0];
    for(i=1;i<n;i++) {
        dp[i] = (d[i] - ld[i-1]*dp[i-1])/(md[i]-ld[i-1]*cp[i-1]);
        dp2[i] = (u[i] - ld[i-1]*dp2[i-1])/(md[i]-ld[i-1]*cp2[i-1]);
    }

    sol[n-1] = dp[n-1];
    sol2[n-1] = dp2[n-1];
    for(i=n-2;i>=0;i--) {
        sol[i] = dp[i] - cp[i]*sol[i+1];
        sol2[i] = dp2[i] - cp2[i]*sol2[i+1];
    }

    double fac1 = 0;
    double fac2 = 0;
    for(i=0;i<n;i++) {
        fac1 += sol[i]*w[i];
        fac2 += sol2[i]*w[i]; 
    }
    fac1 /= (1+fac2);
    
    for(i=0;i<n;i++) {
        sol[i]  -= fac1*sol2[i];
    }

    free(cp); free(bp); free(dp);
    free(cp2); free(bp2); free(dp2);
    free(sol2);
    return;
}
int main(void) {
    int i; 
    int n = 10;
    double *md, *ld, *ud, *w, *u,*rhs,*sol,*ans;
    printf("allocating.\n");
    sol = (double *)malloc(sizeof(double)*n);
    rhs = (double *)malloc(sizeof(double)*n);
    md = (double *)malloc(sizeof(double)*n);
    u = (double *)malloc(sizeof(double)*n);
    w = (double *)malloc(sizeof(double)*n);
    ans = (double *)malloc(sizeof(double)*n);
    ld = (double *)malloc(sizeof(double)*(n-1));
    ud = (double *)malloc(sizeof(double)*(n-1));

    printf("reading.\n");
    FILE *f = fopen("test_data.dat","r");
     fread(md,sizeof(double),n,f);
    for(i=0;i<n;i++) {
        printf("%lg\n",md[i]);
    }

     fread(ld,sizeof(double),n-1,f);
     fread(ud,sizeof(double),n-1,f);
     fread(rhs,sizeof(double),n,f);
    fread(u,sizeof(double),n,f);
    fread(w,sizeof(double),n,f);
     fread(ans,sizeof(double),n,f);
     fclose(f);
    printf("solving.\n");
    trisolve_sm( ld,  md,  ud,  rhs, sol, u,  w,  n);
    double err;
    for(i=0;i<n;i++) {
        err = fabs( (sol[i]-ans[i])/ans[i]);
        printf("%d\t%lg\t%lg\t%e\n",i,ans[i],sol[i],err);
    }
    return 1;


}
