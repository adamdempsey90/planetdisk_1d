#include "pdisk.h"
void trisolve_smn(double *ld, double *md, double *ud, double *d,double *sol, double *wn, int nsol, int n) {
    /* Thomas algorithm using Shermin-Morrison formula for an extra matrix A + u w^T, col, at index icol
     */
     /* only works for nsol = 3 for now */


    int i,j;
    double *cp, *bp, *dp;
    double *cp2;
    double *sol2;
     cp =  (double *)malloc(sizeof(double)*(n-1)*nsol);
    bp =  (double *)malloc(sizeof(double)*n*nsol);
    dp =  (double *)malloc(sizeof(double)*n*nsol);
    
    cp2 = (double *)malloc(sizeof(double)*nsol*(nsol-1));

     sol2 =  (double *)malloc(sizeof(double)*n*nsol);
  
     int indx,cindx;
    for(j=0;j<nsol;j++) {
        indx = j*n;
        cindx = j*(n-1);
        for(i=0;i<n-1;i++) { 
            cp[i+cindx] = 0;
            bp[i+indx] = 1;
            dp[i+indx] = 0;
        }
        bp[n-1 + indx] = 0;
        dp[n-1 + indx] = 0;
    

    
        cp[0+cindx] = ud[0]/md[0];


        for(i=1;i<n-1;i++) {
            cp[i+cindx] = ud[i]/(md[i]- ld[i-1]*cp[i-1+cindx]);
        }
        dp[0+indx] = d[0+indx]/md[0];
        for(i=1;i<n;i++) {
            dp[i+indx] = (d[i+indx] - ld[i-1]*dp[i-1+indx])/(md[i]-ld[i-1]*cp[i-1+cindx]);
        }

        sol2[n-1+indx] = dp[n-1+indx];
        for(i=n-2;i>=0;i--) {
            sol2[i+indx] = dp[i+indx] - cp[i+cindx]*sol2[i+1+indx];
        }
    }

/* Sol contains the solutions to the nsol rhs 
 * Now we combine them using the nsol weights
 */

    int k;
    double res;

    for(i=0;i<nsol-1;i++) {
        for(j=0;j<nsol;j++) {
            /* Compute w_i dot sol2_j */
            res = 0;
#ifdef OPENMP
#pragma omp parallel for private(i,j,k) reduction(+:res)
#endif
            for(k=0;k<n;k++) {
                res += wn[k + n*i] * sol2[k+j*n]; 
            }
            cp2[j+nsol*i] = res;
        }
    }
    double fac1,fac2,fac3;
    fac1=cp2[0]/(1+cp2[1]);
    fac2=cp2[2]/(1+cp2[1]);
    fac3=(-cp2[nsol]+fac1*cp2[nsol+1])/(1+cp2[nsol+2]-fac2*cp2[nsol+1]);
    
    
#ifdef OPENMP
#pragma omp parallel for private(i)
#endif
    for(i=0;i<n;i++) {
        sol[i] = sol2[i]- (fac1+fac2*fac3)*sol2[i+n] + fac3*sol2[i+n*2];

    }

    free(cp); free(bp); free(dp);
    free(cp2);
    free(sol2);
    return;
}

void trisolve_sm(double *ld, double *md, double *ud, double *d,double *sol,double *u, double *w, int icol, int n) {
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
#ifdef OPENMP
#pragma omp parallel for private(i) reduction(+:fac1,fac2) 
#endif
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
void trisolve(double *ld, double *md, double *ud, double *d,double *sol,int n) {
    int i;
    double *cp, *bp, *dp;
    MALLOC_SAFE(( cp =  (double *)malloc(sizeof(double)*(n-1))));
    MALLOC_SAFE(( bp =  (double *)malloc(sizeof(double)*n)));
    MALLOC_SAFE(( dp =  (double *)malloc(sizeof(double)*n)));

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

