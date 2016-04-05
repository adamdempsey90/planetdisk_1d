void test_matvec(void) {
   double c[10] = { 0.5,  0.32,  0.24,  0.81,  0.83,  0.25,  0.89,  0.66,  0.56,  0.};
   double md[10] = { 0.41,  0.12,  0.97,  0.52,  0.32,  0.28,  0.42,  0.58,  0.5 ,  0.51};  
   double ld[9] = {0.89,  0.2 ,  0.98,  0.9 ,  0.63,  0.19,  0.93,  0.79,  0.84};
   double ud[9] = {0.02,  0.19,  0.8 ,  0.4 ,  0.04,  0.11,  0.35,  0.79,  0.57};
   double a[10] = { 0.62,  0.98,  0.04,  0.85,  0.13,  0.87,  0.84,  0.77,  0.23,  0.43};
   double mvans[10]  = { 0.8314,  1.509 ,  0.9848,  1.8384,  1.1346,  1.5608,  1.4923,
               2.4229,  1.0314,  0.9004};
   double tans[10] = {1.90111459393,
    -13.9728491757,
    1.60394690777,
    1.84842666824,
    -4.30762459275,
    13.6213967067,
    -7.7289780403,
    4.42315829331,
    6.68673135109,
    -11.0134398724};

matvec(ld,md,ud,c,a,10);
   int i;
   double err = 0;
   printf("MatVec\n");
   for(i=0;i<10;i++) {
       printf("%.5lg\t%.5lg\n", mvans[i],a[i]);
       err += fabs( (mvans[i] - a[i])/mvans[i]);
    }
   printf("Matvec L1 error: %.2e\n", err);

    

   trisolve(ld,md,ud,c,a,10);

   err = 0;
   printf("TriSolve\n");
   for(i=0;i<10;i++) {
       printf("%.5lg\t%.5lg\n",tans[i],a[i]);
        err += fabs( (tans[i]-a[i])/tans[i]);
   }
    printf("TriSolve L1 error: %.2e\n", err);
   return;
}


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
