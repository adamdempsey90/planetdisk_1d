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


