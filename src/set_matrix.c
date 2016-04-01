#include "pdisk.h"

void set_matrix(void) {
    int i;
    matrix.size = NR;
    if (params.nonlocal_torque && !params.shock_dep) {
        matrix.nsol = 3;

        MALLOC_SAFE(( matrix.u =  (double *)malloc(sizeof(double)*NR*matrix.nsol))); 
        MALLOC_SAFE(( matrix.w =  (double *)malloc(sizeof(double)*NR*(matrix.nsol-1)))); 
    }
    else {
        matrix.nsol = 1;
        MALLOC_SAFE(( matrix.u =  (double *)malloc(sizeof(double)*NR))); 
        MALLOC_SAFE(( matrix.w =  (double *)malloc(sizeof(double)*NR))); 

    }
    
    matrix.icol = 0;
    MALLOC_SAFE(( matrix.ld = (double *)malloc(sizeof(double)*(NR-1))));
    MALLOC_SAFE(( matrix.md = (double *)malloc(sizeof(double)*NR)));
    MALLOC_SAFE(( matrix.ud =  (double *)malloc(sizeof(double)*(NR-1))));
    MALLOC_SAFE(( matrix.fm =  (double *)malloc(sizeof(double)*NR)));  

    printf("Initializing matrices...\n"); 
    for(i=0;i<NR-1;i++) {
        matrix.ld[i] = 0;
        matrix.md[i] = 0;
        matrix.ud[i] = 0;
        matrix.fm[i] = 0;
        matrix.u[i] = 0;
        matrix.w[i] = 0;
                
    }
    matrix.u[NR-1] = 0;
    matrix.w[NR-1] = 0;
    matrix.md[NR-1] = 1;
    matrix.fm[NR-1] = 0;
    return;
}

void free_matrix(void) {
    free(matrix.md); free(matrix.ld); free(matrix.ud); free(matrix.fm);
    free(matrix.u); free(matrix.w);
    return;
}
