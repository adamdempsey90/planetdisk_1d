#include "pdisk.h"
int main(int arc, char *argv[]) {
    int n = atoi(argv[1]);
    double *md, *ld, *ud, *w, *u,*rhs,*sol,*sol1,*sol2;

    md = (double *)malloc(sizeof(double)*n);
    ld = (double *)malloc(sizeof(double)*(n-1));
    ud = (double *)malloc(sizeof(double)*(n-1));
    rhs = (double *)malloc(sizeof(double)*n);
    u = (double *)malloc(sizeof(double)*n*2);
    w  = (double *)malloc(sizeof(double)*n*2);
    sol = (double *)malloc(sizeof(double)*n);
    sol1 = (double *)malloc(sizeof(double)*n);
    sol2 = (double *)malloc(sizeof(double)*n);

    FILE *f = fopen("test_data.dat","r");
    fread(md,sizeof(double),n,f);
     fread(ld,sizeof(double),n-1,f);
     fread(ud,sizeof(double),n-1,f);
     fread(rhs,sizeof(double),n,f);
    fread(w,sizeof(double),n*2,f);
    fread(u,sizeof(double),n*2,f);
     fclose(f);
    trisolve(ld,md,ud,rhs,sol,n);
    trisolve_sm( ld,  md,  ud,  rhs, sol1, u,  w,  n);
    trisolve_sm2( ld,  md,  ud,  rhs, sol2, u,  w,  n);
    f = fopen("test_results.dat","w");
    fwrite(sol,sizeof(double),n,f);
    fwrite(sol1,sizeof(double),n,f);
    fwrite(sol2,sizeof(double),n,f);
    
    return 1;

}
