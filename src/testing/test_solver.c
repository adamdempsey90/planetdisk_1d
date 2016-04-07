#include "pdisk.h"
int main(int arc, char *argv[]) {
    int i;
    int n = atoi(argv[1]);
    double *md, *ld, *ud, *w, *u,*rhs,*sol,*sol1,*sol2;
    double *fm, *sol3, *sol4;
    double a, b;

    md = (double *)malloc(sizeof(double)*n);
    ld = (double *)malloc(sizeof(double)*(n-1));
    ud = (double *)malloc(sizeof(double)*(n-1));
    rhs = (double *)malloc(sizeof(double)*n);
    u = (double *)malloc(sizeof(double)*n*2);
    w  = (double *)malloc(sizeof(double)*n*2);
    fm  = (double *)malloc(sizeof(double)*n);
    sol = (double *)malloc(sizeof(double)*n);
    sol1 = (double *)malloc(sizeof(double)*n);
    sol2 = (double *)malloc(sizeof(double)*n);
    sol3 = (double *)malloc(sizeof(double)*n);
    sol4 = (double *)malloc(sizeof(double)*n);
    FILE *f = fopen("test_data.dat","r");
    fread(md,sizeof(double),n,f);
     fread(ld,sizeof(double),n-1,f);
     fread(ud,sizeof(double),n-1,f);
     fread(rhs,sizeof(double),n,f);
     fread(fm,sizeof(double),n,f);
    fread(w,sizeof(double),n*2,f);
    fread(u,sizeof(double),n*2,f);
    fread(&a,sizeof(double),1,f);
    fread(&b,sizeof(double),1,f);
     fclose(f);
    trisolve(ld,md,ud,rhs,sol,n);
    trisolve_sm( ld,  md,  ud,  rhs, sol1, u,  w,  n);
    trisolve_sm2( ld,  md,  ud,  rhs, sol2, u,  w,  n);
    for(i=0;i<n;i++) {
        sol3[i] = fm[i];
        sol4[i] = fm[i];
    }
    matvec(ld,md,ud,rhs,sol3,a,b,n);
    matvec_full(ld,md,ud,u,w,rhs,sol4,a,b,n,2);
    f = fopen("test_results.dat","w");
    fwrite(sol,sizeof(double),n,f);
    fwrite(sol1,sizeof(double),n,f);
    fwrite(sol2,sizeof(double),n,f);
    fwrite(sol3,sizeof(double),n,f);
    fwrite(sol4,sizeof(double),n,f);
    
    return 1;

}
