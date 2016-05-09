#include "pdisk.h"


int main(int argc, char *argv[]) {
    int nd = atoi(argv[1]);
    int nx = atoi(argv[2]);

    double *xd, *yd, *x, *y;

    xd = (double *)malloc(sizeof(double)*nd);
    yd = (double *)malloc(sizeof(double)*nd);
    x = (double *)malloc(sizeof(double)*nx);
    y = (double *)malloc(sizeof(double)*nx);

    FILE *f;
    int i;
    f = fopen("data.dat","r");
    for(i=0;i<nd;i++) {
        fscanf(f,"%lg\t%lg\n",&xd[i],&yd[i]);

    }
    fclose(f);

    f = fopen("xvals.dat","r");
    for(i=0;i<nx;i++) {
        fscanf(f,"%lg\n",&x[i]);

    }
    fclose(f);

    printf("Testing linear interpolation...\n");
    linear_interpolation(x,y,xd,yd,nx,nd);

    f = fopen("yvals_lint.dat","w");
    for(i=0;i<nx;i++) {
        fprintf(f,"%.16lg\n",y[i]);
    }
    fclose(f);

    for(i=0;i<nx;i++) y[i] = 0;

    printf("Testing cubic spline interpolation...\n");

    cubic_spline_interpolation(x,y,xd,yd,nx,nd);

    printf("Outputting cubic spline interpolation...\n");
    f = fopen("yvals_splint.dat","w");
    for(i=0;i<nx;i++) {
        fprintf(f,"%.16lg\n",y[i]);
    }
    fclose(f);
    
    return 0;

}
