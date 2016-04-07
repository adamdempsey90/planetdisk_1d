#include "pdisk.h"
void set_weights( double *w,double *u,double a, int ia, int n) {
    int i;

    double hp = params.h * a;
    double xd = planet.xd;

    double leftr = (xd-.5)*hp;
    double rightr = (xd+.5)*hp;

    double normL,normR;

/* Lower then upper */
//#pragma omp parallel for private(i) shared(a,rc,dr,ia)
   for(i=0;i<ia;i++) {
        
        w[i] = dlr*dTr_kana(rc[i],a);  // w lower
        w[i+n] = 0;                       // w upper

//        normL = 1/hp ? in_region(rmin[i],a,leftr,rightr) : 0; 
//        normR = 1/hp ? in_region(rmin[i+1],a,leftr,rightr) : 0;
        normL = dep_func(rmin[i],a,xd,hp);
        normR = dep_func(rmin[i+1],a,xd,hp);
   
        if (i==0) {
            u[i+n] = -2*sqrt(rmin[i+1])*normR;
        }
        else {
            u[i+n] = -2*sqrt(rmin[i+1])*normR;  // u lower 
            u[i+n] -= -2*sqrt(rmin[i])*normL;
        }
        //printf("%lg\t%lg\n",(rmin[i]-a)/(params.h*a),dTr_nl(rmin[i+1],a,ia,TRUE));
        u[i+n*2] = 0;                     // u upper
    }


//#pragma omp parallel for private(i) shared(a,rc,dr)
    for(i=ia;i<n;i++) { 
        w[i+n] = dlr*dTr_kana(rc[i],a);  // w upper
        w[i] = 0;                           // w lower
//        normL = 1/hp ? in_region(rmin[i],a,leftr,rightr) : 0; 
//        normR = 1/hp ? in_region(rmin[i+1],a,leftr,rightr) : 0;
        normL = dep_func(rmin[i],a,xd,hp);
        normR = dep_func(rmin[i+1],a,xd,hp);
        if (i==n-1) {
            u[i+n*2] = 2*sqrt(rmin[i])*normL;
        }
        else {
            u[i+n*2] =  -2*sqrt(rmin[i+1])*normR;
            u[i+n*2] -= -2*sqrt(rmin[i])*normL;   // u upper
        }
        u[i+n] = 0;                           // u lower
    }

    if (params.flux_bc) {
        u[n] = 0; //-2*sqrt(rmin[1])*dTr_nl(rmin[1],i,a,TRUE);
        u[2*n+NR-1] =0;// 2*sqrt(rmin[NR-1])*dTr_nl(rmin[NR-1],i,a,TRUE);
    }
    return;
}

