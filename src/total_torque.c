#include "pdisk.h"

double calc_inner_torque(double a, double *y) {
    int i;
    double res = 0;
//#ifdef _OPENMP
//#pragma omp parallel for private(i) reduction(+:res) 
//#endif
    for(i=0;rc[i]<a;i++) {
        if (rc[i] >= a) {
            printf("BAD\n");
        }
        res += dr[i]*dTr_ex(rc[i],a)*y[i];
    }
    return res;
}
double calc_outer_torque(double a, double *y) {
    int i;
    double res = 0;
//#ifdef _OPENMP
//#pragma omp parallel for private(i) reduction(+:res) 
//#endif
    for(i=NR-1;rc[i]>a;i--) {
        res += dr[i]*dTr_ex(rc[i],a)*y[i];
    }
    return res;
}

double calc_total_torque(double a, double *y) {
    int i;
    double res = 0;

#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:res) 
#endif
    for(i=0;i<NR;i++) {
        res += dr[i]*dTr_ex(rc[i],a)*y[i];  // Factor of 1/2 from torque normalization
    }

    return res;

}


void set_torque(double a, double *y, double *res) {


     set_torque_nl(a,y,res,FALSE);
    
    return;

}



