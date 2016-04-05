#include "pdisk.h"
double dTr(double x,double a) {

        double nlfac=0;
        double xd = planet.xd;
        double hp = params.h*a;
        double left_region = (xd-1)*hp/2;
        double right_region = (xd+1)*hp/2;
        if ( (fabs(x-a) > left_region) && (fabs(x-a) < right_region)) {
       //     printf("%lg inside %lg and %lg\n",(x-a)/hp,left_region,right_region);
            nlfac = 1./hp;
        }
        return nlfac;


}
/*
double dTr_kana(double x, double a) {
    
    double norm = planet.eps*.798*pow(planet.mp*params.mth,2)*pow(a,4);
    double dist = x-a;

    if (fabs(dist) < params.h*a*1.3) {
        return 0;
    }
    else {
        if (dist<0) norm *= -1;
        return norm/(x*pow(fabs(x-a),4));
    }


}
*/
