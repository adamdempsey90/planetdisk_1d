#include "pdisk.h"

double dTr_dep(double x,double a) {

#ifdef NONLOCAL
    return dTr_ex(x,a);
#else
    return dTr_ex(x,a);
#endif
}




