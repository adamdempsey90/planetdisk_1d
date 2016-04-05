double calc_total_torque(double a, double *y) {
    int i;
    double res = 0;

#ifdef OPENMP
#pragma omp parallel for private(i) reduction(+:res) 
#endif
    for(i=0;i<NR;i++) {
        res += weights[i]*dTr(rc[i],a)*y[i];  // Factor of 1/2 from torque normalization
    }

    return res*dlr;

}

