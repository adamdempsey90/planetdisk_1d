#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define OPENMP


void trisolve(double *ld, double *md, double *ud, double *d,double *sol,int n);
void trisolve_sm(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n);
void trisolve_sm2(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n);
void matvec(double *ld, double *md, double *ud, double *y, double *res, double a, double b,int n);
void matvec_full(double *ld, double *md, double *ud,double *u, double *w, double *y, double *res, double a, double b,int n,int nw);
