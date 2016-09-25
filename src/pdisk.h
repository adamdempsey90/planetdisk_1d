#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include <string.h>
#include "defines.h"


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif


#define ZEROTORQUE
#define NR params.nr
#define TRUE 1
#define FALSE 0 
#define WRITE_MATRIX
//#define INITIAL_NOISE
#define MAXITERATIONS 300
#define MAXSTRLEN 300
#define HDF5_INSERT_ERROR(status)  if (status < 0) printf("HDF5 error at line %d\n",__LINE__);
#define MALLOC_SAFE(ptr) if (ptr == NULL) printf("Malloc error at line %d!\n",__LINE__);

#define BCMDOTIN 0
#define BCMDOTOUT 1
#define BCLAMIN 2
#define BCLAMOUT 3
#define BCGRADIN 4
#define BCGRADOUT 5
#define BCMIXEDIN 6
#define BCMIXEDOUT 7


typedef struct Parameters {
    double alpha,h,ri,ro,gamma, mach, mth, mvisc, tvisc; 
    double nu0;
    int nr,nt;
    double dt,tend;
    int nvisc;
    int logtime;
    int planet_torque, move_planet,move_planet_implicit;
    int read_initial_conditions;
    int explicit_stepper;
    int start_ss;
    int hs_visc;
    double bc_val[4];
    int bc_type[2];
    double bc_mdot;
    double cfl;
    int flux_bc;
    double release_time;
    char outputname[MAXSTRLEN];
    char torque_file[MAXSTRLEN];
    int nonlocal_torque;
    int shock_dep;
    int forced_torque;
} Parameters;


typedef struct param_t {
    int nr;
    double ri;
    double ro;
    double alpha;
    double gamma;
    double h;
    double bc_lam_inner;
    double bc_lam_outer;
    double bc_mdot;
    int flux_bc;
    double dt;
    double cfl;
    double nvisc;
    int nt;
    double release_time;
    int start_ss;
    int read_initial_conditions;
    int planet_torque;
    int explicit_stepper;
    int move_planet;
    int move_planet_implicit;
    int gaussian;
    int symmetric_torque;
    int nonlocal_torque;
    int shock_dep;
    int forced_torque;
    int hs_visc;
    double one_sided;
    double a;
    double mp;
    double G1;
    double beta;
    double delta;
    double c;
    double eps;
    double xd;
} param_t; 

typedef struct Planet {
    double a, omp, delta, G1, beta, mp,vs;
    double K,q,xd,wd;
    double rh, dep;
    double c,eps;
    double onesided;
    int gaussian;
    int symmetric_torque;
    int nonlocal_torque;
    int shock_dep;
    double T0;
} Planet;

typedef struct TridDiagMat {
    int size;
    int nsol;
    int icol;
    double *md, *ld, *ud, *fm;
    double *u, *w;
} TriDiagMat;

typedef struct Field {
    double *sol;
    double *sol_mdot;
    double *times;
    double *lami;
    double *mdoti;
    double *avals;
    double *vs; 
    double *torque;
    double *vs_ss;
    double *mdot_ss;
    double *ivals_ss;
    double *kvals_ss;
    double *sol_ss;
    double *lamp;
    double *lam0;
    double *efficiency;
    double *nu_grid;
    double *dTr;
    double *dep_func;
    double *grid_torque;
    double *grid_torquec;
} Field;


typedef struct SteadyStateField {
    double *lam;
    double *lamp;
    double mdot;
    double a;
    double vs;
    double *lam0;
    double mdot0;
    double *ivals;
    double *kvals;

} SteadyStateField;


double *rc, *rmin, *lam, *dr;
double *taumin, *tauc;
double *lrc, *lrmin;
double *mass, *ones, *mdot;
double *weights;
double dlr;
Parameters params;
Planet planet; 
TriDiagMat matrix; 
Field fld;
SteadyStateField fld_ss;




void advance_system(double dt, double *t, double tend);
void allocate_field( Field *tmpfld);
void free_field( Field *tmpfld);
void set_boundary(void);
double dTr_dep(double x,double a);
double dTr_ex(double x,double a);
double cs(double x);
double nu(double x);
double scaleH(double x);
double smoothing(double x, double x0, double w);
void write_hdf5_file(void);
void init_lam_from_file(void);
void init_lam(void);
void crank_nicholson_step(double dt, double aplanet, double *y);
void matvec(double *ld, double *md, double *ud, double *y, double *res, double a, double b,int n);
void matvec_full(double *ld, double *md, double *ud,double *u, double *w, double *y, double *res, double a, double b,int n,int nw);
double dotprod(double *v1, double *v2, int n);
void outer(double *v1, double *v2, double *res, int n);
void set_mdot(int planet_torque);
double calc_drift_speed(double a,double *y);
void move_planet(double dt, double *y, double *vs, double *a);
void read_param_file(char *parfile);
void set_grid(void);
void free_grid(void);
void set_matrix(void);
void free_matrix(void);
void set_uw(double *u, double *w, double a, int n);
void allocate_steady_state_field(SteadyStateField *tmpfld);
void free_steady_state_field(SteadyStateField *tmpfld);
void steadystate_config(SteadyStateField *tmpfld, double a);
void trisolve(double *ld, double *md, double *ud, double *d,double *sol,int n);
void trisolve_sm(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n);
void trisolve_sm2(double *ld, double *md, double *ud, double *d, double *sol, double *u, double *w,int n);
double calc_total_torque(double a, double *y);
void set_torque_nl(double a, double *y, double *res);
void set_torque(double a, double *y, double *res);
double calc_inner_torque(double a, double *y);
double calc_outer_torque(double a, double *y);
void zero_array(double *v, int n);
void steadystate_config_nl(SteadyStateField *tmpfld,double a);
void read_torque_file(Field *tmpfld,char *trq_file_name);
void set_planet_deposition(void);
double get_outer_bc_mdot(double);
double get_inner_bc_mdot(double);
