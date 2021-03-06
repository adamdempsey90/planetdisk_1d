#include "pdisk.h"

extern void init_linearwaves(int, int,int,double,double);

int main(int argc, char *argv[]) {
    
    int i,j;
    char parfile[MAXSTRLEN];
    if (argc < 2) {
        strcpy(parfile,"params.in");
    }
    else {
        strcpy(parfile,argv[1]);
    }
/*
#ifdef _OPENMP
    printf("Compilied code with _OPENMP, using %d threads\n",OMP_NUM_THREADS);
#endif
*/
    printf("Reading parameters from %s...\n",parfile);
    if (argc > 2) {
        read_param_file(parfile,argc-2,&argv[2]);
    }
    else {
        read_param_file(parfile,0,NULL);
    }
    
        printf("Setting up grid...\n");
    allocate_field(&fld);
    set_grid();
    printf("Initializing lambda...\n");
    if (params.read_initial_conditions == TRUE) {
        init_lam_from_file();
    }
    else {
        init_lam();
   }

    for(i=0;i<NR;i++) {
        fld.nu_grid[i] = nu(rc[i]);
    }
    allocate_steady_state_field(&fld_ss);

    
       
    set_matrix(); 
    printf("Set matrix\n");

    if (params.linear_torque) {
        init_linearwaves(NR, 1, 50, params.ri,params.ro);
    }

    double total_mass = 0;
    for(i=0;i<NR;i++) total_mass += dr[i]*lam[i];

    double dt = params.dt;
    int nt = params.nt;

    for(i=0;i<NR*nt;i++) {
        fld.sol[i] = 0;
        fld.sol_mdot[i] = 0;
    }
   
#ifdef NONLOCAL
    steadystate_config_nl(&fld_ss,planet.a);
#else
    steadystate_config(&fld_ss,planet.a);
#endif
    fld.vs_ss[0] = fld_ss.vs;
    fld.mdot_ss[0] = fld_ss.mdot;
    fld.efficiency[0] = fld_ss.mdot/fld_ss.mdot0;

    if (params.start_ss) {
        for(i=0;i<NR;i++) {
            lam[i] = fld_ss.lam[i];
            mdot[i] = fld_ss.mdot;
            fld.sol[i] = lam[i];
            fld.sol_mdot[i] = mdot[i];
        }
    }

    set_mdot(TRUE);
    set_torque(planet.a,lam,fld.torque);
    for(i=0;i<NR;i++) {
        fld.sol[i] = lam[i];
        fld.lami[i] = lam[i];
        fld.mdoti[i] = mdot[i];
        fld.sol_mdot[i] = mdot[i];
        fld.sol_ss[i] = fld_ss.lam[i];
        fld.lamp[i] = fld_ss.lamp[i];
        fld.lam0[i] = fld_ss.lam0[i];
        fld.ivals_ss[i] = fld_ss.ivals[i];
        fld.kvals_ss[i] = fld_ss.kvals[i];
        fld.dTr[i] = dTr_ex(rc[i],planet.a);
        fld.dep_func[i] = dep_func(rc[i],planet.a,planet.xd,planet.wd);
        //explicit_step_func(planet.a,lam,fld.mdL,fld.mdR);
    }


    planet.vs = calc_drift_speed(planet.a,lam);
    fld.avals[0] = planet.a; 
    fld.vs[0] = planet.vs;
    fld.times[0] = 0;
    for(i=0;i<nt-1;i++) {
        if (params.nvisc) {
            if (params.logtime) {
                fld.times[i+1] = pow(10,i * log10(params.tend*params.tvisc) /((double) (nt-1)));
            }
            else {
                fld.times[i+1] = (i+1)*params.tend*params.tvisc /((double) (nt-1));
            }
        }
        else {
            if (params.logtime) {
                fld.times[i+1] = pow(10,i * log10(params.tend) /((double) (nt-1)));
            }
            else {
                fld.times[i+1] = (i+1)*params.tend /((double) (nt-1));
            }
        }
    }

    printf("Viscous time is %.1e...\n", params.tvisc); 
    printf("Starting Integration...\n"); 
    //printf("dt = %.3e\n",set_dt());
    double t = 0;
    int planet_has_left = FALSE;

    if (params.one_step) {
        steady_state_step(planet.a,lam);
        set_mdot(params.planet_torque);     
        fld.avals[1] = planet.a;
        fld.vs[1] = planet.vs;
        set_torque(planet.a,lam,&fld.torque[1*NR]);
        fld.vs_ss[1] = fld_ss.vs;
        fld.mdot_ss[1] = fld_ss.mdot;
        fld.efficiency[i] = fld_ss.mdot/fld_ss.mdot0;

        for(j=0;j<NR;j++) {
            fld.sol[j + NR*1] = lam[j];
            fld.sol_mdot[j + NR*1] = mdot[j];
            fld.sol_ss[j + NR*1] = fld_ss.lam[j];
            fld.lamp[j+NR*1] = fld_ss.lamp[j];
            fld.lam0[j+NR*1] = fld_ss.lam0[j];
            fld.ivals_ss[j + NR*1] = fld_ss.ivals[j];
            fld.kvals_ss[j + NR*1] = fld_ss.kvals[j];
            fld.dTr[j + NR*1] = dTr_ex(rc[j],planet.a);
            fld.dep_func[j + NR*1] = dep_func(rc[j],planet.a,planet.xd,planet.wd);
        }
    }
    else {
        for(i=1;i<nt;i++) {
            //printf("t = %.2f\n", fld.times[i]);
            while (t < fld.times[i]) {
                advance_system(dt, &t, fld.times[i]);
                if ((planet.a > params.ro) || (planet.a < params.ri)) {
                    if (!planet_has_left) {
                        printf("\nPlanet has left the domain!\n");
                        planet_has_left= TRUE;
                    }
                }
            }
            printf("\r t = %.2e = %.2e tvisc\t%02d%% complete...", t,t/params.tvisc,(int)(100* i/((float)nt)));
            fflush(stdout);
            set_mdot(params.planet_torque);     
            fld.avals[i] = planet.a;
            fld.vs[i] = planet.vs;
            set_torque(planet.a,lam,&fld.torque[i*NR]);
/*
#ifdef NONLOCAL
        steadystate_config_nl(&fld_ss,planet.a);
#else
        steadystate_config(&fld_ss,planet.a);
#endif
*/
            fld.vs_ss[i] = fld_ss.vs;
            fld.mdot_ss[i] = fld_ss.mdot;
            fld.efficiency[i] = fld_ss.mdot/fld_ss.mdot0;

            for(j=0;j<NR;j++) {
                fld.sol[j + NR*i] = lam[j];
                fld.sol_mdot[j + NR*i] = mdot[j];
                fld.sol_ss[j + NR*i] = fld_ss.lam[j];
                fld.lamp[j+NR*i] = fld_ss.lamp[j];
                fld.lam0[j+NR*i] = fld_ss.lam0[j];
                fld.ivals_ss[j + NR*i] = fld_ss.ivals[j];
                fld.kvals_ss[j + NR*i] = fld_ss.kvals[j];
                fld.dTr[j + NR*i] = dTr_ex(rc[j],planet.a);
                fld.dep_func[j + NR*i] = dep_func(rc[j],planet.a,planet.xd,planet.wd);
            }
            //explicit_step_func(planet.a,lam,&fld.mdL[i*NR],&fld.mdR[i*NR]);
        
        }
    }
    printf("\n"); 
    printf("Outputting results...\n");
    write_hdf5_file();
    
    free_field(&fld);
    free_grid();
    free_matrix();
    free_steady_state_field(&fld_ss);
    return 1;
}

