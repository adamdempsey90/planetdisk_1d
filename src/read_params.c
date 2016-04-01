#include "pdisk.h"
void set_params(char *parfile) {
    read_input_file(parfile);
    
    params.mach = 1/params.h;
    params.nu0 = params.alpha * params.h*params.h;
    params.mth = params.h*params.h*params.h;
    params.mvisc = sqrt(27.*M_PI/8  * params.alpha * params.mach);
    params.tvisc = params.ro*params.ro/nu(params.ro); 
    planet.rh = pow( planet.mp * params.mth/3.,1./3) * planet.a;
    planet.omp = pow(planet.a,-1.5);
    planet.dep = params.h*planet.a;
    planet.T0 = 2*M_PI*planet.a*planet.mp*planet.mp*params.mth/params.h;
    planet.vs = 0;

    printf("Parameters:\n\tnr = %d\n\talpha = %.1e\n\th = %.2f\n\t(ri,ro) = (%lg,%lg)\n\tMach = %.1f\n\tm_th = %.2e\n\tm_visc = %.2e\n\tt_visc=%.2e\n",
            params.nr,
            params.alpha,
            params.h,
            params.ri,
            params.ro,
            params.mach,
            params.mth,
            params.mvisc,
            params.tvisc);

     printf("Planet properties:\n \
            \ta = %lg\n \
            \tmass = %lg mth = %.2e Mstar = %.2e mvisc\n \
            \tK = %.2e\n \
            \trh = %lg = %lg h\n",
            planet.a, 
            planet.mp,planet.mp*params.mth,planet.mp/params.mvisc,
            planet.mp*params.h/params.alpha,
            planet.rh, planet.rh/params.h);
     printf("\tplanet_torque = %s\n \
             \tmove_planet = %s\n \
             \tmove_planet_implicit = %s\n \
             \tgaussian = %s\n \
             \tread_initial_conditions = %s\n \
             \toutputname = %s\n",
             params.planet_torque ? "TRUE" : "FALSE", 
             params.move_planet ? "TRUE" : "FALSE", 
             params.move_planet_implicit ? "TRUE" : "FALSE", 
             planet.gaussian ? "TRUE" : "FALSE", 
             params.read_initial_conditions ? "TRUE" : "FALSE", 
             params.outputname);
    return;
}

void read_input_file(char *fname) {
    char garbage[100],tmpstr[MAXSTRLEN];
    char *gchar;
    int read_res;
    FILE *f;




    f = fopen(fname,"r");

    if (f==NULL) printf("\n\nERROR Can't Find Input File!\n\n");
  
	gchar=fgets(garbage,sizeof(garbage),f);	// Grid Parameters

    read_res=fscanf(f,"nr = %d \n",&NR);
    read_res=fscanf(f,"ri = %lg \n",&params.ri);
    read_res=fscanf(f,"ro = %lg \n",&params.ro);   
    gchar=fgets(garbage,sizeof(garbage),f); //	Disk Parameters:
    read_res=fscanf(f,"alpha = %lg \n",&params.alpha); 
    read_res=fscanf(f,"gamma = %lg \n",&params.gamma);
    read_res=fscanf(f,"h = %lg \n",&params.h);
   
    gchar=fgets(garbage,sizeof(garbage),f);	//Boundary Conditions:

    read_res=fscanf(f,"bc_lam_inner = %lg \n",&params.bc_lam[0]);
    read_res=fscanf(f,"bc_lam_outer = %lg \n",&params.bc_lam[1]);
    read_res=fscanf(f,"bc_mdot = %lg \n",&params.bc_mdot);
    read_res=fscanf(f,"flux_bc = %s \n",tmpstr);
    set_bool(tmpstr,&params.flux_bc);
    
    gchar=fgets(garbage,sizeof(garbage),f); //	Time Parameters:
    read_res=fscanf(f,"dt = %lg \n",&params.dt);
    read_res=fscanf(f,"cfl = %lg \n",&params.cfl);
    read_res=fscanf(f,"nvisc = %lg \n",&params.nvisc);
    read_res=fscanf(f,"nt = %d \n",&params.nt);
    read_res=fscanf(f,"release_time = %lg \n",&params.release_time);
    read_res=fscanf(f,"start_ss = %s \n",tmpstr);
    set_bool(tmpstr,&params.start_ss);
    read_res=fscanf(f,"read_initial_conditions = %s \n",tmpstr);
    set_bool(tmpstr,&params.read_initial_conditions);
  
    gchar=fgets(garbage,sizeof(garbage),f);	// Planet Properties:

    read_res=fscanf(f,"planet_torque = %s \n",tmpstr);
    set_bool(tmpstr,&params.planet_torque);
    read_res=fscanf(f,"explicit_stepper = %s \n",tmpstr);
    set_bool(tmpstr,&params.explicit_stepper);
    read_res=fscanf(f,"move_planet = %s \n",tmpstr);
    set_bool(tmpstr,&params.move_planet);
    read_res=fscanf(f,"move_planet_implicit = %s \n",tmpstr);
    set_bool(tmpstr,&params.move_planet_implicit);
    read_res=fscanf(f,"gaussian = %s \n",tmpstr);
    set_bool(tmpstr,&planet.gaussian);
    
    read_res=fscanf(f,"symmetric_torque = %s \n",tmpstr);
    set_bool(tmpstr,&planet.symmetric_torque);

    read_res=fscanf(f,"nonlocal_torque = %s \n",tmpstr);
    set_bool(tmpstr,&params.nonlocal_torque);
    read_res=fscanf(f,"shock_dep = %s \n",tmpstr);
    set_bool(tmpstr,&params.shock_dep);
    read_res=fscanf(f,"hs_visc = %s \n",tmpstr);
    set_bool(tmpstr,&params.hs_visc);

    read_res=fscanf(f,"one_sided = %lg \n",&planet.onesided);
    read_res=fscanf(f,"a  = %lg \n",&planet.a);
    read_res=fscanf(f,"mp = %lg \n",&planet.mp);
    read_res=fscanf(f,"G1 = %lg \n",&planet.G1);
    read_res=fscanf(f,"beta = %lg \n",&planet.beta);
    read_res=fscanf(f,"delta = %lg \n",&planet.delta);
    read_res=fscanf(f,"c = %lg \n",&planet.c);
    read_res=fscanf(f,"eps = %lg \n",&planet.eps);
    read_res=fscanf(f,"xd = %lg \n",&planet.xd);
    
    read_res=fscanf(f,"outputname = %s\n",tmpstr); 
    sprintf(params.outputname,"outputs/%s",tmpstr);
    fclose(f);
    printf("Outputting results to %s...\n",params.outputname);

/*
    rm_sub_string(fname,"params.in");
    char *dirname = fname;
    sprintf(outputname,"%s%s.hdf5",dirname,inputstr);
*/
// Leave me alone compiler
    read_res += 1; garbage[0] = gchar[0];


  return;

}

void set_bool(char *buff, int *val) {
    if ( (!strncmp(buff,"TRUE",MAXSTRLEN)) || (!strncmp(buff,"true",MAXSTRLEN)) || (!strncmp(buff,"True",MAXSTRLEN)) || (!strncmp(buff,"T",MAXSTRLEN)) || (!strncmp(buff,"T",MAXSTRLEN)) || (!strncmp(buff,"1",MAXSTRLEN))) {
        *val = TRUE;
            }
    else {
        *val = FALSE;
            }
    return;

}

