#include "pdisk.h"
#include <ctype.h>

#define PRINT_DOUBLE(NAME,VAL) printf("\t%s = %lg\n",NAME,VAL)
#define PRINT_INT(NAME,VAL) printf("\t%s = %d\n",NAME,VAL)
#define PRINT_STR(NAME,VAL) printf("\t%s = %s\n",NAME,VAL)
#define FPRINT_DOUBLE(F,NAME,VAL) fprintf(f,"%s = %lg\n",NAME,VAL)
#define FPRINT_INT(F,NAME,VAL) fprintf(f,"%s = %d\n",NAME,VAL)
#define FPRINT_STR(F,NAME,VAL) fprintf(f,"%s = %s\n",NAME,VAL)

void set_var(char *name,int int_val, double double_val, int bool_val, char *str_val) {
    if (strcmp(name,"nr") == 0) {
        params.nr = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"ri") == 0) {	
        params.ri = double_val;
        PRINT_DOUBLE(name,double_val);

    }else if (strcmp(name,"ro") == 0) {	
        params.ro = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"alpha") == 0) {	
        params.alpha = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"gamma") == 0) {	
        params.gamma = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"h") == 0) {	
        params.h = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"inner_bc") == 0) {	
        if (strcmp(str_val,"MDOT") == 0) {
            params.bc_type[0] = BCMDOTIN;
        }
        else if (strcmp(str_val,"GRAD") == 0) {
            params.bc_type[0] = BCGRADIN;
        }
        else if (strcmp(str_val,"LAM") == 0) {
            params.bc_type[0] = BCLAMIN;
        }
        else {
            params.bc_type[0] = BCMIXEDIN;
        }
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"outer_bc") == 0) {	
        if (strcmp(str_val,"MDOT") == 0) {
            params.bc_type[1] = BCMDOTOUT;
        }
        else if (strcmp(str_val,"GRAD") == 0) {
            params.bc_type[1] = BCGRADOUT;
        }
        else if (strcmp(str_val,"LAM") == 0) {
            params.bc_type[1] = BCLAMOUT;
        }
        else {
            params.bc_type[1] = BCMIXEDOUT;
        }
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"bc_val_inner_0") == 0) {	
        params.bc_val[0] = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"bc_val_inner_1") == 0) {	
        params.bc_val[1] = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"bc_val_outer_0") == 0) {	
        params.bc_val[2] = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"bc_val_outer_1") == 0) {	
        params.bc_val[3] = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"dt") == 0) {	
        params.dt = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"cfl") == 0) {	
        params.cfl = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"tend") == 0) {	
        params.tend = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"nvisc") == 0) {	
        params.nvisc = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"logtime") == 0) {	
        params.logtime = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"nt") == 0) {	
        params.nt = int_val;
        PRINT_INT(name,int_val);

    }
    else if (strcmp(name,"release_time") == 0) {	
        params.release_time = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"start_ss") == 0) {	
        params.start_ss = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"read_initial_conditions") == 0) {	
        params.read_initial_conditions = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"mp") == 0) {	
        planet.mp = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"planet_torque") == 0) {	
        params.planet_torque = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"move_planet") == 0) {	
        params.move_planet = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"nonlocal_torque") == 0) {	
        planet.nonlocal_torque = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"shock_dep") == 0) {	
        planet.shock_dep = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"symmetric_torque") == 0) {	
        planet.symmetric_torque = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"a") == 0) {	
        planet.a = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"eps") == 0) {	
        planet.eps = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"outputname") == 0) {	
        sprintf(params.outputname,"outputs/%s",str_val);
        PRINT_STR(name,str_val);

    }

    return;
}

void read_param_file(char *fname) {
    FILE *f;

    char tok[20] = "\t :=>";

    char line[100],name[100],strval[100];
    char *data;
    double temp;
    int status;
    int int_val;
    int bool_val;
    char testbool;
    unsigned int i;

    f= fopen(fname,"r");

    while (fgets(line,100,f)) {
       // printf("%s\n",line);
        status = sscanf(line,"%s",name);

      //  printf("%s\n",name);
        if (name[0] != '#' && status == 1) {
        
             data = line + (int)strlen(name);
             sscanf(data + strspn(data,tok),"%lf",&temp);
             sscanf(data + strspn(data,tok),"%s",strval);
             //printf("%lf\t%s\n",temp,strval);
            int_val = (int)temp;
            testbool = toupper(strval[0]);
            if (testbool == 'Y') bool_val = TRUE;
            else bool_val = FALSE;
            
            for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
            
            set_var(name,int_val,temp,bool_val,strval);

        }
    }

    params.mach = 1/params.h;
    params.nu0 = params.alpha * params.h*params.h;
    params.mth = params.h*params.h*params.h;
    params.mvisc = sqrt(27.*M_PI/8  * params.alpha * params.mach);
    params.tvisc = params.ro*params.ro/nu(params.ro); 

    planet.q = planet.mp * params.mth;
    planet.rh = pow( planet.q/3.,1./3) * planet.a;
    planet.omp = pow(planet.a,-1.5);
    planet.dep = params.h*planet.a;
    planet.T0 = 2*M_PI*planet.a*planet.mp*planet.mp*params.mth/params.h;
    planet.vs = 0;
    planet.K = planet.q*planet.q*pow(params.h,-5.0)/params.alpha;

    set_planet_deposition();

    printf("Planet:\n\tK = %lg\n\txd = %lg\n\twd = %lg\n",planet.K,planet.xd,planet.wd);
    


    return;
}

