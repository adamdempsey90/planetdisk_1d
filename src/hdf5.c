#include "pdisk.h"

void write_hdf5_double(double *data, hsize_t *dims, int ndims, hid_t group_path, char *name) {
  hid_t dspc_id, dset_id;
 
  dspc_id = H5Screate_simple(ndims,dims,NULL);
  dset_id = H5Dcreate(group_path,name,H5T_NATIVE_DOUBLE,dspc_id,H5P_DEFAULT);

    HDF5_INSERT_ERROR( H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data) );
	
    HDF5_INSERT_ERROR( H5Sclose(dspc_id));
   HDF5_INSERT_ERROR( H5Dclose(dset_id));



  return;


}
void write_hdf5_params(hid_t *params_id) {
  hid_t memtype,dspc_id, dset_id;
  hsize_t dims[1] = {1};




    param_t out_par;
    out_par.nr = params.nr;
    out_par.ri = params.ri;
    out_par.ro = params.ro;
    out_par.alpha= params.alpha;
    out_par.gamma= params.gamma;
    out_par.h = params.h;
    out_par.bc_lam_inner =  params.bc_lam[0];
    out_par.bc_lam_outer = params.bc_lam[1];
    out_par.bc_mdot = params.bc_mdot;
    out_par.flux_bc = params.flux_bc;
    out_par.dt= params.dt;
    out_par.cfl = params.cfl;
    out_par.nvisc = params.nvisc;
    out_par.nt = params.nt;
    out_par.release_time = params.release_time;
    out_par.start_ss = params.start_ss;
    out_par.read_initial_conditions = params.read_initial_conditions;
    out_par.planet_torque = params.planet_torque;
    out_par.explicit_stepper = params.explicit_stepper;
    out_par.move_planet = params.move_planet;
    out_par.move_planet_implicit = params.move_planet_implicit;
    out_par.gaussian = planet.gaussian;
    out_par.symmetric_torque = planet.symmetric_torque;
    out_par.nonlocal_torque = params.nonlocal_torque;
    out_par.shock_dep = params.shock_dep;
    out_par.hs_visc = params.hs_visc;
    out_par.one_sided = planet.onesided;
    out_par.a = planet.a;
    out_par.mp = planet.mp;
    out_par.G1 = planet.G1;
    out_par.beta = planet.beta;
    out_par.delta = planet.delta;
    out_par.c = planet.c;
    out_par.eps = planet.eps;
    out_par.xd = planet.xd;

    memtype = H5Tcreate (H5T_COMPOUND, sizeof (param_t));
     HDF5_INSERT_ERROR(H5Tinsert (memtype, "nr", HOFFSET (param_t, nr), H5T_NATIVE_INT));
 
    HDF5_INSERT_ERROR(H5Tinsert (memtype, "ri", HOFFSET (param_t, ri), H5T_NATIVE_DOUBLE));

    HDF5_INSERT_ERROR(H5Tinsert (memtype, "ro", HOFFSET (param_t, ro), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "alpha", HOFFSET (param_t, alpha), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "gamma", HOFFSET (param_t, gamma), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "h", HOFFSET (param_t, h), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "bc_lam_inner", HOFFSET (param_t, bc_lam_inner), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "bc_lam_outer", HOFFSET (param_t, bc_lam_outer), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "bc_mdot", HOFFSET (param_t, bc_mdot), H5T_NATIVE_DOUBLE));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "flux_bc", HOFFSET (param_t, flux_bc), H5T_NATIVE_INT));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "dt", HOFFSET (param_t, dt), H5T_NATIVE_DOUBLE));
   HDF5_INSERT_ERROR(H5Tinsert (memtype, "cfl", HOFFSET (param_t, cfl), H5T_NATIVE_DOUBLE));

   HDF5_INSERT_ERROR(H5Tinsert (memtype, "nvisc", HOFFSET (param_t, nvisc), H5T_NATIVE_DOUBLE));

      HDF5_INSERT_ERROR(H5Tinsert (memtype, "nt", HOFFSET (param_t, nt), H5T_NATIVE_INT));
 
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "release_time", HOFFSET (param_t, release_time), H5T_NATIVE_DOUBLE));
 
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "start_ss", HOFFSET (param_t, start_ss), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "read_initial_conditions", HOFFSET (param_t, read_initial_conditions), H5T_NATIVE_INT));
 
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "planet_torque", HOFFSET (param_t, planet_torque), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "explicit_stepper", HOFFSET (param_t, explicit_stepper), H5T_NATIVE_INT));
 
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "move_planet", HOFFSET (param_t, move_planet), H5T_NATIVE_INT));
 
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "move_planet_implicit", HOFFSET (param_t, move_planet_implicit), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "gaussian", HOFFSET (param_t, gaussian), H5T_NATIVE_INT));
 
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "symmetric_torque", HOFFSET (param_t, symmetric_torque), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "nonlocal_torque", HOFFSET (param_t, nonlocal_torque), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "shock_dep", HOFFSET (param_t, shock_dep), H5T_NATIVE_INT));
 HDF5_INSERT_ERROR(H5Tinsert (memtype, "hs_visc", HOFFSET (param_t, hs_visc), H5T_NATIVE_INT));

    HDF5_INSERT_ERROR(H5Tinsert (memtype, "one_sided", HOFFSET (param_t, one_sided), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "a", HOFFSET (param_t, a), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "mp", HOFFSET (param_t, mp), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "G1", HOFFSET (param_t, G1), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "beta", HOFFSET (param_t,beta), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "delta", HOFFSET (param_t, delta), H5T_NATIVE_DOUBLE));
      HDF5_INSERT_ERROR(H5Tinsert (memtype, "c",HOFFSET (param_t,c), H5T_NATIVE_DOUBLE));
        HDF5_INSERT_ERROR(H5Tinsert (memtype, "eps",HOFFSET (param_t,eps), H5T_NATIVE_DOUBLE));
        HDF5_INSERT_ERROR(H5Tinsert (memtype, "xd",HOFFSET (param_t,xd), H5T_NATIVE_DOUBLE));

  printf("%lg\n\n\n",out_par.xd); 


  dspc_id = H5Screate_simple(1,dims,NULL);
  dset_id = H5Dcreate(*params_id,"Parameters",memtype,dspc_id,H5P_DEFAULT);

  HDF5_INSERT_ERROR(H5Dwrite(dset_id,memtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,&out_par));

  HDF5_INSERT_ERROR(H5Sclose(dspc_id));
  HDF5_INSERT_ERROR(H5Dclose(dset_id));

    return;
}

void write_hdf5_file(void) {
  
  printf("Outputting Results to %s...\n",params.outputname);
  
  hid_t file_id = H5Fcreate(params.outputname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t root_id = H5Gcreate(file_id,"/Migration",0);
  hid_t mesh_id = H5Gcreate(root_id,"Mesh",0);
  hid_t solution_id = H5Gcreate(root_id,"Solution",0);
  hid_t matrix_id = H5Gcreate(root_id,"Matrix",0);
  hid_t params_id = H5Gcreate(root_id,"Parameters",0);
    hid_t steadystate_id = H5Gcreate(root_id,"SteadyState",0);


    hsize_t dims1[1] = {NR};
    hsize_t dims1_t[1]= {params.nt};
    hsize_t dims1_small[1] = {NR-1};
    hsize_t dims2[2] = {params.nt,NR};

// Write Mesh data
    write_hdf5_double(rc,dims1,1,mesh_id,"rc");
    write_hdf5_double(dr,dims1,1,mesh_id,"dr");
      write_hdf5_double(rmin,dims1,1,mesh_id,"rmin");
      write_hdf5_double(fld.lami,dims1,1,mesh_id,"lami");
      write_hdf5_double(fld.mdoti,dims1,1,mesh_id,"mdoti");
      write_hdf5_double(fld.nu_grid,dims1,1,mesh_id,"nu_grid");
    write_hdf5_double(tauc,dims1,1,mesh_id,"tauc");
    write_hdf5_double(taumin,dims1,1,mesh_id,"taumin");
// Write Matrix
   write_hdf5_double(matrix.md,dims1,1,matrix_id,"md");
    write_hdf5_double(matrix.ld,dims1_small,1,matrix_id,"ld");
    write_hdf5_double(matrix.ud,dims1_small,1,matrix_id,"ud");
   write_hdf5_double(matrix.fm,dims1,1,matrix_id,"fm");

// Write Solution
    write_hdf5_double(fld.sol,dims2,2,solution_id,"lam");
    write_hdf5_double(fld.torque,dims2,2,solution_id,"torque");
    write_hdf5_double(fld.sol_mdot,dims2,2,solution_id,"mdot");
    write_hdf5_double(fld.times,dims1_t,1,solution_id,"times");

     write_hdf5_double(fld.avals,dims1_t,1,solution_id,"avals");
     write_hdf5_double(fld.vs,dims1_t,1,solution_id,"vs");
  
// Steady State Solution
    write_hdf5_double(fld.sol_ss,dims2,2,steadystate_id,"lam_ss");
    write_hdf5_double(fld.lamp,dims2,2,steadystate_id,"lamp");
    write_hdf5_double(fld.lam0,dims2,2,steadystate_id,"lam0");
    write_hdf5_double(fld.mdot_ss,dims1_t,1,steadystate_id,"mdot_ss");
    write_hdf5_double(fld.vs_ss,dims1_t,1,steadystate_id,"vs_ss");
    write_hdf5_double(fld.efficiency,dims1_t,1,steadystate_id,"eff");
    write_hdf5_double(fld.ivals_ss,dims2,2,steadystate_id,"ivals_ss");
    write_hdf5_double(fld.kvals_ss,dims2,2,steadystate_id,"kvals_ss");

    write_hdf5_params(&params_id);

    HDF5_INSERT_ERROR(H5Gclose(mesh_id));
  HDF5_INSERT_ERROR(H5Gclose(matrix_id));
  HDF5_INSERT_ERROR(H5Gclose(solution_id));
 HDF5_INSERT_ERROR(H5Gclose(params_id));
 HDF5_INSERT_ERROR(H5Gclose(steadystate_id));

  HDF5_INSERT_ERROR(H5Gclose(root_id));
  HDF5_INSERT_ERROR(H5Fclose(file_id));


  return;


}
