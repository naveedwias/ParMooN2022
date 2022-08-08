/** ****************************************************************************
 *
 * @name  snaps_pod_rom_test
 * @brief A test program to test snapshot Collector, POD-basis computation
 *
 * @autor Baptiste Moreau
 * @date  18.03.2019
*******************************************************************************/

#include <AuxParam2D.h>
#include <ConvDiff.h>
#include <Database.h>
#include <Domain.h>
#include <MainUtilities.h>
#include <ParMooN_repository_info.h>
#include <SnapshotsCollector.h>
#include <TCD_POD.h>
#include <TCD_ROM.h>
#include <TimeConvectionDiffusion.h>
#include <TimeDiscretizations.h>
#include <TimeDiscRout.h>
#include "TimeNavierStokes.h"
#include <TNSE_POD.h>
#include <TNSE_ROM.h>
#include "ParMooN.h"

#include <cmath>
#include <list>



/** ***************************************************************************/
ParameterDatabase set_database_tcd()
{
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Example2D::default_example_database());
  db.merge(LocalAssembling2D::default_local_assembling_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  db["example"]         = 0;
  db["problem_type"]    = 2;
  db["reynolds_number"] = 1;

  db["time_discretization"] = "backward_euler";
  db["time_start"]          = 0.;
  db["time_step_length"]    = 0.01;
  db["time_end"]            = 0.1;

  db["write_snaps"]       = true;
  db["compute_POD_basis"] = true;
  db["ROM_method"]        = true;

  // declaration of databases
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles"});
  db.add("refinement_n_initial_steps", (size_t) 2,"");

  // declaration of parameters used for snapshots, POD and ROM
  db.add("snaps_directory", "snapshots_tmp", "");
  db.add("snaps_basename", "snaps_tmp", "");
  db.add("snaps_time_derivative", false, "");
  db.add("snaps_time_derivative_projection", false, "");
  db.add("steps_per_snap", (size_t) 5, "");
  db.add("pod_directory", "pod_basis_tmp", "");
  db.add("pod_basename", "pod_tmp", "");
  db.add("pod_rank", (size_t) 0, "");
  db.add("pod_fluctuations_only", true, "");
  db.add("pod_inner_product", "L2", "");
  db.add("rom_init_regularized", true, "");
  db.add("differential_filter_width", 0.00625, "");

  // test direct solver with Galerkin
  Output::print("\n\nTesting Galerkin\n");
  db.add("solver_type", "direct", "", {"direct", "petsc"});
  db["space_discretization_type"] = "galerkin";

  return db;
}

/** ***************************************************************************/
void test_mean_tcd(std::vector<double> mean_snaps)
{
  const int nb_mean = 25;
  double    eps     = 1e-10;

  double mean_exp[nb_mean] = {1.63237083058, -0.0188781592017, 0.00270594410213,
                              0.0204667777831, -1.63524406776, -0.0188781592017,
                              1.63237083058, 0.0204667777831, -1.62796917941,
                              0., 0., 0., 0., 0., 0., -8.1643119943e-17, 0.,
                              -9.99839855107e-33, 1.99967971022e-32,
                              8.1643119943e-17, 8.1643119943e-17,
                              -9.99839855107e-33, 0., 0., -8.1643119943e-17};


  for(int i=0 ; i<nb_mean ; i++)
  {
      if( std::abs(mean_snaps[i] - mean_exp[i]) > eps )
      {
        ErrThrow("test snapshot with tcd2d: the average of snapshots is not "
                 "correct ", std::setprecision(12), mean_snaps[i],
                 " != ", std::setprecision(12), mean_exp[i]);
      }
  }
}

/** ***************************************************************************/
void test_eigs_tcd(const double* eigs_pod)
{
  const int nb_eigs = 2;
  double    eps     = 1e-10;

  double eigs_exp[nb_eigs] = {0.0276895090976, 1.13837961698e-06};

  for(int i=0 ; i<nb_eigs ; i++)
  {
      if( std::abs(eigs_pod[i] - eigs_exp[i]) > eps )
      {
        ErrThrow("test POD with tcd2d: the eigenvalue is not correct ",
                 std::setprecision(12), eigs_pod[i],
                 " != ", std::setprecision(12), eigs_exp[i]);
      }
  }
}

/** ***************************************************************************/
void test_sols_tcd(TCD_ROM<2>& tcd_rom)
{
  const int nb_red = 2;
  const int nb_val = 25;
  double    eps    = 1e-10;

  std::vector<double> red_sol = tcd_rom.get_sol_r();
  std::vector<double> full_sol(nb_val);
  tcd_rom.ROM::get_full_solution(red_sol, &full_sol[0]);

  double r_sol_ref[nb_red] = {0.196019680498, -0.000840001193343};

  double f_sol_ref[nb_val] = {2.22067477936, -0.0320697148294,
                              0.00500075392263, 0.035291704813,
                              -2.22569248147, -0.0320697148294,
                              2.22067477936, 0.0352917048131,
                              -2.21240365636, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0};

  for(int i=0 ; i<nb_red ; i++)
  {
      if( std::abs(red_sol[i] - r_sol_ref[i]) > eps )
      {
        ErrThrow("test ROM with tcd2d: the reduced solution is not correct ",
                 std::setprecision(12), red_sol[i],
                 " != ", std::setprecision(12), r_sol_ref[i]);
      }
  }

  for(int i=0 ; i<nb_val ; i++)
  {
      if( std::abs(full_sol[i] - f_sol_ref[i]) > eps )
      {
        ErrThrow("test ROM with tcd2d: the full solution is not correct ",
                 std::setprecision(12), full_sol[i],
                 " != ", std::setprecision(12), f_sol_ref[i]);
      }
  }
}



/** ***************************************************************************/
ParameterDatabase set_database_tnse(std::string inner_prod)
{
  const std::string path_to_repo = parmoon::source_directory;
  const std::string path_to_meshes = path_to_repo
                                     + "/data/mesh/flow_around_cylinder2D/";
  const std::string filename_prm = path_to_meshes + "flow_around_cylinder.PRM";
  const std::string filename_geo = path_to_meshes
                                 + "flow_around_cylinder_quad_structured.mesh";

    ParameterDatabase db = ParameterDatabase::parmoon_default_database();
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());

    db["example"]         = 3;
    db["problem_type"]    = 6;
    db["reynolds_number"] = 1000;

    db["time_discretization"] = "backward_euler";
    db["time_start"]          = 0.;
    db["time_step_length"]    = 0.01;
    db["time_end"]            = 0.1;

    db["write_snaps"]       = true;
    db["compute_POD_basis"] = true;
    db["ROM_method"]        = true;

    // declaration of databases
    db.add("boundary_file", filename_prm, "");
    db.add("geo_file", filename_geo, "");
    db.add("refinement_n_initial_steps", (size_t) 0,"");

    // declaration of parameters used for snapshots, POD and ROM
    db.add("snaps_directory", "snapshots_tmp", "");
    db.add("snaps_basename", "snaps_tmp", "");
    db.add("snaps_time_derivative", false, "");
    db.add("snaps_time_derivative_projection", false, "");
    db.add("steps_per_snap", (size_t) 5, "");
    db.add("pod_directory", "pod_basis_tmp", "");
    db.add("pod_basename", "pod_tmp", "");
    db.add("pod_rank", (size_t) 0, "");
    db.add("pod_fluctuations_only", true, "");
    db.add("pod_inner_product", inner_prod, "");
    db.add("rom_init_regularized", false, "");

    // test direct solver with Galerkin
    Output::print("\n\nTesting Galerkin\n");
    db.add("solver_type", "direct", "", {"direct", "petsc"});
    db["space_discretization_type"] = "galerkin";

    return db;
}

/** ***************************************************************************/
void test_mean_tnse(void)//std::vector<double> mean_snaps)
{
  const int nb_mean = 25;
  double    eps     = 1e-10;

double mean_snaps[nb_mean] = {};
  double mean_exp[nb_mean] = {};


  for(int i=0 ; i<nb_mean ; i++)
  {
      if( std::abs(mean_snaps[i] - mean_exp[i]) > eps )
      {
        ErrThrow("test snapshot with tnse2d: the average of snapshots is not "
                 "correct ", std::setprecision(12), mean_snaps[i],
                 " != ", std::setprecision(12), mean_exp[i]);
      }
  }
}

/** ***************************************************************************/
void test_eigs_tnse(void)//const double* eigs_pod)
{
  const int nb_eigs = 2;
  double    eps     = 1e-10;
double eigs_pod[nb_eigs] = {};
  double eigs_exp[nb_eigs] = {};

  for(int i=0 ; i<nb_eigs ; i++)
  {
      if( std::abs(eigs_pod[i] - eigs_exp[i]) > eps )
      {
        ErrThrow("test POD with tnse2d: the eigenvalue is not correct ",
                 std::setprecision(12), eigs_pod[i],
                 " != ", std::setprecision(12), eigs_exp[i]);
      }
  }
}

/** ***************************************************************************/
void test_sols_tnse(TNSE_ROM<2>& tnse_rom)
{
  const int nb_red = 2;
  const int nb_val = 25;
  double    eps    = 1e-10;

  std::vector<double> red_sol ;//= tnse_rom.get_sol_r();
  std::vector<double> full_sol(nb_val);
//   tnse_rom.ROM::get_full_solution(red_sol, &full_sol[0]);

  double r_sol_ref[nb_red] = {};

  double f_sol_ref[nb_val] = {};

  for(int i=0 ; i<nb_red ; i++)
  {
      if( std::abs(red_sol[i] - r_sol_ref[i]) > eps )
      {
        ErrThrow("test ROM with tnsed: the reduced solution is not correct ",
                 std::setprecision(12), red_sol[i],
                 " != ", std::setprecision(12), r_sol_ref[i]);
      }
  }

  for(int i=0 ; i<nb_val ; i++)
  {
      if( std::abs(full_sol[i] - f_sol_ref[i]) > eps )
      {
        ErrThrow("test ROM with tnse2d: the full solution is not correct ",
                 std::setprecision(12), full_sol[i],
                 " != ", std::setprecision(12), f_sol_ref[i]);
      }
  }
}


/** ***************************************************************************/
int remove_tmp_test_files(ParameterDatabase db, std::vector<std::string> variable_suffix={""})
{
  int length = variable_suffix.size();
  std::vector<int> status(length, 0);

  for(unsigned int i = 0 ; i < (unsigned)length ; ++i)
  {
    std::string snaps_dir_name = db["snaps_directory"];
    std::string pod_dir_name   = db["pod_directory"];
    std::string snaps_name    = snaps_dir_name + "/"
                              + db["snaps_basename"].get<std::string>()
                              + variable_suffix[i];
    std::string pod_basename  = pod_dir_name + "/"
                              + db["pod_basename"].get<std::string>()
                              + variable_suffix[i] + "."
                              + db["pod_inner_product"].get<std::string>() + ".";
    if(db["pod_fluctuations_only"])
    {
      pod_basename += "fluc.";
    }
    std::string pod_name_pod  = pod_basename + "pod";
    std::string pod_name_mean = pod_basename + "mean";
    std::string pod_name_eigs = pod_basename + "eigs";

    status[i] |= (std::remove(snaps_name.c_str()) == 0) ? 0 : 1;
    if(i+1==(unsigned)length)
      status[i] |= (std::remove(snaps_dir_name.c_str()) == 0) ? 0 : 1;

    status[i] |= (std::remove(pod_name_pod.c_str()) == 0) ? 0 : 1;
    status[i] |= (std::remove(pod_name_mean.c_str()) == 0) ? 0 : 1;
    status[i] |= (std::remove(pod_name_eigs.c_str()) == 0) ? 0 : 1;
    if(i+1==(unsigned)length)
      status[i] |= (std::remove(pod_dir_name.c_str()) == 0) ? 0 : 1;

  }
  return std::accumulate(status.begin(),status.end(),0);;
}

/** ***************************************************************************/
void test_pod_rom_tcd(void)
{
//   TDatabase Database;
  ParameterDatabase db = set_database_tcd();
  TDatabase::ParamDB->ANSATZ_ORDER=1;
  TDomain domain(db);

  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);

  TimeConvectionDiffusion<2> tcd(domain, db);

  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();

  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();

  tcd.assemble_initial_time();

  //----------------------------------------------------------------------------
  // solve the PDE using FEM and acquire snapshots
  //----------------------------------------------------------------------------
  if(db["write_snaps"])
  {
    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps(db);

    // store initial condition as snapshot
    snaps.write_data(tcd.get_solution().block(0),
                     tcd.get_solution().length(0));

    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      SetTimeDiscParameters(1);

      tss.current_time_ += tss.get_step_length();;
      TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();

      Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
      tcd.assemble();
      tcd.solve();

      // write the snapshots
      if(db["write_snaps"])
      {
        snaps.write_data(tcd.get_solution().block(0),
                         tcd.get_solution().length(0),
                         tss.current_step_);
      }
    }
  }

  //----------------------------------------------------------------------------
  // compute the POD basis from snapshots
  //----------------------------------------------------------------------------
  auto collections            = domain.get_grid_collections();
  TCollection& cellCollection = *collections.front();
  TCD_POD<2> tcd_pod(cellCollection, db);

  if(db["compute_POD_basis"])
  {
    tcd_pod.compute_pod_basis();
  }

  // ===========================================================================
  // Solve the PDE using ROM instead of FEM
  // ===========================================================================
  double start_time = db["time_start"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  TCD_ROM<2> tcd_rom(domain, db);

  if( db["ROM_method"] )
  {
    tss.current_time_ = start_time;
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble matrices and right hand side at start time
    tcd_rom.assemble_initial_time();
//     tcd_rom.output();

    // time iteration
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;
      tss.set_time_disc_parameters();
      tss.current_time_ += tss.get_step_length();
      double tau = db["time_step_length"];

      TDatabase::TimeDB->CURRENTTIME += tau;
      Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);
      SetTimeDiscParameters(1);

      tcd_rom.assemble();
      tcd_rom.solve();

//      tcd_rom.output();
    }
  } // if ROM

  //----------------------------------------------------------------------------
  // test snapshots average and eigenvalues from POD
  // test reduced and full solutions from ROM
  // delete temporary files used for this test
  //----------------------------------------------------------------------------
  test_mean_tcd(tcd_pod.get_pod_c()->get_snaps_avr());
  test_eigs_tcd(tcd_pod.get_pod_c()->get_pod_eigs());
  test_sols_tcd(tcd_rom);

  if(remove_tmp_test_files(db) != 0)
  {
    ErrThrow("snaps_pod_rom_test: TCD Error deleting temporary files.");
  }
//  TDatabase::destroy();
}


/** ***************************************************************************/
void test_pod_rom_tnse(void)
{
  unsigned int dim = 2;
  std::string inner_prod = "L2"; //"H1";
  ParameterDatabase db = set_database_tnse(inner_prod);
  TDatabase::ParamDB->VELOCITY_SPACE = 2;
//   TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = 4;
  TDomain domain(db);

  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);

  // set some parameters for time stepping
  SetTimeDiscParameters(0);

  TimeNavierStokes<2> tnse(domain, db);

  TimeDiscretization& tss = tnse.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.set_time_disc_parameters();

  TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();

  tnse.assemble_initial_time();

  //----------------------------------------------------------------------------
  // solve the PDE using FEM and acquire snapshots
  //----------------------------------------------------------------------------
  if(db["write_snaps"])
  {
    // initialize snapshot writer (only needed if parmoon_db["write_snaps"])
    // the solutions are stored in a special collector which can be later read
    // by the computation of POD the basis
    SnapshotsCollector snaps_u(db, "u");
    SnapshotsCollector snaps_p(db, "p");

    // store initial condition as snapshot
    snaps_u.write_data(tnse.get_solution().block(0),
                       tnse.get_solution().length(0)*dim);
    snaps_p.write_data(tnse.get_solution().block(dim),
                       tnse.get_solution().length(dim));

    
//     LoopInfo<double> loop_info_time("time loop");
//     loop_info_time.print_time_every_step = true;
//     loop_info_time.verbosity_threshold = 1;
    int linear_iteration=0;

    TDatabase::TimeDB->CURRENTTIME = tss.get_start_time();
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;

      // set the time parameters
      tss.set_time_disc_parameters();
      // tau may change depending on the time discretization (adaptive time)
      double tau = tss.get_step_length();
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tau;
//       Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);

      tnse.assemble_matrices_rhs(0);

//       LoopInfo<Residuals> loop_info("nonlinear");
//       loop_info.print_time_every_step = true;
//       loop_info.verbosity_threshold = 1;
      for(unsigned int i=0;; i++)
      {
        if(tnse.stop_it(i))
        {
//           loop_info.finish(i,tnse2d.get_residuals());
          linear_iteration +=i;
//           loop_info_time.print(linear_iteration, tnse2d.get_full_residual());
          break;
        }
//         else
//           loop_info.print(i, tnse2d.get_residuals());

        tnse.solve();
        tnse.assemble_matrices_rhs(i+1);
      }

      if(db["write_snaps"])
      {
        // solution_m1 and solution_m2 were already updated
        snaps_u.write_data(tnse.get_solution().block(0),
                           tnse.get_solution().length(0)*dim,
                           tss.current_step_,
                           tnse.get_solution_m2().block(0),
                           tau);
        snaps_p.write_data(tnse.get_solution().block(dim),
                           tnse.get_solution().length(dim),
                           tss.current_step_,
                           tnse.get_solution_m2().block(dim),
                           tau);
      }
    }
  } // if(db["write_snaps"])

  //----------------------------------------------------------------------------
  // compute the POD basis from snapshots
  //----------------------------------------------------------------------------
  auto collections            = domain.get_grid_collections();
  TCollection& cellCollection = *collections.front();
  TNSE_POD<2> tnse_pod(cellCollection, db);

  if(db["compute_POD_basis"])
  {
    tnse_pod.compute_pod_basis(db);
  }

  // ===========================================================================
  // Solve the PDE using ROM instead of FEM
  // ===========================================================================
  double start_time = db["time_start"];
  TDatabase::TimeDB->CURRENTTIME = start_time;
  TNSE_ROM<2> tnse_rom(domain, db);

  if( db["ROM_method"] )
  {
    tss.current_time_ = start_time;
    tss.current_step_ = 0;
    tss.set_time_disc_parameters();

    // assemble matrices and right hand side at start time
    tnse_rom.assemble_initial_time();

    
    int linear_iteration=0;
    while(!tss.reached_final_time_step())
    {
      tss.current_step_++;

      // set the time parameters
      tss.set_time_disc_parameters();
      // tau may change depending on the time discretization (adaptive time)
      double tau = tss.get_step_length();
      tss.current_time_ += tss.get_step_length();
      // this is used at several places, e.g., in the example file etc.
      TDatabase::TimeDB->CURRENTTIME += tau;

      bool is_velocity = true;
//       LoopInfo<Residuals> loop_info("nonlinear");
//       loop_info.print_time_every_step = true;
//       loop_info.verbosity_threshold = 1;
      for(unsigned int i=0;; i++)
      {
        tnse_rom.assemble_matrices_rhs(i);

        if(tnse_rom.stop_it(i))
        {
//           loop_info.finish(i, tnse_rom.get_residuals());
          linear_iteration +=i;
//           loop_info_time.print(linear_iteration, tnse_rom.get_full_residual());
          break;
        }
//         else
//         {
//           loop_info.print(i, tnse_rom.get_residuals());
//         }

        tnse_rom.solve(is_velocity);
//         tnse_rom.assemble_matrices_rhs(i+1);
      }
      tnse_rom.solve(!is_velocity);
    }
  } // if ROM

  //----------------------------------------------------------------------------
  // test snapshots average and eigenvalues from POD
  // test reduced and full solutions from ROM
  // delete temporary files used for this test
  //----------------------------------------------------------------------------
//   test_mean_tnse();//read value in file tnse_pod.get_pod_c()->get_snaps_avr());
//   test_eigs_tnse();//read value in filetnse_pod.get_pod_c()->get_pod_eigs());
//   test_sols_tnse(tnse_rom);

  if(remove_tmp_test_files(db, {".u", ".p"}) != 0)
  {
    ErrThrow("snaps_pod_rom_test: TNSE Error deleting temporary files.");
  }
}

/** ***************************************************************************/
int main(int, char**)
{
  parmoon::parmoon_initialize();

  //============================================================================
  // test with TCD2D
  //============================================================================
  test_pod_rom_tcd();

  //============================================================================
  // test with TNSE2D
  //============================================================================
  test_pod_rom_tnse();

  parmoon::parmoon_finalize();
}
