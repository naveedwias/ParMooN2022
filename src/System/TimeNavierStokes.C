#include "TimeNavierStokes.h"
#include "Database.h"

#include "LocalProjection.h"
#include "Hotfixglobal_AssembleNSE.h"
#include "GridTransfer.h"
#include "Multigrid.h"
#include "Variational_multiscale.h"
#include "MainUtilities.h"
#include "PointwiseAssemblyData.h"
#include "NonNewtonianViscosity.h"

#ifdef __2D__
#include "Upwind.h"

#include "Matrix2D.h"
#include "SquareMatrix2D.h"
#include "Assemble2D.h"
#include "Assemble_DG.h"
#include "AuxParam2D.h"
#include "BoundaryAssembling2D.h"
#else
#include "BoundaryAssembling3D.h"
#include "Upwind3D.h"
#include "Matrix3D.h"
#include "SquareMatrix3D.h"
#include "Assemble3D.h"
#include "AuxParam3D.h"
#include "BoundaryAssembling3D.h"
#endif
#ifdef _MPI
#include "ParFECommunicator3D.h"
#include "MumpsWrapper.h"
#endif

/* ************************************************************************** */
template <int d>
ParameterDatabase TimeNavierStokes<d>::default_tnse_database(bool complete)
{
  Output::print<5>("creating a default TNSE parameter database");
  // we use a parmoon default database because this way these parameters are
  // available in the default tnse database as well.
  ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.set_name("TimeNavierStokes parameter database");
  if(complete)
  {
    db.merge(TDomain::default_domain_parameters());
    if(d == 3)
      db.add_nested_database(TDomain::default_sandwich_grid_parameters());
  }
  db.merge(Example_TimeNSE::default_example_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  db.merge(LocalAssembling<d>::default_local_assembling_database());
  db.merge(CheckpointIO::default_checkpoint_io_database());
  if(complete)
  {
    db.merge(Solver<>::default_solver_database());
  }
  
  db.add("time_discretization_nonlinear_term", "fully_implicit", 
         "Type of time discretization for the nonlinear term. The default "
         "(fully_implicit) requires a nonlinear loop in each time step. The "
         "implicit-explicit (imex) scheme only requires one linear solve for "
         "each time step with a changing matrix. The fully-explicit approach "
         "puts the nonlinear term completely on the right-hand side and "
         "therefore requires one linear solve per time step with a constant "
         "matrix.",
         {"fully_implicit", "imex", "fully_explicit"});
  
  db.add("extrapolation_type", "linear",
         "This parameter is used whenever you use imex or the fully-explicit "
         "approach to discretize the nonlinear term in time. These schemes "
         "require a velocity solution which can be either the one from the "
         "previous time step (constant) or a linear combination from the "
         "previous two time steps (linear).",
         {"constant", "linear"});
  
  db["problem_type"] = 6; // means time-dependent Navier-Stokes
  return db;
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::System_per_grid::System_per_grid(
  const Example_TimeNSE& example, const TCollection& coll,
  std::pair<int, int> order)
 : velocity_space(new FESpace(&coll, "u", example.get_bc(0), order.first)),
   pressure_space(new FESpace(&coll, "p", example.get_bc(d), order.second))
{
  switch(TDatabase::ParamDB->NSTYPE)
  {
#ifdef __2D__
    case 1:
      matrix = BlockFEMatrix::NSE2D_Type1(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type1(velocity_space, pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE2D_Type2(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type2(velocity_space, pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE2D_Type3(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type3(velocity_space, pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE2D_Type4(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE2D_Type14(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE2D_Type14(velocity_space, pressure_space);
      break;
#else
    case 1:
      matrix = BlockFEMatrix::NSE3D_Type1(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type1(velocity_space, pressure_space);
      break;
    case 2:
      matrix = BlockFEMatrix::NSE3D_Type2(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type2(velocity_space, pressure_space);
      break;
    case 3:
      matrix = BlockFEMatrix::NSE3D_Type3(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type3(velocity_space, pressure_space);
      break;
    case 4:
      matrix = BlockFEMatrix::NSE3D_Type4(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type4(velocity_space, pressure_space);
      break;
    case 14:
      matrix = BlockFEMatrix::NSE3D_Type14(velocity_space, pressure_space);
      mass_matrix = BlockFEMatrix::Mass_NSE3D_Type14(velocity_space, pressure_space);
      break;
#endif
    default:
      ErrThrow("NSTYPE: ", TDatabase::ParamDB->NSTYPE, " is not known");
  }
  rhs = BlockVector(matrix, true);
  solution = BlockVector(matrix, false);
  u = FEVectFunct(velocity_space, "u", solution.block(0), d);
  p = FEFunction(pressure_space, "p", this->solution.block(d));
  solution_m1 = BlockVector(matrix, false);
  u_m1 = FEVectFunct(velocity_space, "u", solution_m1.block(0), d);
  p_m1 = FEFunction(pressure_space, "p", this->solution_m1.block(d));
  solution_m2 = BlockVector(matrix, false);
  u_m2 = FEVectFunct(velocity_space,"u", solution_m2.block(0), d);
  p_m2 = FEFunction(pressure_space, "p", this->solution_m2.block(d));

  tau_m1 = TDatabase::TimeDB->TIMESTEPLENGTH;
  tau_m2 = TDatabase::TimeDB->TIMESTEPLENGTH;

  time_avg_sol = BlockVector(matrix, false);
  u_time_avg = FEVectFunct(velocity_space, "u_t_avg", time_avg_sol.block(0), d);
  p_time_avg = FEFunction(pressure_space, "p_t_avg", 
                          this->time_avg_sol.block(d));

  combined_old_sols = BlockVector(matrix, false);
  comb_old_u = FEVectFunct(velocity_space, "u", combined_old_sols.block(0), d);
  extrapolate_sol = BlockVector(matrix, false);
  extrapolate_u = FEVectFunct(velocity_space, "u", extrapolate_sol.block(0), d);
  extrapolate_p = FEFunction(pressure_space, "p", extrapolate_sol.block(d));

  persistent_data = std::make_shared<PointwiseAssemblyData>();

#ifdef _MPI
  //print some information
  velocity_space->get_communicator().print_info();
  pressure_space->get_communicator().print_info();
#endif
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::System_per_grid::System_per_grid(
  const System_per_grid& other)
 : velocity_space(other.velocity_space), pressure_space(other.pressure_space),
   matrix(other.matrix), mass_matrix(other.mass_matrix), rhs(other.rhs), solution(other.solution),
   solution_m1(other.solution_m1), solution_m2(other.solution_m2),
   time_avg_sol(other.time_avg_sol), combined_old_sols(other.combined_old_sols),
   extrapolate_sol(other.extrapolate_sol)
{
  // the fe functions must be newly created, because copying would mean
  // referencing the BlockVectors in 'other'.
  u = FEVectFunct(velocity_space, "u", solution.block(0), d);
  p = FEFunction(pressure_space, "p", solution.block(d));
  u_m1 = FEVectFunct(velocity_space, "u", solution_m1.block(0), d);
  p_m1 = FEFunction(pressure_space, "p", this->solution_m1.block(d));
  u_m2 = FEVectFunct(velocity_space,"u", solution_m2.block(0), d);
  p_m2 = FEFunction(pressure_space, "p", this->solution_m2.block(d));
  u_time_avg = FEVectFunct(velocity_space, "u_t_avg", time_avg_sol.block(0), d);
  p_time_avg = FEFunction(pressure_space, "p_t_avg",
                          this->time_avg_sol.block(d));
  comb_old_u = FEVectFunct(velocity_space, "u", combined_old_sols.block(0), d);
  extrapolate_u = FEVectFunct(velocity_space, "u", extrapolate_sol.block(0), d);
  extrapolate_p = FEFunction(pressure_space, "p", extrapolate_sol.block(d));

  persistent_data = std::make_shared<PointwiseAssemblyData>(*other.persistent_data);
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::TimeNavierStokes(const TDomain& domain,
                                      const ParameterDatabase& param_db)
: TimeNavierStokes<d>(domain, param_db, Example_TimeNSE(param_db))
{
}

/* ************************************************************************** */
template <int d>
TimeNavierStokes<d>::TimeNavierStokes(const TDomain& domain,
                                      const ParameterDatabase& param_db,
                                      const Example_TimeNSE& ex)
 : db(default_tnse_database()), systems(), outputWriter(param_db), example(ex),
   solver(param_db), defect(), old_residuals(), initial_residual(1e10),
   solve_count(0), time_stepping_scheme(param_db), is_rhs_and_mass_matrix_nonlinear(false),
   checkpoint_io(param_db), time_average_io(checkpoint_io), Lines()
{
  db.merge(param_db);
  this->check_and_set_parameters();

  std::pair<int,int> velo_pres_order(TDatabase::ParamDB->VELOCITY_SPACE,
                                     TDatabase::ParamDB->PRESSURE_SPACE);
  if(db["space_discretization_type"].is("jump_stab"))
  {
    TDatabase::ParamDB->NSTYPE = 4;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 6;
  }
  // set the velocity and pressure spaces
  // this function returns a pair which consists of
  // velocity and pressure order
  this->get_velocity_pressure_orders(velo_pres_order);

  bool usingMultigrid = this->solver.is_using_multigrid();

  auto collections = domain.get_grid_collections();
  TCollection *coll = collections.front(); // finest grid collection

  // create finite element space, functions, matrices, rhs and solution
  // at the finest grid
  this->systems.emplace_back(example, *coll, velo_pres_order);

#ifdef __3D__

  // prepare the spaces and matrices for vms method
  if(db["space_discretization_type"].is("vms_projection"))
  {
    this->prepare_vms(domain);
  }

#endif

  if (usingMultigrid)
  {
    // Construct multigrid object
    auto mg = solver.get_multigrid();
    size_t n_multigrid_levels = mg->get_n_geometric_levels();
    size_t n_grids = collections.size();

    if (n_multigrid_levels > n_grids)
    {
      ErrThrow("Wrong number of grids for multigrid! expecting ",
               n_multigrid_levels, " geometric grids but got", n_grids,".");
    }

    // remove not needed coarser grid from list of collections
    for(size_t i = n_multigrid_levels; i < n_grids; ++i)
    {
      collections.pop_back();
    }

    if (mg->is_using_mdml())
    {
      // change the discretization on the coarse grids to lowest order
      // non-conforming(-1). The pressure space is chosen automatically(-4711).

      velo_pres_order = { -1, -4711 };
      this->get_velocity_pressure_orders(velo_pres_order);
    }
    else
    {
      // for standard multigrid, pop the finest collection - it was already
      // used to construct a space before the "if(usingMultigrid)" clause
      // and will not (as in mdml) be used a second time with a different discretization

      collections.pop_front();
    }

    // Construct systems per grid and store them, finest level first
    std::list<BlockFEMatrix*> matrices;

    // matrix on finest grid is already constructed
    matrices.push_back(&systems.back().matrix);
    for (auto coll: collections)
    {
      systems.emplace_back(example, *coll, velo_pres_order);

      // prepare input argument for multigrid object
      matrices.push_front(&systems.back().matrix);
    }

    // initialize the multigrid object with all the matrices on all levels
    mg->initialize(matrices);
  }

  // initial solution on finest grid - read-in or interpolation
  initialize_solution();

  // the defect has the same structure as the rhs (and as the solution)
  this->defect.copy_structure(this->systems.front().rhs);
  this->rhs_from_time_disc.copy_structure(this->systems.front().rhs);

  outputWriter.add_fe_vector_function(&this->get_velocity());
  outputWriter.add_fe_function(&this->get_pressure());

  if (db["output_compute_time_average"])
  {
    outputWriter.add_fe_vector_function(&systems.front().u_time_avg);
    outputWriter.add_fe_function(&systems.front().p_time_avg);

    std::string base_name = db["write_solution_binary_file"];
    std::string base_name_for_reading = db["initial_solution_file"];
    db["write_solution_binary_file"].set(base_name + "_time_average", false);
    db["initial_solution_file"].set(base_name_for_reading + "_time_average",
                                    false);
    time_average_io = CheckpointIO(db);
    db["write_solution_binary_file"] = base_name; // reset
    db["initial_solution_file"] = base_name_for_reading; // reset
    try
    {
      time_average_io.read(systems.front().time_avg_sol);
    }
    catch(const std::runtime_error& err)
    {
      Output::warn("TNSE", "Reading time averaged solution failed",
                   "Now I start with a zero time average. The error was ",
                   err.what());
    }
  }

  // print out the information (cells, dofs, etc)
  this->output_problem_size_info();

  this->errors.fill(0.);

  if (db["output_along_line"] )
  {
    Lines = LinesEval<d>(domain, param_db);
  }

  stream_function_space
     = std::make_shared<FESpace>(coll, "s", example.get_bc(0), 1);
  n_psi = stream_function_space->get_n_dof();
  psi.resize(n_psi, 0.);
  stream_function
     = std::make_shared<FEFunction>(stream_function_space,"s", psi.data());

  // add to the wrapper
  outputWriter.add_fe_function(stream_function.get());

  zero_vorticity = -4711;
  if (db["compute_vorticity_divergence"])
  {
    vorticity_space
      =std::make_shared<FESpace>(coll, "v", example.get_bc(0), 1);
    n_vort_dofs = vorticity_space->get_n_dof();
    vorticity.resize(2*n_vort_dofs, 0.);
    vorticity_funct 
      = std::make_shared<FEFunction>(vorticity_space, "v", 
                                     vorticity.data()+n_vort_dofs);
    divergence
      = std::make_shared<FEFunction>(vorticity_space, "d", vorticity.data());
    outputWriter.add_fe_function(vorticity_funct.get());
    outputWriter.add_fe_function(divergence.get());
  }

  if (db["compute_wall_shear_stress"])
  {
#ifdef __3D__

    wall_shear_stress_space
      = std::make_shared<FESpace>(coll, "tau_w", example.get_bc(0), 1);

    n_wall_shear_stress_dofs = wall_shear_stress_space->get_n_dof();

    wall_shear_stress.resize(d * n_wall_shear_stress_dofs, 0.0);

    wall_shear_stress_funct = std::make_shared<FEVectFunct>(
      wall_shear_stress_space, "tau_w", wall_shear_stress.data(), d);

    outputWriter.add_fe_vector_function(wall_shear_stress_funct.get());

#else
    Output::root_warn("TNSE", "Wall shear stress computation is currently only implemented in 3D.");
#endif
  }

  if (db["compute_time_derivative"])
  {
    time_derivative = BlockVector(systems.front().matrix);

    u_time_derivative_funct = std::make_shared<FEVectFunct>(
      systems.front().velocity_space, "dt_u", time_derivative.block(0), d);

    p_time_derivative_funct = std::make_shared<FEFunction>(
      systems.front().velocity_space, "dt_p", time_derivative.block(d));

    outputWriter.add_fe_vector_function(u_time_derivative_funct.get());
    outputWriter.add_fe_function(p_time_derivative_funct.get());
  }

  if (db["compute_turbulent_kinetic_energy"])
  {
  #ifdef __3D__

    if (!db["space_discretization_type"].is("smagorinsky"))
    {
      Output::root_warn("TNSE", "Turbulent kinetic energy computation is "
        "currently only implemented for LES discretizations "
        "(\"smagorinsky\").");
    }
    else
    {
      tke_space = systems.front().velocity_space;

      n_tke_dofs = tke_space->get_n_dof();

      tke.resize(n_tke_dofs, 0.0);

      tke_funct = std::make_shared<FEFunction>(tke_space, "k", tke.data());

      outputWriter.add_fe_function(tke_funct.get());
    }
#else
    Output::root_warn("TNSE", "Turbulent kinetic energy computation is "
      "currently only implemented in 3D.");
#endif
  }

  if (db["compute_effective_viscosity"])
  {
  #ifdef __3D__

    if (db["viscosity_mode"].is("newtonian"))
    {
      Output::root_warn("TNSE", "Effective viscosity computation is "
        "only functional for non-newtonian viscosities.");
    }
    else
    {
      nu_eff_space = systems.front().velocity_space;

      n_nu_eff_dofs = nu_eff_space->get_n_dof();

      nu_eff.resize(n_nu_eff_dofs, 0.0);

      nu_eff_funct = std::make_shared<FEFunction>(nu_eff_space, "nu_eff", nu_eff.data());

      outputWriter.add_fe_function(nu_eff_funct.get());
    }
#else
    Output::root_warn("TNSE", "Effective viscosity computation is "
      "currently only implemented in 3D.");
#endif
  }

  if (db.contains("output_boundary_flow_rates"))
  {
    std::vector<size_t> flow_bc = db["output_boundary_flow_rates"];

    for (unsigned int i = 0; i < flow_bc.size(); i++)
    {
      unsigned int bd = flow_bc[i];
      std::ostringstream ostr;
      ostr << "Q_" << bd;

      outputWriter.add_global_output_variable(ostr.str(), [&, bd] () -> double
        {
          return systems.front().u.compute_flux(bd);
        });
    }
  }

  if (db["output_compute_errors"])
  {
    last_kinetic_energy = 0.0;

    outputWriter.add_global_output_variable("E_k", [&] () -> double
        {
          return last_kinetic_energy;
        });

    last_tke = 0.0;

    if (db["space_discretization_type"].is("smagorinsky")
      || db["space_discretization_type"].is("residual_based_vms")
      || db["space_discretization_type"].is("rbvms_time"))
    {
      outputWriter.add_global_output_variable("k_sgs", [&] () -> double
          {
            return last_tke;
          });
    }

    last_max_velocity = 0.0;

    outputWriter.add_global_output_variable("u_max", [&] () -> double
        {
          return last_max_velocity;
        });
  }

  // initialize WK boundaries
  unsigned int n_windkessel_bd = param_db.try_get_value("n_windkessel_bd", 0);

  if (n_windkessel_bd)
  {
    // Windkessel BC
    std::vector<size_t> windkessel_id = param_db["windkessel_id"];

    if (windkessel_id.size() < n_windkessel_bd)
    {
      ErrThrow("Not enough data for ", n_windkessel_bd, " WK boundary conditions");
    }

    std::vector<double> windkessel_Rp = db.contains("windkessel_Rp") ? db["windkessel_Rp"] : std::vector<double>();
    std::vector<double> windkessel_Rd = db.contains("windkessel_Rd") ? db["windkessel_Rd"] : std::vector<double>();
    std::vector<double> windkessel_C = db.contains("windkessel_C") ? db["windkessel_C"] : std::vector<double>();
    std::vector<double> windkessel_pi = db.contains("windkessel_initial_distal_pressure") ? db["windkessel_initial_distal_pressure"] : std::vector<double>();

    if (windkessel_Rp.size() < n_windkessel_bd)
    {
      Output::root_warn("TNSE", "windkessel_Rp has fewer than ", n_windkessel_bd, " elements. Remainder will default to zero.");  
    }

    if (windkessel_Rd.size() < n_windkessel_bd)
    {
      Output::root_warn("TNSE", "windkessel_Rd has fewer than ", n_windkessel_bd, " elements. Remainder will default to zero.");  
    }

    if (windkessel_C.size() < n_windkessel_bd)
    {
      Output::root_warn("TNSE", "windkessel_C has fewer than ", n_windkessel_bd, " elements. Remainder will default to zero.");  
    }

    if (windkessel_pi.size() < n_windkessel_bd)
    {
      Output::root_warn("TNSE", "windkessel_initial_distal_pressure has fewer than ", n_windkessel_bd, " elements. Remainder will default to zero.");  
    }

    for (unsigned int k = 0; k < n_windkessel_bd; k++)
    {
      double Rp_k = k < windkessel_Rp.size() ? windkessel_Rp[k] : 0.0;
      double Rd_k = k < windkessel_Rd.size() ? windkessel_Rd[k] : 0.0;
      double C_k = k < windkessel_C.size() ? windkessel_C[k] : 0.0;
      double pi_k = k < windkessel_pi.size() ? windkessel_pi[k] : 0.0;

      Output::root_info("TNSE", "Initialize WK: ", Rp_k, " " , C_k, " ", Rd_k, " p_distal ", pi_k);

      WindkesselBoundary wk(Rp_k, C_k, Rd_k, pi_k);
      wk_boundary.push_back(wk);
      wk_p_out.push_back(wk.get_pressure());
    }

    if (db.try_get_value("output_windkessel_pressures", false))
    {
      for (unsigned int k = 0; k < n_windkessel_bd; k++)
      {
        unsigned int bd = windkessel_id[k];

        std::ostringstream ostr;
        ostr << "WK_P_" << bd;

        outputWriter.add_global_output_variable(ostr.str(), [&, k] () -> double
        {
          return wk_p_out[k];
        });
      }
    }

    Output::root_info("TNSE", "I have initialized ", wk_boundary.size(), " Windkessel boundaries");
  }

  bool diagnostics = db["nonlinloop_diagnostics"];
  std::ofstream diag_f;

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_rank == 0)
  {
#endif
    if (diagnostics)
    {
      diag_f = std::ofstream(outputWriter.get_base_name() + "_nonlinloop.csv", std::ios::trunc);
    }
#ifdef _MPI
  }
#endif

  if (diagnostics)
  {
    diag_f << "Time,Tau,Damping,Residual\n";
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::check_and_set_parameters()
{
  if (!db["problem_type"].is(6))
  {
    if (db["problem_type"].is(0))
    {
      db["problem_type"] = 6;
    }
    else
    {
      Output::root_warn<2>("TNSE", "The parameter problem_type doesn't correspond to "
                      "TimeNavierStokes. It is now reset to the correct value "
                      "for TimeNavierStokes (=6).");
      db["problem_type"] = 6;
    }
  }

  if (db["time_discretization"].is("forward_euler"))
  {
    ErrThrow("time discretization: ", db["time_discretization"],
             " is not supported");
  }

  if (db["time_discretization_nonlinear_term"].is("imex") 
     || db["time_discretization_nonlinear_term"].is("fully_explicit"))
  {
    // turn off nonlinear damping, because in this case there is no nonlinear
    // iteration, only one step
    if (!db["nonlinloop_damping_factor"].is(1.0))
    {
      Output::root_info("TNSE", "resetting nonlinloop_damping_factor to 1.0, "
                   "because time_discretization_nonlinear_term is set to ",
                   db["time_discretization_nonlinear_term"], ".");
    }

    db["nonlinloop_damping_factor"] = 1.0;

    if (!db["nse_nonlinear_form"].is("convective"))
    {
      ErrThrow("using ", db["time_discretization_nonlinear_term"], 
               " for the 'time_discretization_nonlinear_term' currently "
               "requires the 'nse_nonlinear_form' to be 'convective' rather "
               "than ", db["nse_nonlinear_form"], ".");
    }
  }

  // set the discretization parameters
  // standard Galerkin
  if (db["space_discretization_type"].is("galerkin") || 
      db["space_discretization_type"].is("jump_stab"))
  {
    space_disc_global = 1;
    /// set scaling factor for B, BT's block
    time_stepping_scheme.n_scale_block = 2*d;
    time_stepping_scheme.b_bt_linear_nl = "linear";
  }

  if (db["space_discretization_type"].is("supg"))
  {
    // supg: NOTE: only tested with BDF2 so far
    if(!db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }
    if(TDatabase::ParamDB->NSTYPE < 4)
    {
      ErrThrow("supg method is only supported for type 4, and 14 ");
    }
    space_disc_global = 2;
    /// set scaling factor for B, BT's block
    // depends on how to deal the nonlinearity in the 
    // test function: fully implicit case
    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    time_stepping_scheme.n_scale_block = d;
    if(TDatabase::ParamDB->NSTYPE == 14)
      time_stepping_scheme.n_scale_block = 2*d+1;
  }

  // Smagorinsky
  if (db["space_discretization_type"].is("smagorinsky"))
  {
    //space_disc_global = 4;
    time_stepping_scheme.n_scale_block = 6;
    time_stepping_scheme.b_bt_linear_nl = "linear";
  }

  if (db["space_discretization_type"].is("local_projection"))
  {
    if(d == 3)
      ErrThrow("local_projection in 3D not implemented");
    space_disc_global = 14;
  }

  // the only case where one have to re-assemble the right hand side
  if (db["space_discretization_type"].is("supg")
     && db["time_discretization"].is("bdf_two"))
  {
    is_rhs_and_mass_matrix_nonlinear = true;
  }

  if (db["space_discretization_type"].is("residual_based_vms")
    || db["space_discretization_type"].is("rbvms_time"))
  {
    is_rhs_and_mass_matrix_nonlinear = true;

    time_stepping_scheme.b_bt_linear_nl = "nonlinear";
    time_stepping_scheme.n_scale_block = d;
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::get_velocity_pressure_orders(
  std::pair<int, int> &velo_pres_order) 
{
  int velocity_order = velo_pres_order.first;
  int pressure_order = velo_pres_order.second;
  int order = 0;
  switch(velocity_order)
  {
    case 1: case 2: case 3: case 4: case 5:
    case 12: case 13: case 14: case 15:
      if(velocity_order > 10)
        order = velocity_order-10;
      else
        order = velocity_order;
      break;
    case -1: case -2: case -3: case -4: case -5: case -101:
      order = velocity_order;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("velocity_order ", velocity_order, " not supported in 3D");
      order = velocity_order;
      break;
    // conforming fe spaces with bubbles on triangles
    case 22: case 23: case 24:
      order = velocity_order;
      break;
      // discontinuous spaces
    case -11: case -12: case -13:
      order = velocity_order*10;
      break;
  }
  TDatabase::ParamDB->VELOCITY_SPACE = order;
  velo_pres_order.first = order;
  switch(pressure_order)
  {
    case -4711:
    {
      switch(velocity_order)
      {
        case -1: case -2: case -3: case -4:
          // nonconforming pw (bi)linear velo/ pw constant pressure
          // conforming pw (bi)linear velo/ pw constant pressure (not stable !!!)
          pressure_order = -velocity_order-1;
          break;
        case 1: // discontinuous space
          pressure_order = 0;
          break;
        case 2: case 3: case 4: case 5:
        // standard conforming velo and continuous pressure
          pressure_order = velocity_order-1;
          break;
          // discontinuous pressure spaces
          // standard conforming velo and discontinuous pressure
          // this is not stable on triangles !!!
        case 12: case 13: case 14: case 15:
        case -11: case -12: case -13: case -14:
          pressure_order = -(velocity_order-1)*10;
          break;
        case 22: case 23: case 24:
          pressure_order = -(velocity_order-11)*10;
          break;
        case 100: case 201: case 302: case 403: case 504:
          pressure_order = -(velocity_order%100 + 10)*10;
          break; 
      }
      break;
    }
    case 1: case 2: case 3: case 4: case 5:
      // pressure order is chosen correctly
      break;
    // discontinuous spaces
    case -11: case -12: case -13: case -14:
      pressure_order = pressure_order*10;
      break;
    case 100: case 201: case 302: case 403: case 504:
      if(d == 3)
        ErrThrow("pressure_order ", pressure_order, " not supported in 3D");
      // pressure order is chosen correctly
      break;
    default:
      ErrThrow("pressure space is not chosen properly ", pressure_order);
  }
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  velo_pres_order.second = pressure_order;

  Output::root_info("TNSE", "velocity space order: ", setw(6), velo_pres_order.first);
  Output::root_info("TNSE", "pressure space order: ", setw(6), velo_pres_order.second);
}

/* ************************************************************************** */
template <int d>
bool TimeNavierStokes<d>::initialize_solution()
{
  bool solution_read = checkpoint_io.read(systems.front().solution);

  if (solution_read)
  {
    if (this->time_stepping_scheme.get_start_time() == 0.)
    {
      Output::root_warn<1>("Initial Solution",
        "Restarting from existing solution but initial time is 0! This is "
        "probably not what you want! If for example your BC or RHS are "
        "time-dependent, you will apply initial BC or RHS to an already "
        "developed flow, instead of continuing your simulation! Set your "
        "time_start to the time of the binary file you are re-starting from "
        "(and don't forget to set continue_output_after_restart to true, and "
        "also read_metis parameters if needed). Or ignore this warning if you "
        "know what you are doing and are aware of the consequences (rhs and bc "
        "reset to t = 0, output restart from 0 and probably overwrites outputs "
        "from old simulation).");
    }
  }
  else
  {
    Output::root_info("Initial Solution", "Interpolating initial solution from example.");
    for (System_per_grid& s : this->systems)
    {
      for (int i = 0; i < d; ++i)
      {
        auto ui = s.u.GetComponent(i);
        ui->Interpolate(example.get_initial_cond(i));
      }

      s.p.Interpolate(example.get_initial_cond(d));

      if (s.matrix.pressure_projection_enabled())
      {
        s.p.project_into_L20();
      }
    }
  }

  return solution_read;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_initial_time()
{
  solve_count = 0;

  if (systems.size() > 1) // using multigrid
  {
    this->restrict_function();
  }

  for (auto &s: this->systems)
  {
    call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesMass);
  }

  bool reset_to_vms = false;

  for (auto &s : this->systems)
  {
    s.solution_m1 = s.solution;
    s.solution_m2 = s.solution;

    s.tau_m1 = TDatabase::TimeDB->TIMESTEPLENGTH;
    s.tau_m2 = TDatabase::TimeDB->TIMESTEPLENGTH;

    if (!db["time_discretization_nonlinear_term"].is("fully_explicit"))
    {
      call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesAll);
    }
    else
    {
      call_assembling_routine(s, LocalAssembling_type::NavierStokesAll);
    }

    // update matrices for local projection stabilization
    if (db["space_discretization_type"].is("local_projection"))
    {
      update_matrices_lps(s);
    }

    adjust_assembled_rhs();

    // copy nonactives
    s.solution.copy_nonactive(s.rhs);

    if (db["space_discretization_type"].is("vms_projection"))
    {
      std::vector<std::shared_ptr<FEMatrix>> blocks
                      = s.matrix.get_blocks_uniquely();

      // update mass matrix of projection
      LumpMassMatrixToDiagonalMatrix<d>(matrices_for_turb_mod.at(6));

      // update stiffness matrix
      VMS_ProjectionUpdateMatrices<d>(blocks, matrices_for_turb_mod);

      // reset flag for projection-based VMS method such that Smagorinsky LES method
      // is used on coarser grids
      space_disc_global = 4;
      db["space_discretization_type"] = "smagorinsky";
      reset_to_vms = true;
    }

    /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
     * MPI-case in order to be consistently saved. It is necessary that the vector
     * is consistently saved because it is the only way to ensure that its
     * multiplication with an inconsistently saved matrix (multiplication which
     * appears in the defect and rhs computations) give the correct results. When
     * we call copy_nonactive in MPI-case, we have to remember the following: it
     * can happen that some slave ACTTIVE DoFs are placed in the block of
     * NON-ACTIVE DoFs (because they are at the interface between processors).
     * Doing copy_nonactive changes then the value of these DOFs,although they are
     * actually active. That's why we have to update the values so that the vector
     * becomes consistent again.
     */
#ifdef _MPI
    for (int i = 0; i < d; ++i)
    {
      double *ui = s.solution.block(i);
      s.velocity_space->get_communicator().consistency_update(ui, 3);
    }
    double *p = s.solution.block(d);
    s.pressure_space->get_communicator().consistency_update(p, 3);
#endif
  }

  // reset DISCTYPE to VMS_PROJECTION to be correct in the next assembling
  if (reset_to_vms)
  {
    space_disc_global = 9;
    db["space_discretization_type"] = "vms_projection";
  }

  // the call to assembleslip is necessary here, in order to get
  // the correct old_rhs, i.e., zeros on the slip dofs
  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    this->modify_slip_bc(true, true);
  }

  auto& s = this->systems.front();

  // copy the current right hand side vector to the old_rhs
  base_rhs = s.rhs;
  base_rhs_m1 = s.rhs;
  base_rhs_m2 = s.rhs;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_matrices_rhs(unsigned int it_counter)
{
  bool first_nonlinear_step = (it_counter == 0);

  if (first_nonlinear_step)
  {
    solve_count = 0;
  }

  if (!first_nonlinear_step && linearized())
  {
    // no further assembling needed here
    return;
  }

  bool fully_explicit =
    this->db["time_discretization_nonlinear_term"].is("fully_explicit");

  // note: for fully_explicit==true, the matrix is constant time, still we 
  // call reset_linear_matrices and later prepare_system_matrix which will
  // produce the same matrix. The reason is that we need the individual parts
  // (mass and laplacian) to build the right-hand side.
  reset_residuals();

  if (first_nonlinear_step)
  {
    // Important: We have to descale the matrices, since they are scaled
    // before the solving process. Only A11 and A22 matrices are
    // reset and assembled again but the A12 and A21 are scaled, so
    // for the next iteration we have to descale
    if (time_stepping_scheme.current_step_ > 1)
    {
      for (System_per_grid & s: this->systems)
      {
        time_stepping_scheme.reset_linear_matrices(s.matrix, s.mass_matrix);
      }
    }

    // initialize the rhs from the time discretization
    System_per_grid& s = this->systems.front();

    // only assemble the right-hand side
    call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesRhs,
                            first_nonlinear_step);

    // copy the non active to the solution vector
    // since the rhs vector will be passed to the solver
    // and is modified with matrix vector multiplication
    // which also uses the non-actives
    adjust_assembled_rhs();
    s.solution.copy_nonactive(s.rhs);

    // all matrices from the previous time step are available
    unsigned int n_sols = time_stepping_scheme.n_old_solutions();
    std::vector<BlockVector> oldsolutions(n_sols);

    oldsolutions[0] = s.solution_m1;
    if (oldsolutions.size() == 2)
    {
      oldsolutions[1] = s.solution_m2;
    }

    // modification of the matrices due to slip b.c
    if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
    {
      this->modify_slip_bc(true, true);
    }

    // one needs two right hand sides only for the crank-Nicolson
    // and fractional step theta schemeis
    std::vector<BlockVector> all_rhs(2);
    all_rhs[0] = s.rhs; // current rhsw
    all_rhs[1] = base_rhs_m1;

    if (s.tau_m1 != TDatabase::TimeDB->TIMESTEPLENGTH)
    {
      // KLUDGE: inter-/extrapolate.
      // rhs_{-1}^* = rhs + tau / tau_{-1} * (rhs_{-1} - rhs)

      // TODO: replace with reassembly of previous time's rhs.

      double tau_ratio = TDatabase::TimeDB->TIMESTEPLENGTH / s.tau_m1;

      all_rhs[1].scale(tau_ratio);
      all_rhs[1].add_scaled(all_rhs[0], 1.0 - tau_ratio);
    }

    base_rhs = s.rhs;

    // NOTE: scale the B blocks only at the first iteration
    for (System_per_grid& sys: this->systems)
    {
      time_stepping_scheme.scale_descale_all_b_blocks(sys.matrix);
    }

    time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                                                    all_rhs, oldsolutions);

    if (fully_explicit)
    {
      // assemble the nonlinear term with given functions for the velocity
      s.rhs.reset();
      call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesExplNL,
                              first_nonlinear_step);
      all_rhs[0].addScaledActive(s.rhs,
                                 -time_stepping_scheme.get_step_length());
    }

    all_rhs[0].copy_nonactive(s.solution);
    rhs_from_time_disc = all_rhs[0];

    /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
     * MPI-case in order to be consistently saved. It is necessary that the vector
     * is consistently saved because it is the only way to ensure that its
     * multiplication with an inconsistently saved matrix (multiplication which
     * appears in the defect and rhs computations) give the correct results. When
     * we call copy_nonactive in MPI-case, we have to remember the following: it
     * can happen that some slave ACTTIVE DoFs are placed in the block of
     * NON-ACTIVE DoFs (because they are at the interface between processors).
     * Doing copy_nonactive changes then the value of these DOFs,although they are
     * actually active. That's why we have to update the values so that the vector
     * becomes consistent again.
     */
#ifdef _MPI
    for (int i = 0; i < d; ++i)
    {
      double *ui = s.solution.block(i);
      s.velocity_space->get_communicator().consistency_update(ui, 3);
      s.velocity_space->get_communicator().consistency_update(rhs_from_time_disc.block(i), 3);
    }
    double *p = s.solution.block(d);
    s.pressure_space->get_communicator().consistency_update(p, 3);
    s.pressure_space->get_communicator().consistency_update(rhs_from_time_disc.block(d), 3);
#endif // _MPI
  }

  // assemble the nonlinear matrices
  if (systems.size() > 1)
  {
    this->restrict_function();
  }

  bool do_upwinding = false;
  if (!fully_explicit)
  {
    bool reset_to_vms = false;

    for (System_per_grid & s: systems)
    {
      bool mdml = this->solver.is_using_multigrid()
                && this->solver.get_multigrid()->is_using_mdml();

      bool on_finest_grid = &systems.front() == &s;
      do_upwinding = (db["space_discretization_type"].is("upwind")
                      || (mdml && !on_finest_grid));

      if (do_upwinding)
      {
        /**
        * This will assemble the linear blocks of velocity-velocity 
        * coupling used only in the mdml case. One can extend this 
        * also for the assembling of B-blocks. This is done separately 
        * here because in the LocalAssembling_type::TimeNavierStokesAll
        * we assemble the B-blocks as well which were scaled with the 
        * corresponding factor from the time stepping schemes and therefore
        * needs not to be re-assemble here once more. 
        */
        call_assembling_routine(s, LocalAssembling_type::NavierStokesLinear);
      }
      else
      {
        if (TDatabase::ParamDB->NSTYPE == 14
          && (db["space_discretization_type"].is("supg")
           || db["space_discretization_type"].is("residual_based_vms")
           || db["space_discretization_type"].is("rbvms_time")))
        {
          call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesAll,
            first_nonlinear_step);

          adjust_assembled_rhs();
          s.solution.copy_nonactive(s.rhs);
        }
        else if (is_rhs_and_mass_matrix_nonlinear && !first_nonlinear_step)
        {
          call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesAll);

          adjust_assembled_rhs();
          s.solution.copy_nonactive(s.rhs);
        }
        else
        {
          call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesNL);
        }

        // assemble the mass matrix in each nonlinear iteration 
        if (db["space_discretization_type"].is("supg")
          || db["space_discretization_type"].is("residual_based_vms")
          || db["space_discretization_type"].is("rbvms_time"))
        {
          call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesMass);
        }
      }

      // update matrices with local projection term
      if (db["space_discretization_type"].is("local_projection"))
      {
        update_matrices_lps(s);
      }

      if (db["space_discretization_type"].is("vms_projection"))
      {
        std::vector<std::shared_ptr<FEMatrix>> blocks
           = s.matrix.get_blocks_uniquely();

        // update mass matrix of projection
        LumpMassMatrixToDiagonalMatrix<d>(matrices_for_turb_mod.at(6));

        // update stiffness matrix
        VMS_ProjectionUpdateMatrices<d>(blocks, matrices_for_turb_mod);

        // reset flag for projection-based VMS method such that Smagorinsky LES method
        // is used on coarser grids 
        space_disc_global = 4;
        db["space_discretization_type"].set("smagorinsky");

        reset_to_vms = true;
      }
    }

    // reset   DISCTYPE to VMS_PROJECTION to be correct in the next assembling
    if (reset_to_vms)
    {
      space_disc_global = 9;
      db["space_discretization_type"].set("vms_projection");
    }
  }

  // slip boundary modification of matrices
  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
  {
    if (is_rhs_and_mass_matrix_nonlinear)
    {
      this->modify_slip_bc(true, true);
    }
    else
    {
      this->modify_slip_bc();
    }
  }

  // prepare the right-hand side if it is nonlinear, assembling already 
  // done together with the matrices
  if (is_rhs_and_mass_matrix_nonlinear)
  {
    this->assemble_rhs_nonlinear();
  }

  // also prepare the system matrices for the solver
  for (System_per_grid& s: this->systems)
  {
    // call the preparing method
    time_stepping_scheme.prepare_system_matrix(s.matrix, s.mass_matrix);

    if (db["space_discretization_type"].is("supg")
      || db["space_discretization_type"].is("residual_based_vms")
      || db["space_discretization_type"].is("rbvms_time"))
    {
      time_stepping_scheme.scale_nl_b_blocks(s.matrix);
    }
  }

  Output::print<5>("Assembling of matrices and right hand side is done");
}

template <int d>
void TimeNavierStokes<d>::adjust_assembled_rhs()
{
  Output::root_info("TNSE", "No RHS adjustment necessary.");
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::assemble_rhs_nonlinear()
{
  System_per_grid& s = this->systems.front();

  base_rhs = s.rhs;

  // all matrices from the previous time step are available
  unsigned int n_sols = time_stepping_scheme.n_old_solutions();

  std::vector<BlockVector> oldsolutions(n_sols);

  oldsolutions[0] = s.solution_m1;
  if (oldsolutions.size() == 2)
  {
    oldsolutions[1] = s.solution_m2;
  }

  // one needs two right hand sides only for the crank-Nicolson
  // and fractional step theta schemes
  std::vector<BlockVector> rhs_(2);
  rhs_[0] = base_rhs; // current rhs
  rhs_[1] = base_rhs_m1; // old right hand side is needed for the Crank-Nicolson time stepping

  if (s.tau_m1 != TDatabase::TimeDB->TIMESTEPLENGTH)
  {
    // KLUDGE: inter-/extrapolate.
    // rhs_{-1}^* = rhs + tau / tau_{-1} * (rhs_{-1} - rhs)

    // TODO: replace with reassembly of previous time's rhs.

    double tau_ratio = TDatabase::TimeDB->TIMESTEPLENGTH / s.tau_m1;

    rhs_[1].scale(tau_ratio);
    rhs_[1].add_scaled(rhs_[0], 1.0 - tau_ratio);
  }

  // prepare the right hand side for the solver
  time_stepping_scheme.prepare_rhs_from_time_disc(s.matrix, s.mass_matrix,
                                                  rhs_, oldsolutions);
  rhs_from_time_disc = rhs_[0];

  // TODO: this should completely ruin the iteration because base_rhs_m1
  // (previously named old_rhs) should be the previous time step's base rhs
  // but gets overwritten on every step of the nonlinear iteration!
  // figure out what's actually being done here and what *should* be
  // base_rhs_m1 = s.rhs;

  // copy the non-actives
  rhs_from_time_disc.copy_nonactive(s.rhs);
  s.solution.copy_nonactive(s.rhs);

  /** After copy_nonactive, the solution vectors needs to be Comm-updated in 
   * MPI-case in order to be consistently saved. It is necessary that the vector
   * is consistently saved because it is the only way to ensure that its
   * multiplication with an inconsistently saved matrix (multiplication which
   * appears in the defect and rhs computations) give the correct results. When
   * we call copy_nonactive in MPI-case, we have to remember the following: it
   * can happen that some slave ACTTIVE DoFs are placed in the block of
   * NON-ACTIVE DoFs (because they are at the interface between processors).
   * Doing copy_nonactive changes then the value of these DOFs,although they are
   * actually active. That's why we have to update the values so that the vector
   * becomes consistent again.
   */
#ifdef _MPI
  for (int i = 0; i < d; ++i)
  {
    double *ui = s.solution.block(i);
    s.velocity_space->get_communicator().consistency_update(ui, 3);

    s.velocity_space->get_communicator().consistency_update(rhs_from_time_disc.block(i), 3);
  }

  double *p = s.solution.block(d);
  s.pressure_space->get_communicator().consistency_update(p, 3);
  s.pressure_space->get_communicator().consistency_update(rhs_from_time_disc.block(d), 3);
#endif // _MPI
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::call_assembling_routine(
  TimeNavierStokes<d>::System_per_grid& s, LocalAssembling_type type,
  bool first_nonlinear_step, bool do_upwinding, bool assemble_dirichlet_rows)
{
  Output::root_info("TNSE", "Assembling: ", type, ", ", first_nonlinear_step);

  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;

  // set arrays of spaces for matrices and rhs
  std::vector<const FESpace*> spaces_mat;
  std::vector<const FESpace*> spaces_rhs;
  std::vector<const FEFunction*> fefunctions;

  // call to routine to set arrays
  set_arrays(s, spaces_mat, spaces_rhs, fefunctions, first_nonlinear_step);

  // prepare matrices and rhs for assembling
  std::vector<SquareMatrixD*> sqMatrices;
  std::vector<MatrixD*> rectMatrices;
  std::vector<double*> rhs_array;

  // call the routine to prepare the matrices
  set_matrices_rhs(s, type, sqMatrices, rectMatrices, rhs_array);

  std::array<BoundaryConditionFunction*, d + 1> boundCondition;
  std::array<BoundaryValuesFunction*, d + 1> boundValues;

  for (int i = 0; i < d; ++i)
  {
    boundCondition[i] = s.velocity_space->get_boundary_condition();
  }

  boundCondition[d] = s.pressure_space->get_boundary_condition();

  for (int i = 0; i < d + 1 ; ++i)
  {
    boundValues[i] = example.get_bd(i);
  }

  // local assembling settings
  LocalAssembling<d> la(this->db, type, fefunctions,
                       this->example.get_coeffs(), space_disc_global);

  la.SetPersistentData(s.persistent_data);

#ifdef __3D__
    Assemble3D(
#else
    Assemble2D(
#endif
               spaces_mat.size(), spaces_mat.data(), sqMatrices.size(),
               sqMatrices.data(), rectMatrices.size(), rectMatrices.data(),
               rhs_array.size(), rhs_array.data(), spaces_rhs.data(),
               boundCondition.data(), boundValues.data(), la,
               assemble_dirichlet_rows);

    if(db["space_discretization_type"].is("jump_stab") && 
      type == LocalAssembling_type::TimeNavierStokesNL)
    {
      AddJumpStabilizationCIP(1, spaces_mat.data(), 1, sqMatrices.data(), 
            boundCondition.data(), boundValues.data(), db, la);
    }
  /**
     @todo currently the function call_assembling_routine is
     called several times throughout the code, with different
     LocalAssembling types. We need to know the type here
     to decide whether the boundary terms should be assembled.
     The following part should be checked in general cases. 
     Another option would be to pass the type to the bd function, i.e.
     assemble_boundary_terms(type)
  **/
  if ((type == LocalAssembling_type::TimeNavierStokesAll
    || type == LocalAssembling_type::NavierStokesAll
    || type == LocalAssembling_type::TimeNavierStokesRhs
    || type == LocalAssembling_type::TimeNavierStokesExplNL)
    && !assemble_dirichlet_rows)
  {
    // boundary terms (weak BC, Neumann BC)
    assemble_boundary_terms(first_nonlinear_step);
  }

  if (do_upwinding)
  {
    switch (TDatabase::ParamDB->NSTYPE)
    {
      case 1:
      case 2:
      {
        // do upwinding with one matrix
#ifdef __2D__
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
#else
        // the inverse of the example's diffusion coefficient
        double one_over_nu = 1./example.get_nu();
        UpwindForNavierStokes3D(sqMatrices[0], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
#endif
        Output::print<3>("UPWINDING DONE");
        break;
      }

      case 3:
      case 4:
      case 14:
      {
        // do upwinding with d matrices
#ifdef __2D__
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[0],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
        UpwindForNavierStokes(la.GetCoeffFct(), sqMatrices[1],
                              la.get_fe_function(0),
                              la.get_fe_function(1));
#else
        // the inverse of the example's diffusion coefficient
        double one_over_nu = 1./example.get_nu();
        UpwindForNavierStokes3D(sqMatrices[0], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
        UpwindForNavierStokes3D(sqMatrices[1], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
        UpwindForNavierStokes3D(sqMatrices[2], la.get_fe_function(0),
                                la.get_fe_function(1), la.get_fe_function(2),
                                one_over_nu);
#endif
        Output::print<3>("UPWINDING DONE");
        break;
      }
    } // endswitch
  }

  for (int i = 0; i < d; i++)
  {
    delete fefunctions[i];
  }

  if (type == LocalAssembling_type::TimeNavierStokesMass)
  {
    int num_matrix_denormals = s.mass_matrix.clean_denormals();
    if (num_matrix_denormals > 0)
    {
      Output::root_warn("TNSE", "Mass matrix had ", num_matrix_denormals,
        " non-normal components!");
    }
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::set_arrays(TimeNavierStokes<d>::System_per_grid& s,
                                     std::vector<const FESpace*>& spaces,
                                     std::vector<const FESpace*>& spaces_rhs,
                                     std::vector<const FEFunction*>& functions,
                                     bool first_nonlinear_step)
{
  spaces.resize(2);
  spaces[0] = s.velocity_space.get();
  spaces[1] = s.pressure_space.get();

  spaces_rhs.resize(d + 1);

  for (int i = 0; i < d; ++i)
  {
    spaces_rhs[i] = s.velocity_space.get();
  }

  spaces_rhs[d] = s.pressure_space.get();

  // standard for all methods.

  functions.resize(d + 1);

  if (first_nonlinear_step)
  {
    // at the beginning of a new time step, use the solution(s) from previous 
    // time step(s)
    s.extrapolate_sol.reset();

    if (db["extrapolation_type"].is("constant"))
    {
      s.extrapolate_sol = s.solution_m1;
    }
    else if (db["extrapolation_type"].is("linear"))
    {
      // u_0 = u_{-1} + tau / tau_{-1} (u_{-1} - u_{-2})

      double tau_ratio = TDatabase::TimeDB->TIMESTEPLENGTH / s.tau_m1;

      s.extrapolate_sol = s.solution_m1;
      s.extrapolate_sol.scale(1.0 + tau_ratio);
      s.extrapolate_sol.add_scaled(s.solution_m2, -tau_ratio);
    }
    else
    {
      ErrThrow("Unknown extrapolation type: ",
               db["extrapolation_type"]);
    }

    for (int i = 0; i < d; ++i)
    {
      functions[i] = s.extrapolate_u.GetComponent(i).release();
    }

    functions[d] = &s.extrapolate_p;
  }
  else
  {
    // take the solution from the previous nonlinear iteration
    for (int i = 0; i < d; ++i)
    {
      functions[i] = s.u.GetComponent(i).release();
    }

    functions[d] = &s.p;
  }

  if (db["space_discretization_type"].is("vms_projection"))
  {
    spaces.resize(4);

    // projection space
    spaces[2] = projection_space_.get();

    // label space that indicates local projection space
    spaces[3] = label_for_local_projection_space_.get();

    // append function for the labels of the local projection
    functions.resize(d + 2);
    functions[4] = label_for_local_projection_fefct.get();
  }
  else if (db["space_discretization_type"].is("residual_based_vms")
    || db["space_discretization_type"].is("rbvms_time"))
  {
    functions.resize(2 * d + 1);

    for (int i = 0; i < d; ++i)
    {
      functions[d + 1 + i] = s.u_m1.GetComponent(i).release();
    }
  }

  bool is_linearized = linearized();
  if (is_linearized)
  {
    if ((db["space_discretization_type"].is("galerkin") ||    
         db["space_discretization_type"].is("local_projection") ||
         db["space_discretization_type"].is("jump_stab"))
      && !first_nonlinear_step)
    {
      s.extrapolate_sol.reset();

      if (db["extrapolation_type"].is("constant"))
      {
        s.extrapolate_sol = s.solution_m1;
      }
      else if(db["extrapolation_type"].is("linear"))
      {
        // u_0 = u_{-1} + tau / tau_{-1} (u_{-1} - u_{-2})

        double tau_ratio = TDatabase::TimeDB->TIMESTEPLENGTH / s.tau_m1;

        s.extrapolate_sol = s.solution_m1;
        s.extrapolate_sol.scale(1.0 + tau_ratio);
        s.extrapolate_sol.add_scaled(s.solution_m2, -tau_ratio);
      }
      else
      {
        ErrThrow("Unknown extrapolation type: ",
                 db["extrapolation_type"]);
      }

      for (int i = 0; i < d; ++i)
      {
        delete functions[i];
        functions[i] = s.extrapolate_u.GetComponent(i).release();
      }
    }

    //NOTE: correctly copied from the turbulent branch
    // supg: NOTE: only tested with BDF2 so far

    if (db["space_discretization_type"].is("supg")
      && !db["time_discretization"].is("bdf_two"))
    {
      ErrThrow("supg method is only implemented for BDF2 time stepping scheme");
    }

    if (db["space_discretization_type"].is("supg"))
    {
      if (TDatabase::ParamDB->NSTYPE < 4)
      {
        ErrThrow("SUPG is implemented for equal-order only so far !!!");
      }

      functions.resize(2 * d);

      // BDF2: constant time step only!
      s.extrapolate_sol.reset();
      s.extrapolate_sol = s.solution_m1;
      s.extrapolate_sol.scale(2.);
      s.extrapolate_sol.add_scaled(s.solution_m2, -1.);

      for (int i = 0; i < d; ++i)
      {
        delete functions[i];
        functions[i] = s.extrapolate_u.GetComponent(i).release();
      }

      // combination of previous time solutions for assembling the right-hand
      // side, this is used for the pressure part.
      s.combined_old_sols.reset();

      // copy and scale the solution at previous time step with factor 2
      s.combined_old_sols = s.solution_m1;
      s.combined_old_sols.scale(2.);

      // subtract with right factor the solution at pre-previous solution
      s.combined_old_sols.add_scaled(s.solution_m2, -0.5);

      for (int i = 0; i < d; ++i)
      {
        functions[i + d] = s.comb_old_u.GetComponent(i).release();
      }
    }
  }
  else
  {
    // For the standard methods or symmetric stabilization schemes
    // there is no time derivative involved in combination with the 
    // pressure or velocity as in the supg case. So nothing to do 
    // for those schemes.
    if (db["space_discretization_type"].is("supg"))
    {
      functions.resize(2 * d);

      if (time_stepping_scheme.pre_stage_bdf)
      {
        for(int i = 0; i < d; ++i)
        {
          functions[i + d] = s.u_m1.GetComponent(i).release();
        }
      }
      else
      {
        s.combined_old_sols.reset();

        // copy and scale the solution at previous time step with factor 2
        s.combined_old_sols = s.solution_m1;
        s.combined_old_sols.scale(2.);

        // subtract with right factor the solution at pre-previous solution
        s.combined_old_sols.add_scaled(s.solution_m2, -0.5);

        for(int i = 0; i < d; ++i)
        {
          functions[i + d] = s.comb_old_u.GetComponent(i).release();
        }
      }
    }
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::set_matrices_rhs(
  TimeNavierStokes<d>::System_per_grid& s, LocalAssembling_type type,
  std::vector<SquareMatrixD*>& sqMat, std::vector<MatrixD*>& reMat,
  std::vector<double*>& rhs_array)
{
  sqMat.resize(d * d + 1, nullptr); // maximum number of square matrices (type 14)
  reMat.resize(2 * d, nullptr);

  // right hand side: for NSTYPE: 1,2 and 3, size is 2
  rhs_array.resize(d + 1, nullptr);

  auto blocks = s.matrix.get_blocks_uniquely();
  auto mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
  int nstype = TDatabase::ParamDB->NSTYPE;

  switch (type)
  {
    case LocalAssembling_type::TimeNavierStokesMass:

      if (db["space_discretization_type"].is("residual_based_vms")
        || db["space_discretization_type"].is("rbvms_time"))
      {
        for (int i = 0; i < d; i++)
        {
          for (int j = 0; j < d; j++)
          {
            sqMat[i * d + j] = reinterpret_cast<SquareMatrixD*>(
              mass_blocks.at(i * (d + 1) + j).get());
          }
        }

        for (int i = 0; i < d; i++)
        {
          reMat[i] = reinterpret_cast<MatrixD*>(mass_blocks.at((d + 1) * d + i).get());
          reMat[i + d] = reinterpret_cast<MatrixD*>(mass_blocks.at((d + 1) * i + d + 1).get());
        }
      }
      else
      {
        if (nstype == 1 || nstype == 2)
        {
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(0).get());
        }
        else
        {
          for (int i = 0; i < d; ++i)
          {
            sqMat[i] = reinterpret_cast<SquareMatrixD*>(mass_blocks.at(i * (d + 2)).get());
          }
        }

        if (db["space_discretization_type"].is("vms_projection"))
        {
          sqMat[d * d] = reinterpret_cast<SquareMatrixD*>(matrices_for_turb_mod.at(6).get());
        }
      }
      break; //TimeNavierStokesMass

    case LocalAssembling_type::NavierStokesLinear:
      if (nstype == 1 || nstype == 2)
      {
        sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
      }
      else
      {
        for(int i = 0, j = 0; i < d * d; ++i, ++j)
        {
          if (i % d == 0 && i > 0)
          {
            j++;
          }

          sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
        }
      }
      break;

    case LocalAssembling_type::TimeNavierStokesAll:
    case LocalAssembling_type::NavierStokesAll:
      for (int i = 0; i <= d; ++i)
      {
        rhs_array[i] = s.rhs.block(i);
      }
      s.rhs.reset();

      switch (nstype)
      {
        case 1:
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());

          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[1 + i].get());
          }
          break;

        case 2:
          sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());

          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d + 1 + i].get());
          }

          for(int i = 0; i < d; ++i)
          {
            reMat[i + d] = reinterpret_cast<MatrixD*>(blocks[1 + i].get());
          }
          break;

        case 3:
          for (int i = 0, j = 0; i < d * d; ++i, ++j)
          {
            if (i % d == 0 && i > 0)
            {
              j++;
            }

            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }

          for(int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d + (d + 1) * i].get());
          }
          break;

        case 4:
          for (int i = 0, j = 0; i < d * d; ++i, ++j)
          {
            if (i % d == 0 && i > 0)
            {
              j++;
            }

            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }

          if (db["space_discretization_type"].is("vms_projection"))
          {
            sqMat[d * d + 1] = reinterpret_cast<SquareMatrixD*>(blocks[(d + 1) * (d + 1) - 1].get());
          }

          for(int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d * (d + 1) + i].get());
          }

          for (int i = 0; i < d; ++i)
          {
            reMat[d + i] = reinterpret_cast<MatrixD*>(blocks[d + (d + 1) * i].get());
          }

          if (db["space_discretization_type"].is("vms_projection"))
          {
            reMat.resize(2 * d * 2);

            for (int i = 0; i < 2 * d; i++)
            {
              reMat[2 * d + i] = reinterpret_cast<MatrixD*>(matrices_for_turb_mod.at(i).get());
            }
          }
          break;

        case 14:
          for (int i = 0, j = 0; i < d * d; ++i, ++j)
          {
            if (i % d == 0 && i > 0)
            {
              j++;
            }

            sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[j].get());
          }

          sqMat[d * d] = reinterpret_cast<SquareMatrixD*>(blocks[(d + 1) * (d + 1) - 1].get());

          for (int i = 0; i < d; ++i)
          {
            reMat[i] = reinterpret_cast<MatrixD*>(blocks[d * (d + 1) + i].get());
          }

          for (int i = 0; i < d; ++i)
          {
            reMat[d + i] = reinterpret_cast<MatrixD*>(blocks[i * (d + 1) + d].get());
          }
          break;
      }
      break;//TimeNavierStokesAll

    case LocalAssembling_type::TimeNavierStokesNL:
    case LocalAssembling_type::NavierStokesNL:
    {
      std::vector<std::vector<size_t>> cells;
      for (size_t i = 0; i < d; ++i)
      {
        for (size_t j = 0; j < d; ++j)
        {
          cells.push_back({{i, j}});
        }
      }

      auto blocks = s.matrix.get_blocks_uniquely(cells);

      if (nstype == 1 || nstype == 2)
      {
        sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks[0].get());
      }
      else
      {
        for (int i = 0; i < d * d; ++i)
        {
          sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[i].get());
        }
      }

      if (db["space_discretization_type"].is("residual_based_vms")
        || db["space_discretization_type"].is("rbvms_time"))
      {
        bool t;

        sqMat[d * d] = reinterpret_cast<SquareMatrixD*>(
          s.matrix.get_block(d, d, t).get());

        for (int i = 0; i < d; ++i)
        {
          reMat[i] = reinterpret_cast<MatrixD*>(
            s.matrix.get_block(d, i, t).get());

          reMat[d + i] = reinterpret_cast<MatrixD*>(
            s.matrix.get_block(i, d, t).get());
        }
      }

      if (db["space_discretization_type"].is("vms_projection"))
      {
        reMat.resize(2 * d * 2);
        for (int i = 0; i < d; i++)
        {
          reMat[2 * d + i] = reinterpret_cast<MatrixD*>(matrices_for_turb_mod.at(i).get());
        }
      }

      if (db["space_discretization_type"].is("supg") && nstype == 4)
      {
        auto blocks = s.matrix.get_blocks_uniquely();

        for (int i = 0; i < d; ++i)
        {
          reMat[d + i] = reinterpret_cast<MatrixD*>(blocks[d + (d + 1) * i].get());
        }

        for (int i = 0; i < d + 1; ++i)
        {
          rhs_array[i] = s.rhs.block(i);
        }

        s.rhs.reset();
      }
      break;//TimeNavierStokesNL
    }

    case LocalAssembling_type::TimeNavierStokesRhs:
    case LocalAssembling_type::TimeNavierStokesExplNL:
      // no matrices to be assembled
      for (int i = 0; i < d + 1; ++i)
      {
        rhs_array[i] = s.rhs.block(i);
      }

      s.rhs.reset();
      break;//TimeNavierStokesRhs

    default:
      ErrThrow("This LocalAssembling_type (", type, ") is not supported here.");
  }

  for (auto* sm : sqMat)
  {
    if (sm != nullptr)
    {
      sm->reset();
    }
  }

  for (auto rm : reMat)
  {
    if (rm != nullptr)
    {
      rm->reset();
    }
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::restrict_function()
{
  // assembling requires an approximate velocity solution on every grid
  for (int block = 0; block < d; ++block)
  {
    std::vector<const FESpace*> spaces;
    std::vector<double*> u_entries;
    std::vector<size_t> u_ns_dofs;

    for (auto &s: systems)
    {
      spaces.push_back(s.velocity_space.get());
      u_entries.push_back(s.solution.block(block));
      u_ns_dofs.push_back(s.solution.length(block));
    }

    GridTransfer::RestrictFunctionRepeatedly(spaces, u_entries, u_ns_dofs);
  }
}

/* ************************************************************************** */
template <int d>
bool TimeNavierStokes<d>::stop_it(unsigned int it_counter)
{
  bool iteration_failed;
  return stop_it(it_counter, false, iteration_failed);
}

template <int d>
bool TimeNavierStokes<d>::stop_it(unsigned int it_counter, bool &iteration_failed)
{
  return stop_it(it_counter, true, iteration_failed);
}

template <int d>
bool TimeNavierStokes<d>::stop_it(unsigned int it_counter, bool reset_on_failure, bool &iteration_failed)
{
  compute_residuals();

  const double norm_of_residual = this->get_full_residual();

  if (it_counter == 0)
  {
    initial_residual = norm_of_residual;
  }

  // check if minimum number of iterations was performed already
  size_t min_it = db["nonlinloop_minit"];
  if (it_counter < min_it)
  {
    iteration_failed = false;
    return false;
  }

  // hold the residual from 10 iterations ago
  const double very_old_norm_of_residual = old_residuals.front().fullResidual;
  size_t max_it = db["nonlinloop_maxit"];
  double convergence_speed = db["nonlinloop_slowfactor"];

  double limit = db["nonlinloop_epsilon"];
  double max_limit = db["nonlinloop_max_epsilon"];

  bool hard_limit = db["nonlinloop_hard_epsilon"];
  bool linearized_scheme = this->linearized();

  // check various stopping criteria

  if (db["nonlinloop_scale_epsilon_with_size"])
  {
    limit *= std::sqrt(this->get_size());
    Output::root_info("TNSE", "stopping tolerance for nonlinear iteration ", limit);
  }

  bool epsilon_reached = norm_of_residual <= limit || old_residuals.back().fullResidualMax <= max_limit;

  bool final_iteration = epsilon_reached
    || it_counter == max_it
    || (linearized_scheme && it_counter > 0);

  iteration_failed = ((norm_of_residual >= convergence_speed * very_old_norm_of_residual)
    || (!epsilon_reached && it_counter >= max_it)) && hard_limit;

  if (iteration_failed)
  {
    Output::root_info<3>("TNSE", "Iteration failed: ",
      "\nresidual: ", norm_of_residual, " / ", convergence_speed * very_old_norm_of_residual,
      " / ", limit,
      "\niterations: ", it_counter, " / ", max_it);
  }

  if (!iteration_failed && final_iteration)
  {
    // check kinetic energy and velocity criteria, if applicable

    double threshold = db.try_get_value("nonlinloop_slow_difference_threshold", 0.0);
    double ratio_energy = db.try_get_value("nonlinloop_slow_difference_ratio_energy", -1.0);
    double ratio_difference = db.try_get_value("nonlinloop_slow_difference_ratio_difference", -1.0);

    double velocity_threshold
      = db.try_get_value("nonlinloop_slow_difference_velocity_threshold", 0.0);
    double ratio_max_velocity
      = db.try_get_value("nonlinloop_slow_difference_ratio_max_velocity", -1.0);

    if (ratio_max_velocity > 0.0)
    {
      // MPI: velocity is consistency level 3, so no need to check masters

      const FEVectFunct& velocity_old = systems.front().u_m1;
      const FEVectFunct& velocity_new = systems.front().u;

      double max_velocity_old = 0.0;
      double max_velocity_new = 0.0;

      int n_dof = systems.front().velocity_space->get_n_dof();

      const double* vo_data = velocity_old.GetValues();
      const double* vn_data = velocity_new.GetValues();

      for (int i = 0; i < n_dof; i++)
      {
        double vo = 0.0;
        double vn = 0.0;

        for (int j = 0; j < d; j++)
        {
          int k = i + j * n_dof;

          vo += vo_data[k] * vo_data[k];
          vn += vn_data[k] * vn_data[k];
        }

        if (vo > max_velocity_old)
        {
          max_velocity_old = vo;
        }

        if (vn > max_velocity_new)
        {
          max_velocity_new = vn;
        }
      }

      max_velocity_old = std::sqrt(max_velocity_old);
      max_velocity_new = std::sqrt(max_velocity_new);

#ifdef _MPI
      MPI_Allreduce(MPI_IN_PLACE, &max_velocity_old, 1, MPI_DOUBLE,
        MPI_MAX, MPI_COMM_WORLD);

      MPI_Allreduce(MPI_IN_PLACE, &max_velocity_new, 1, MPI_DOUBLE,
        MPI_MAX, MPI_COMM_WORLD);
#endif

      if (max_velocity_old > velocity_threshold)
      {
        if (std::abs(max_velocity_new - max_velocity_old) > ratio_max_velocity * max_velocity_old)
        {
          Output::root_info<3>("TNSE", "Iteration failed: maximum velocity ratio exceeded.");
          iteration_failed = true;
        }
      }

      Output::root_info<3>("TNSE", "Old max velocity: ", max_velocity_old,
        "\nNew max velocity: ", max_velocity_new,
        "\nDifference ratio: ", (max_velocity_new - max_velocity_old) / max_velocity_old);
    }

    if (ratio_difference > 0.0 || ratio_energy > 0.0)
    {
#ifdef __3D__
      TAuxParam3D aux;
      MultiIndex3D nsValueOnly[1] = { MultiIndex3D::D000 };
#else
      TAuxParam2D aux;
      MultiIndex2D nsValueOnly[1] = { MultiIndex2D::D00 };
#endif

      const FESpace *velocity_space = systems.front().velocity_space.get();
      int n_u_dof = velocity_space->get_n_dof();

      std::vector<double> u_diff(n_u_dof);

      FEFunction velocity_difference(systems.front().velocity_space, "u_diff",
                                     u_diff.data());

      const FEVectFunct& velocity_old = systems.front().u_m1;
      const FEVectFunct& velocity_new = systems.front().u;

      double kinetic_energy_old = 0.0;
      double kinetic_energy_new = 0.0;
      double kinetic_energy_difference = 0.0;

      for (int i = 0; i < d; i++)
      {
        double tmp_ko;
        double tmp_kn;
        double tmp_kd;

        velocity_old.GetComponent(i)->GetErrors(
#ifdef __2D__
          unknown_solution_2d,
#else
          unknown_solution_3d,
#endif
          1, nsValueOnly, 1, L2Error, nullptr, &aux, 1, &velocity_space,
          &tmp_ko);

        velocity_new.GetComponent(i)->GetErrors(
#ifdef __2D__
          unknown_solution_2d,
#else
          unknown_solution_3d,
#endif
          1, nsValueOnly, 1, L2Error, nullptr, &aux, 1, &velocity_space,
          &tmp_kn);

        const double* old_values = velocity_old.GetComponent(i)->GetValues();
        const double* new_values = velocity_new.GetComponent(i)->GetValues();

        for (int j = 0; j < n_u_dof; j++)
        {
          u_diff[j] = new_values[j] - old_values[j];
        }

        velocity_difference.GetErrors(
#ifdef __2D__
          unknown_solution_2d,
#else
          unknown_solution_3d,
#endif
          1, nsValueOnly, 1, L2Error, nullptr, &aux, 1, &velocity_space,
          &tmp_kd);

#ifdef _MPI
        MPI_Allreduce(MPI_IN_PLACE, &tmp_ko, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &tmp_kn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &tmp_kd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // GetErrors returns squared L2 norms in the MPI case
        kinetic_energy_old += tmp_ko;
        kinetic_energy_new += tmp_kn;
        kinetic_energy_difference += tmp_kd;
#else
        // GetErrors returns L2 norms (square roots) in the non-MPI case
        kinetic_energy_old += tmp_ko * tmp_ko;
        kinetic_energy_new += tmp_kn * tmp_kn;
        kinetic_energy_difference += tmp_kd * tmp_kd;
#endif
      }

      // half the total squared velocity magnitude
      kinetic_energy_old *= 0.5;
      kinetic_energy_new *= 0.5;
      kinetic_energy_difference *= 0.5;

      Output::root_info<3>("TNSE", "Old kinetic energy: ", kinetic_energy_old,
        "\nNew kinetic energy: ", kinetic_energy_new,
        "\nDifference kinetic energy: ", kinetic_energy_difference,
        "\nKinetic energy ratio: ", (kinetic_energy_new - kinetic_energy_old) / kinetic_energy_old,
        "\nKinetic energy diff ratio: ", std::sqrt(kinetic_energy_difference / kinetic_energy_old));

      if (kinetic_energy_old > threshold)
      {
        if (ratio_energy > 0.0)
        {
          if (std::abs(kinetic_energy_new - kinetic_energy_old) > ratio_energy * kinetic_energy_old)
          {
            Output::root_info<3>("TNSE", "Iteration failed: kinetic energy ratio exceeded.");
            iteration_failed = true;
          }
        }

        if (ratio_difference > 0.0)
        {
          if (kinetic_energy_difference > ratio_difference * ratio_difference * kinetic_energy_old)
          {
            Output::root_info<3>("TNSE", "Iteration failed: difference kinetic energy ratio exceeded.");
            iteration_failed = true;
          }
        }
      }
    }
  }

  if (iteration_failed && reset_on_failure)
  {
    Output::root_warn("TNSE", "Slow convergence or velocity changing too quickly! res = ", norm_of_residual, ", ", it_counter, " iterations");

    // reset old rhs
    base_rhs = base_rhs_m1;

    // reset solution
    for (System_per_grid& s: this->systems)
    {
      s.solution = s.solution_m1;
    }

    return true;
  }
  else if (final_iteration)
  {
    const ParameterDatabase e_db = example.get_database();
    int n_windkessel_bd = e_db.try_get_value("n_windkessel_bd", 0);

    // advance windkessel bc
    if (n_windkessel_bd)
    {
      std::vector<size_t> windkessel_id = e_db["windkessel_id"];

      size_t r = std::min(windkessel_id.size(), (size_t)n_windkessel_bd);

      for (size_t k = 0; k < r; k++)
      {
        wk_boundary[k].advance();
      }
    }

    // advance old rhs
    base_rhs_m2 = base_rhs_m1;
    base_rhs_m1 = base_rhs;

    // advance solution
    for (System_per_grid& s: this->systems)
    {
      s.solution_m2 = s.solution_m1;
      s.solution_m1 = s.solution;

      s.tau_m2 = s.tau_m1;
      s.tau_m1 = TDatabase::TimeDB->TIMESTEPLENGTH;

      s.persistent_data->Advance();
    }

    adjust_pressure();
    return true;
  }
  else
  {
    return false;
  }
}

/* ************************************************************************** */
template <int d>
Residuals TimeNavierStokes<d>::compute_residuals_unstored()
{
  System_per_grid& s = this->systems.front();
  unsigned int number_u_Dof = s.solution.length(0);

#ifdef _MPI
    // MPI: put solution in consistency level 3
    auto comms = s.matrix.get_communicators();
    for (size_t bl = 0; bl < comms.size() ;++bl)
    {
      comms[bl]->consistency_update(s.solution.block(bl), 3);
    }
#endif

  // copy rhs to defect and compute defect
  this->defect = rhs_from_time_disc;
  s.matrix.apply_scaled_add(s.solution, defect, -1.);

  if (s.matrix.pressure_projection_enabled())
  {
    FEFunction defect_fctn(s.pressure_space, "p_defect",
                           &defect[d * number_u_Dof]);
    defect_fctn.project_into_L20();
  }

  // This is the calculation of the residual, given the defect.
  std::vector<unsigned int> velocity_blocks(d, 0);
  std::iota(std::begin(velocity_blocks), std::end(velocity_blocks), 0);

#ifdef _MPI
  std::vector<const TParFECommunicator3D*> velocity_comms(d, nullptr);
  for(int i = 0; i < d; ++i)
  {
    velocity_comms[i] = comms[i];
  }

  double momentum_residual_square = defect.norm(velocity_blocks, velocity_comms);
  double mass_residual_square = defect.norm({d}, {comms[d]});

  double momentum_residual_max = defect.norm_infty(velocity_blocks, velocity_comms);
  double mass_residual_max = defect.norm_infty({d}, {comms[d]});

#else
  double momentum_residual_square = defect.norm(velocity_blocks);
  double mass_residual_square = defect.norm({d});

  double momentum_residual_max = defect.norm_infty(velocity_blocks);
  double mass_residual_max = defect.norm_infty({d});
#endif

  // the struct 'Residuals' takes squares:
  momentum_residual_square *= momentum_residual_square;
  mass_residual_square *= mass_residual_square;

  return Residuals(momentum_residual_square, mass_residual_square,
    momentum_residual_max, mass_residual_max);
}

/* ************************************************************************** */

template <int d>
void TimeNavierStokes<d>::compute_residuals()
{
  old_residuals.add(compute_residuals_unstored());
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::solve()
{
  ++solve_count;

  bool diagnostics = db["nonlinloop_diagnostics"];
  std::ofstream diag_f;

#ifdef _MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_rank == 0)
  {
#endif
    if (diagnostics)
    {
      diag_f = std::ofstream(outputWriter.get_base_name() + "_nonlinloop.csv", std::ios::app);
    }
#ifdef _MPI
  }
#endif

  Residuals initial_res = compute_residuals_unstored();

  Output::root_info("TNSE", "Initial ", initial_res);

  System_per_grid& s = this->systems.front();

  // store previous solution for damping, it is a pointer so that we can avoid
  // the copy in case of no damping
  double damping = this->db["nonlinloop_damping_factor"];

  bool damping_auto = this->db.try_get_value("nonlinloop_damping_auto", false);
  bool apply_linear_auto_damping = damping_auto;
  bool apply_nonlinear_auto_damping = damping_auto && !linearized();

  if (apply_nonlinear_auto_damping)
  {
    int skip = db["nonlinloop_damping_auto_skip"];
    int maxit = db["nonlinloop_maxit"];

    if (solve_count <= skip
      && solve_count <= maxit)
    {
      apply_nonlinear_auto_damping = false;
    }
  }

  bool apply_damping = apply_linear_auto_damping
    || apply_nonlinear_auto_damping
    || (damping != 1.0);

  std::shared_ptr<BlockVector> old_solution(nullptr);

  if (apply_damping)
  {
    old_solution = std::make_shared<BlockVector>(s.solution);
  }

  bool direct_solver = solver.get_db()["solver_type"].is("direct");

  bool fully_explicit =
    this->db["time_discretization_nonlinear_term"].is("fully_explicit");

  int num_rhs_denormals = rhs_from_time_disc.clean_denormals();
  if (num_rhs_denormals > 0)
  {
    Output::root_warn("TNSE", "RHS had ", num_rhs_denormals, " non-normal components!");
  }

  int num_matrix_denormals = s.matrix.clean_denormals();
  if (num_matrix_denormals > 0)
  {
    Output::root_warn("TNSE", "Matrix had ", num_matrix_denormals, " non-normal components!");
  }

  bool bdf2 = time_stepping_scheme.n_old_solutions() > 1;

#ifdef _MPI
  if (direct_solver)
  {
    MumpsWrapper mumps_wrapper(s.matrix);

    if (solver.get_db()["mumps_distributed"])
    {
      mumps_wrapper.enable_distributed_solution(
        solver.get_db()["max_n_iterations"], solver.get_db()["residual_tolerance"]);
    }

    mumps_wrapper.solve(rhs_from_time_disc, s.solution);
  }
  else
#endif
  {
    if (direct_solver && fully_explicit
       && this->time_stepping_scheme.current_step_ > (bdf2 ? 2 : 1))
    {
      // reuse the already computed factorisation
      solver.solve(rhs_from_time_disc, s.solution);
    }
    else
    {
      solver.solve(s.matrix, rhs_from_time_disc, s.solution);
    }
  }

#ifdef _MPI
  auto comms = s.matrix.get_communicators();
  for (size_t bl = 0; bl < comms.size(); ++bl)
  {
    comms[bl]->consistency_update(s.solution.block(bl), 3);
  }
#endif

  int num_sol_denormals = s.solution.clean_denormals();
  if (num_sol_denormals > 0)
  {
    Output::root_warn("TNSE", "Solution had ", num_sol_denormals, " non-normal components!");
  }

  if (apply_damping)
  {
    if (damping_auto)
    {
      BlockVector new_solution(s.solution);

      double best_damping = 1.0;
      double base_residual = 0.0;
      double best_residual = 0.0;

      if (apply_linear_auto_damping)
      {
        std::array<double, 3> res;

        for (int i = 0; i < 3; i++)
        {
          s.solution = new_solution;
          s.solution.scale((double)(i - 1));
          s.solution.add_scaled(*old_solution, (double)(2 - i));

          if (s.matrix.pressure_projection_enabled())
          {
            s.p.project_into_L20();
          }

          res[i] = compute_residuals_unstored();
          res[i] *= res[i];
        }

        // res = at² + bt + c
        // res[0] = a - b + c
        // res[1] = c
        // res[2] = a + b + c
        // a + b = res[2] - res[1]
        // a - b = res[0] - res[1]
        // 2a = res[0] - 2 res[1] + res[2]

        double a = 0.5 * (res[0] + res[2]) - res[1];
        double b = res[2] - res[1] - a;
        double c = res[1];

        double best_damping;

        if (a > 0.0)
        {
          // 0 = res' = 2at + b
          // t = - b / 2a

          best_damping = -0.5 * b / a;
        }
        else
        {
          Output::root_warn("TNSE", "Non-pd auto damping...?");
          best_damping = 1.0;
        }

        best_residual = std::sqrt(c + best_damping * (b + a * best_damping));
        base_residual = std::sqrt(res[2]);

        if (true)
        {
          Output::root_info("TNSE", "Linear auto damping: ", best_damping);

          new_solution.scale(best_damping);
          new_solution.add_scaled(*old_solution, 1.0 - best_damping);

          base_residual = best_residual;
          best_damping = 1.0;
        }
      }

      if (apply_nonlinear_auto_damping)
      {
        int reassembly_step = -1;

        // first: simple logarithmically spaced search from 0.5 to 2 with
        // (2 * precision + 1) steps. if 0.5 is the smallest, also searches
        // from 0.0625 to 0.5 with the same logarithmic step length.

        int precision = db["nonlinloop_damping_auto_precision"];
        bool binary_search = this->db.try_get_value(
          "nonlinloop_damping_auto_binary_search", false);

        precision = std::max(1, precision);

        double inv_precision = 1.0 / precision;

        bool first = true;
        int best_i = -1;

        BlockVector tmp_rhs(s.rhs);
        BlockVector tmp_base_rhs(base_rhs);

        for (int i = -precision; i <= precision; i++)
        {
          double damping_i = std::pow(2.0, inv_precision * (double)i);

          s.solution = new_solution;
          s.solution.scale(damping_i);
          s.solution.add_scaled(*old_solution, 1.0 - damping_i);

          if (s.matrix.pressure_projection_enabled())
          {
            s.p.project_into_L20();
          }

          base_rhs = tmp_base_rhs;
          s.rhs = tmp_rhs;

          // reassemble, or we'll be optimizing the wrong residual
          assemble_matrices_rhs(reassembly_step);

          double res_i = compute_residuals_unstored();

          Output::root_info<5>("TNSE", "Damping ", damping_i, ": res = ", res_i);

          if (i == 0)
          {
            base_residual = res_i;
          }

          if (first || res_i < best_residual)
          {
            best_damping = damping_i;
            best_i = i;
            best_residual = res_i;

            first = false;
          }
        }

        if (best_i == -precision)
        {
          for (int i = -4 * precision; i < -precision; i++)
          {
            double damping_i = std::pow(2.0, inv_precision * (double)i);

            s.solution = new_solution;
            s.solution.scale(damping_i);
            s.solution.add_scaled(*old_solution, 1.0 - damping_i);

            if (s.matrix.pressure_projection_enabled())
            {
              s.p.project_into_L20();
            }

            base_rhs = tmp_base_rhs;
            s.rhs = tmp_rhs;

            // reassemble, or we'll be optimizing the wrong residual
            assemble_matrices_rhs(reassembly_step);

            double res_i = compute_residuals_unstored();

            Output::root_info<5>("TNSE", "Damping ", damping_i, ": res = ", res_i);

            if (res_i < best_residual)
            {
              best_damping = damping_i;
              best_i = i;
              best_residual = res_i;
            }
          }
        }

        if (binary_search)
        {
          // binary search on a not necessarily monotone function - may not
          // converge to the correct minimum.
          // - initial interval: logarithmically symmetric around the previous
          //   search's result, just short of its neighbors in the previous
          //   search.
          // - on each iteration, replaces the end of the interval with the
          //   larger residual with the interval's center.

          double prec_mul = (double)precision / (double)(precision + 1);

          double d_a = std::pow(2.0, inv_precision * (double)(best_i - prec_mul));
          double d_b = std::pow(2.0, inv_precision * (double)(best_i + prec_mul));

          s.solution = new_solution;
          s.solution.scale(d_a);
          s.solution.add_scaled(*old_solution, 1.0 - d_a);

          if (s.matrix.pressure_projection_enabled())
          {
            s.p.project_into_L20();
          }

          base_rhs = tmp_base_rhs;
          s.rhs = tmp_rhs;

          // reassemble, or we'll be optimizing the wrong residual
          assemble_matrices_rhs(reassembly_step);

          double res_a = compute_residuals_unstored();

          s.solution = new_solution;
          s.solution.scale(d_b);
          s.solution.add_scaled(*old_solution, 1.0 - d_b);

          if (s.matrix.pressure_projection_enabled())
          {
            s.p.project_into_L20();
          }

          base_rhs = tmp_base_rhs;
          s.rhs = tmp_rhs;

          // reassemble, or we'll be optimizing the wrong residual
          assemble_matrices_rhs(reassembly_step);

          double res_b = compute_residuals_unstored();

          while (d_b - d_a > 1e-6)
          {
            double d_c = 0.5 * (d_a + d_b);

            s.solution = new_solution;
            s.solution.scale(d_c);
            s.solution.add_scaled(*old_solution, 1.0 - d_c);

            if (s.matrix.pressure_projection_enabled())
            {
              s.p.project_into_L20();
            }

            base_rhs = tmp_base_rhs;
            s.rhs = tmp_rhs;

            // reassemble, or we'll be optimizing the wrong residual
            assemble_matrices_rhs(reassembly_step);

            double res_c = compute_residuals_unstored();

            if (res_c < best_residual)
            {
              best_residual = res_c;
              best_damping = d_c;
            }

            if (res_a < res_b)
            {
              d_b = d_c;
              res_b = res_c;
            }
            else
            {
              d_a = d_c;
              res_a = res_c;
            }
          }
        }

        base_rhs = tmp_base_rhs;
        s.rhs = tmp_rhs;
      }

      s.solution = new_solution;
      s.solution.scale(best_damping);
      s.solution.add_scaled(*old_solution, 1.0 - best_damping);

      if (s.matrix.pressure_projection_enabled())
      {
        s.p.project_into_L20();
      }

      Output::root_info("TNSE", "Best damping: ", best_damping,
          ", residual: ", best_residual,
          ", at 1.0:  ", base_residual);

      damping = best_damping;
    }
    else
    {
      s.solution.scale(damping);
      s.solution.add_scaled(*old_solution, 1.0 - damping);
    }
  }

  if (!damping_auto && s.matrix.pressure_projection_enabled())
  {
    s.p.project_into_L20();
  }

  if (diagnostics)
  {
    BlockVector tmp_rhs(s.rhs);
    BlockVector tmp_base_rhs(base_rhs);

    // reassemble to get correct residual
    assemble_matrices_rhs(-1);

    double final_res = compute_residuals_unstored();

    base_rhs = tmp_base_rhs;
    s.rhs = tmp_rhs;

#ifdef _MPI
    if (mpi_rank == 0)
    {
#endif
      diag_f.precision(18);

      diag_f << TDatabase::TimeDB->CURRENTTIME << ", ";
      diag_f << TDatabase::TimeDB->CURRENTTIMESTEPLENGTH << ", ";
      diag_f << damping << ", ";
      diag_f << final_res << "\n";
#ifdef _MPI
    }
#endif
  }
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::output()
{
  bool i_am_root = true;

#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == 0);
#endif

  System_per_grid& s = this->systems.front();
  std::array<std::unique_ptr<FEFunction>, d> velocity_components;

  for (int i = 0; i < d; ++i)
  {
    velocity_components[i] = s.u.GetComponent(i);
  }

  if (db["output_compute_minmax"])
  {
    for(int i = 0; i < d; ++i)
    {
      velocity_components[i]->PrintMinMax();
    }

    s.p.PrintMinMax();
  }

  double t = time_stepping_scheme.current_time_;

  if (db["output_compute_errors"])
  {
#ifdef __3D__
    bool compute_tke = db["space_discretization_type"].is("smagorinsky")
      || db["space_discretization_type"].is("residual_based_vms")
      || db["space_discretization_type"].is("rbvms_time");
#else
    bool compute_tke = false;
#endif

    std::vector<std::array<double, 5>> computed_errors(d + 1);

#ifdef __3D__
    TAuxParam3D aux;
    MultiIndex3D nsAllDerivs[d + 1] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                        MultiIndex3D::D010, MultiIndex3D::D001 };
#else
    TAuxParam2D aux;
    MultiIndex2D nsAllDerivs[d + 1] = { MultiIndex2D::D00, MultiIndex2D::D10,
                                        MultiIndex2D::D01 };
#endif

    const FESpace *velocity_space = this->systems.front().velocity_space.get();
    const FESpace *pressure_space = this->systems.front().pressure_space.get();

    double tau = TDatabase::TimeDB->TIMESTEPLENGTH;

    double kinetic_energy = 0.0;
    double sgs_energy = 0.0;
    double max_velocity = 0.0;

    int n_u_dof = systems.front().velocity_space->get_n_dof();
    const double* u_data = systems.front().u.GetValues();

    for (int i = 0; i < n_u_dof; i++)
    {
      double v = 0.0;

      for (int j = 0; j < d; j++)
      {
        int k = i + j * n_u_dof;

        v += u_data[k] * u_data[k];
      }

      if (v > max_velocity)
      {
        max_velocity = v;
      }
    }

    max_velocity = std::sqrt(max_velocity);

#ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, &max_velocity, 1, MPI_DOUBLE,
      MPI_MAX, MPI_COMM_WORLD);
#endif

    // errors in the velocity components
    for (int i = 0; i < d; ++i)
    {
      auto ui = velocity_components[i].get();
      ui->GetErrors(example.get_exact(i), d + 1, nsAllDerivs, d, L2H1Errors,
                    nullptr, &aux, 1, &velocity_space,
                    computed_errors[i].data());

      // computing the L^2 integral of the velocity to be able to compute the 
      // kinetic energy, we use the entries in `computed_errors` which are later
      // used for the pressure.
      ui->GetErrors(
#ifdef __2D__
        unknown_solution_2d,
#else
        unknown_solution_3d,
#endif
       d + 1, nsAllDerivs, d, L2H1Errors, nullptr, &aux, 1, &velocity_space,
       computed_errors[d].data());

#ifdef _MPI
      kinetic_energy += computed_errors[d][0];
#else
      kinetic_energy += 0.5 * computed_errors[d][0] * computed_errors[d][0];
#endif
    }

#ifdef __3D__
    if (compute_tke)
    {
      if (db["space_discretization_type"].is("smagorinsky"))
      {
        sgs_energy = ComputeTurbulentKineticEnergySmagorinsky(systems.front().u,
          db["turbulent_kinetic_energy_scale"]);
      }
      else if (db["space_discretization_type"].is("residual_based_vms"))
      {
        sgs_energy = ComputeTurbulentKineticEnergyRBVMS(systems.front().u,
          systems.front().u_m1, systems.front().p, example.get_nu(),
          RBVMS_Settings(db));
      }
      else if (db["space_discretization_type"].is("rbvms_time"))
      {
        sgs_energy = ComputeTurbulentKineticEnergyRBVMSTime(systems.front().u,
          systems.front().persistent_data);
      }
      else
      {
        Output::root_warn("TNSE", "Cannot compute turbulent kinetic energy "
          "for space_discretization_type = ", db["space_discretization_type"]);
      }
    }
#endif

    // error in divergence
    double div = s.u.GetL2NormDivergenceError(example.get_exact(0),
                                              example.get_exact(1)
#ifdef __3D__
                                            , example.get_exact(2)
#endif
                                             );
    // errors in pressure
    s.p.GetErrors(example.get_exact(d), d + 1, nsAllDerivs, d, L2H1Errors,
                  nullptr, &aux, 1, &pressure_space, computed_errors[d].data());

#ifdef _MPI
    // L2 and H1 (i.e. 2) errors for d velocity components and the pressure, 
    // additionally, the divergence and kinetic energy (+2)

    constexpr int n_send = 2 * (d + 1) + 2;

    double err_red[n_send]; //memory for global (across all processes) error
    double err_send[n_send]; //fill send buffer
    err_send[0] = computed_errors[0][0]; // L2 error of u0
    err_send[1] = computed_errors[0][1]; // H1 error of u0
    err_send[2] = computed_errors[1][0]; // L2 error of u1
    err_send[3] = computed_errors[1][1]; // H1 error of u1
    err_send[4] = computed_errors[2][0]; // L2 error of u2
    err_send[5] = computed_errors[2][1]; // H1 error of u2
    err_send[2 * d] = div*div;
    err_send[2 * d + 1] = computed_errors[d][0]; // L2 error of p
    err_send[2 * d + 2] = computed_errors[d][1]; // H1 error of p
    err_send[2 * (d + 1) + 1] = kinetic_energy;

    MPI_Allreduce(err_send, err_red, n_send, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    for (int i = 0; i < n_send; i++)
    {
      // MPI: sqrt was skipped in GetErrors function - do it here globally!
      err_red[i] = std::sqrt(err_red[i]);
    }

    //fill the reduced errors back where they belong
    computed_errors[0][0] = err_red[0];
    computed_errors[0][1] = err_red[1];
    computed_errors[1][0] = err_red[2];
    computed_errors[1][1] = err_red[3];
    computed_errors[2][0] = err_red[4];
    computed_errors[2][1] = err_red[5];

    div = err_red[2 * d];

    computed_errors[d][0] = err_red[2 * d + 1];
    computed_errors[d][1] = err_red[2 * d + 2];

    kinetic_energy = 0.5 * err_red[2 * (d + 1) + 1] * err_red[2 * (d + 1) + 1];
#endif

    double l2_u = 0;
    double h1_u = 0;
    for(int i = 0; i < d; ++i)
    {
      l2_u += computed_errors[i][0] * computed_errors[i][0];
      h1_u += computed_errors[i][1] * computed_errors[i][1];
    }

    l2_u = std::sqrt(l2_u);
    h1_u = std::sqrt(h1_u);
    double l2_p = computed_errors[d][0];
    double h1_p = computed_errors[d][1];

    double old_l2_u = errors[0];
    double old_div = errors[1];
    double old_h1_u = errors[2];
    double old_l2_p = errors[3];
    double old_h1_p = errors[4];

    errors[0] = l2_u;
    errors[1] = div;
    errors[2] = h1_u;
    errors[3] = l2_p;
    errors[4] = h1_p;

    errors[5] += (l2_u * l2_u + old_l2_u * old_l2_u) * tau * 0.5;
    errors[6] += (div * div + old_div * old_div) * tau * 0.5;
    errors[7] += (h1_u * h1_u + old_h1_u * old_h1_u) * tau * 0.5;
    errors[8] += (l2_p * l2_p + old_l2_p * old_l2_p) * tau * 0.5;
    errors[9] += (h1_p * h1_p + old_h1_p * old_h1_p) * tau * 0.5;

    if (i_am_root)
    {
      using namespace Output;
      stat("TimeNavierStokes", "Measured errors");
      dash(t, " L2(u)              : ", setprecision(14), errors[0]);
      dash(t, " L2(div(u))         : ", setprecision(14), errors[1]);
      dash(t, " H1-semi(u)         : ", setprecision(14), errors[2]);
      dash(t, " L2(p)              : ", setprecision(14), errors[3]);
      dash(t, " H1-semi(p)         : ", setprecision(14), errors[4]);
      dash(t, " L2(0,t,L2(u))      : ", setprecision(14), std::sqrt(errors[5]));
      dash(t, " L2(0,t,L2(div(u))) : ", setprecision(14), std::sqrt(errors[6]));
      dash(t, " L2(0,t,H1-semi(u)) : ", setprecision(14), std::sqrt(errors[7]));
      dash(t, " L2(0,t,L2(p))      : ", setprecision(14), std::sqrt(errors[8]));
      dash(t, " L2(0,t,H1-semi(p)) : ", setprecision(14), std::sqrt(errors[9]));
      dash(t, " kinetic energy     : ", setprecision(14), kinetic_energy);
      dash(t, " max velocity       : ", setprecision(14), max_velocity);
      if (compute_tke)
      {
        dash(t, " sgs kinetic energy : ", setprecision(14), sgs_energy);
      }

      last_tke = sgs_energy;
      last_kinetic_energy = kinetic_energy;
      last_max_velocity = max_velocity;
    }
  }

  if (db["output_compute_time_average"])
  {
    this->time_averaging();
  }

  if (db["output_along_line"])
  {
    // fill a vector with all fe functions to be evaluated using this->Lines
    std::vector<const FEFunction*> fe_functions;

    for (int i = 0; i < d; ++i)
    {
      fe_functions.push_back(velocity_components[i].get());
    }

    fe_functions.push_back(&s.p);

    if (db["output_compute_time_average"])
    {
      for (int i = 0; i < d; ++i)
      {
        fe_functions.push_back(s.u_time_avg.GetComponent(i).release());
      }

      fe_functions.push_back(&s.p_time_avg);
    }

    Lines.write_fe_values(fe_functions, t);

    if (db["output_compute_time_average"])
    {
      for(int i = 0; i < d; ++i)
      {
        delete fe_functions[d + 1 + i];
      }
    }
  }

  if (db["compute_time_derivative"])
  {
    get_time_derivative(time_derivative);
  }

#ifdef __2D__

  int n = s.solution.length(0);
  double *sol = s.solution.get_entries();
  const FESpace *space = s.velocity_space.get();

  if (!db["space_discretization_type"].is("local_projection"))
  {
    StreamFunction(space, sol, sol + n, stream_function_space.get(), psi.data());
  }

  if (db["compute_vorticity_divergence"])
  {
    ComputeVorticityDivergence(velocity_components[0].get(),
                  velocity_components[1].get(), vorticity_space.get(),
                  vorticity_funct->GetValues(), divergence->GetValues());
  }

#else

  if (db["compute_wall_shear_stress"])
  {
    Output::root_info("TNSE", "Computing wall shear stress...");

    if (wall_shear_stress_funct == nullptr)
    {
      Output::root_warn("TNSE", "Wall shear stress computation not set up.");
    }
    else if (db["wall_shear_stress_mode"].is("face"))
    {
      ComputeWallShearStressFace(example.get_nu(), get_velocity(),
        *wall_shear_stress_funct, ViscositySettings(db));
    }
    else if (db["wall_shear_stress_mode"].is("vertex"))
    {
      ComputeWallShearStress(example.get_nu(), get_velocity(),
        *wall_shear_stress_funct, ViscositySettings(db));
    }
    else
    {
      ErrThrow("Unknown wall shear stress computation mode '",
        db["wall_shear_stress_mode"], "'.");
    }
  }

  if (db["compute_turbulent_kinetic_energy"])
  {
    Output::root_info("TNSE", "Computing turbulent kinetic energy...");

    if (tke_funct == nullptr)
    {
      Output::root_warn("TNSE", "Turbulent kinetic energy computation not set up.");
    }
    else if (db["turbulent_kinetic_energy_mode"].is("cell"))
    {
      ComputeTurbulentKineticEnergyCell(get_velocity(), *tke_funct,
        db["turbulent_kinetic_energy_scale"]);
    }
    else if (db["turbulent_kinetic_energy_mode"].is("vertex"))
    {
      ComputeTurbulentKineticEnergy(get_velocity(), *tke_funct,
        db["turbulent_kinetic_energy_scale"]);
    }
    else
    {
      ErrThrow("Unknown turbulent kinetic energy computation mode '",
        db["turbulent_kinetic_energy_mode"], "'.");
    }
  }

  if (db["compute_effective_viscosity"])
  {
    Output::root_info("TNSE", "Computing effective viscosity...");

    if (nu_eff_funct == nullptr)
    {
      Output::root_warn("TNSE", "Effective viscosity computation not set up.");
    }
    else
    {
      ComputeEffectiveViscosity(get_velocity(), *nu_eff_funct,
        ViscositySettings(db));
    }
  }

#endif

  example.do_post_processing(*this, zero_vorticity);

  outputWriter.write(t);

  checkpoint_io.write(s.solution, time_stepping_scheme);

  if (db["output_compute_time_average"])
  {
    time_average_io.write(s.time_avg_sol, time_stepping_scheme);
  }
}

/* ************************************************************************** */
template<int d>
std::array<double, int(10)> TimeNavierStokes<d>::get_errors() const
{
  return this->errors;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::output_problem_size_info() const
{
  int i_am_root = true;
#ifndef _MPI
  auto & velocity_space = *this->systems.front().velocity_space;
  auto & pressure_space = *this->systems.front().pressure_space;

  size_t nDofu = d*velocity_space.get_n_dof();
  size_t nDofp = pressure_space.get_n_dof();
  size_t nTotal = nDofu + nDofp;
  size_t nActive = d*velocity_space.get_n_active();

  auto coll = velocity_space.GetCollection();
  double hmin, hmax;
  coll->GetHminHmax(&hmin, &hmax);
  int n_cells = coll->GetN_Cells();
#else
  int my_rank;
  int root = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  i_am_root = (my_rank == root);
  auto velocity_comm = systems.front().velocity_space->get_communicator();
  auto pressure_comm = systems.front().pressure_space->get_communicator();
  int nDofu  = d*velocity_comm.get_n_global_dof();
  int nDofp  = pressure_comm.get_n_global_dof();
  int nTotal = nDofu + nDofp;
  
  auto coll = systems.front().velocity_space->GetCollection();
  int n_local_master_cells = coll->GetN_OwnCells();
  int n_cells;
  MPI_Reduce(&n_local_master_cells, &n_cells, 1, MPI_DOUBLE, MPI_SUM, root,
             MPI_COMM_WORLD);
  double local_hmin, local_hmax;
  coll->GetHminHmax(&local_hmin, &local_hmax);
  double hmin, hmax;
  MPI_Reduce(&local_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);
  MPI_Reduce(&local_hmax, &hmax, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
#endif

  if(i_am_root)
  {
    Output::print("N_Cells     : ", setw(10), n_cells);
    Output::print("h (min,max) : ", setw(12), hmin ," ", setw(12), hmax);
    Output::print("dof Velocity: ", setw(10), nDofu);
    Output::print("dof Pressure: ", setw(10), nDofp);
    Output::print("dof all     : ", setw(10), nTotal);
#ifndef _MPI
    Output::print("active dof  : ", setw(10), nActive);
#endif
  }
}


/* ************************************************************************** */
template <int d>
bool TimeNavierStokes<d>::linearized()
{
  if(db["time_discretization_nonlinear_term"].is("imex") 
    || db["time_discretization_nonlinear_term"].is("fully_explicit"))
  {
    return true;
  }
  return false;
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::modify_slip_bc(bool BT_Mass, bool slip_A_nl)
{
  // modification of the matrices due to the
  // slip type boundary conditions: If the mass matrices,
  // the off-diagonal A-blocks , and the BT's block,
  // are unchanged during the time iteration, then this modification
  // is done only once in the time loop. However, in the SUPG
  // and residual based VMS method these matrices are also
  // updated during the time steps, so modification of all
  // of them including the right-hand side is necessary.The
  // modification of the diagonal A-blocks are necessary
  // in any case.
  if(TDatabase::ParamDB->NSTYPE < 4)
  {
    ErrThrow("Slip with friction b.c. is only implemented for NSTYPE 4 and 14");
  }
  using SquareMatrixD = typename Template_names<d>::SquareMatrixD;
  using MatrixD = typename Template_names<d>::MatrixD;
  using BoundaryValuesFunction
    = typename Template_names<d>::BoundaryValuesFunction;
  using BoundaryConditionFunction 
    = typename Template_names<d>::BoundaryConditionFunction;
  
  std::vector<const FESpace*> spaces_mat(1);
  std::vector<double*> rhs_array(d);
  std::vector<const FESpace*> rhs_space(d);

  for(System_per_grid& s: this->systems)
  {
    spaces_mat[0] = s.velocity_space.get();
    for(int i = 0; i < d; ++i)
      rhs_space[i] = spaces_mat[0];

    for(int i = 0; i < d; ++i)
      rhs_array[i] = s.rhs.block(i);

    auto blocks = s.matrix.get_blocks_uniquely();

    auto mass_blocks = s.mass_matrix.get_blocks_uniquely(true);
    
    
    std::array<BoundaryConditionFunction*, d+1> boundCondition;
    std::array<BoundaryValuesFunction*, d+1> boundValues;
    for(int i = 0; i < d; ++i)
      boundCondition[i] = s.velocity_space->get_boundary_condition();
    boundCondition[d] = s.pressure_space->get_boundary_condition();
    for(int i = 0; i < d+1; ++i)
      boundValues[i] = example.get_bd(i);

    std::vector<SquareMatrixD*> sqMat;
    std::vector<MatrixD*> reMat;
    sqMat.resize(d);
    
   /* all d*d A blocks at the first time step
    * ------------------
    * a11 a12 a12 b1t
    * a21 a22 a23 b2t
    * a31 a32 a33 b3t
    * b1  b2  b3  c
    * and only the first 2 within the nonlinear loop */
    for(int i = 0; i < d; ++i)
      sqMat[i] = reinterpret_cast<SquareMatrixD*>(blocks[i*(d+2)].get()); // aii

    // if the off-diagonal are not changing within the non-linear loop
    // then dont need to assemble them again
    if(slip_A_nl)
    {
      sqMat.resize(d*d);
#ifdef __2D__
      sqMat[2] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(1).get());//a12
      sqMat[3] = reinterpret_cast<TSquareMatrix2D*>(blocks.at(3).get());//a21
#else // __3D__
      sqMat[3] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(1).get());//a12
      sqMat[4] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(2).get());//a13
      sqMat[5] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(4).get());//a21
      sqMat[6] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(6).get());//a23
      sqMat[7] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(8).get());//a31
      sqMat[8] = reinterpret_cast<TSquareMatrix3D*>(blocks.at(9).get());//a32
#endif   
    }

    // either at the first time step
    // or every time step if M and B's are changing
    reMat.resize(0);
    if(BT_Mass)
    {
#ifdef __2D__
      sqMat.resize(8);
      sqMat[4] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(0).get());//m11
      sqMat[5] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(4).get());//m22
      sqMat[6] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(1).get());//m12
      sqMat[7] = reinterpret_cast<TSquareMatrix2D*>(mass_blocks.at(3).get());//m21
      reMat.resize(2);
      reMat[0] = reinterpret_cast<TMatrix2D*>(blocks.at(2).get()); //the standing B blocks
      reMat[1] = reinterpret_cast<TMatrix2D*>(blocks.at(5).get());
#else // __3D__
      sqMat.resize(18);
      sqMat[9]  = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(0).get()); //m11
      sqMat[10] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(5).get()); //m22
      sqMat[11] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(10).get());//m33
  
      sqMat[12] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(1).get());//m12
      sqMat[13] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(2).get());//m13
      sqMat[14] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(4).get());//m21
      sqMat[15] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(6).get());//m23
      sqMat[16] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(8).get());//m31
      sqMat[17] = reinterpret_cast<TSquareMatrix3D*>(mass_blocks.at(9).get());//m32
      reMat.resize(3);
      reMat[0] = reinterpret_cast<TMatrix3D*>(blocks.at(3).get()); // b1t
      reMat[1] = reinterpret_cast<TMatrix3D*>(blocks.at(7).get()); // b2t
      reMat[2] = reinterpret_cast<TMatrix3D*>(blocks.at(11).get());// b3t
#endif  
    }

    auto u1 = s.u.GetComponent(0);
    auto u2 = s.u.GetComponent(1);
    // update the matrices and right hand side
#ifdef __2D__
    Assemble2DSlipBC(
#else // __3D__
    Assemble3DSlipBC(
#endif      
      spaces_mat.size(), spaces_mat.data(), sqMat.size(), sqMat.data(), 
      reMat.size(), reMat.data(), rhs_array.size(), rhs_array.data(), 
      rhs_space.data(), boundCondition.data(), boundValues.data()
#ifdef __2D__
      , u1.get(), u2.get());
#else
    );
#endif
  }
}

/* ************************************************************************* */
template <int d>
std::unique_ptr<typename TimeNavierStokes<d>::FEFunction>
  TimeNavierStokes<d>::get_velocity_component(int i)
{
  if(i >= 0 && i < d)
    return this->systems.front().u.GetComponent(i);
  else
    ErrThrow("There are only ", d, " velocity components!");
}


/* ************************************************************************* */
template <int d>
const Residuals& TimeNavierStokes<d>::get_residuals() const
{
  return old_residuals.back();
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_impuls_residual() const
{
  return old_residuals.back().momentumResidual;
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_mass_residual() const
{
  return old_residuals.back().massResidual;
}

/* ************************************************************************* */
template <int d>
double TimeNavierStokes<d>::get_full_residual() const
{
  return old_residuals.back().fullResidual;
}

/* ************************************************************************* */
template <int d>
void TimeNavierStokes<d>::reset_residuals()
{
  this->old_residuals = FixedSizeQueue<10, Residuals>();
}

/* ************************************************************************* */
template <int d>
void TimeNavierStokes<d>::update_matrices_lps(System_per_grid &s)
{
#ifdef __2D__
  std::vector<std::shared_ptr<FEMatrix>> blocks;
  blocks = s.matrix.get_blocks_uniquely();
  if(TDatabase::ParamDB->NSTYPE == 3 || TDatabase::ParamDB->NSTYPE == 4)
  {
    //update matrices for local projection stabilization
    std::vector<SquareMatrixD*> sqMat(2);
    sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    sqMat[1] = reinterpret_cast<SquareMatrixD*>(blocks.at(4).get());
    UltraLocalProjection(sqMat[0], false);
    UltraLocalProjection(sqMat[1], false);
  }
  else
  {
    std::vector<SquareMatrixD*> sqMat(1);
    sqMat[0] = reinterpret_cast<SquareMatrixD*>(blocks.at(0).get());
    UltraLocalProjection(sqMat[0], false);
  }
#else // 2D->3D
  (void)s; // suppress compiler warning
  ErrThrow("TimeNavierStokes<d>::update_matrices_lps is only implemented in "
           "2D.");
#endif
}

/* ************************************************************************** */
template <int d>
void TimeNavierStokes<d>::time_averaging()
{

  double t           = time_stepping_scheme.current_time_;
  double tau         = TDatabase::TimeDB->TIMESTEPLENGTH;
  double t0          = db["time_start"];
  double t0_avg      = db["start_time_averaging_at"];
  System_per_grid& s = this->systems.front();

  if( t == t0 )
  {
    return; // in case of restart (i.e. continue_output_after_restart)
  }
  else if( t-tau >= t0_avg )
  {
    s.time_avg_sol.scale(t - t0_avg - tau);
    s.time_avg_sol.add_scaled(s.solution_m2, tau/2.);
    s.time_avg_sol.add_scaled(s.solution, tau/2.);
    s.time_avg_sol.scale(1./(t - t0_avg));
  }
}

/* ************************************************************************** */
template<int d> 
void TimeNavierStokes<d>::adjust_pressure()
{
  System_per_grid& s = this->systems.front();

  int sign = 0;

  if (db["nse_nonlinear_form"].is("rotational"))
  {
    sign = -1;
  }
  else if (db["nse_nonlinear_form"].is("emac"))
  {
    sign = 1;
  }

  if (sign)
  {
    std::array<std::unique_ptr<FEFunction>, d> velocity_components;
    for (int i = 0; i < d; ++i)
    {
      velocity_components[i] = s.u.GetComponent(i);
    }

    typename FEFunction::AnalyticFunction 
    f = [&velocity_components, sign](const TBaseCell* cell, int i,
                                     std::array<double, d> xyz)
        {
          double val = 0;

          for (int c = 0; c < d; ++c)
          {
            double val_ui;
#ifdef __3D__
            velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                   xyz[2], &val_ui);
#else
            velocity_components[c]->FindValueLocal(cell, i, xyz[0], xyz[1],
                                                   &val_ui);
#endif
            val += val_ui*val_ui;
          }

          return sign * 0.5 * val;
        };

    s.p.add(f);

    // project pressure if necessary
    if (s.matrix.pressure_projection_enabled())
    {
      s.p.project_into_L20();
    }
  }
}

/* ************************************************************************* */
template <int d>
void TimeNavierStokes<d>::assemble_boundary_terms(bool first_nonlinear_step)
{
  const ParameterDatabase e_db = example.get_database();
  int n_windkessel_bd = e_db.try_get_value("n_windkessel_bd", 0);
  int n_directional_bd = e_db.try_get_value("n_directional_do_nothing_bd", 0);
  int n_neumann_bd = e_db["n_neumann_bd"];
  int n_nitsche_bd = e_db["n_nitsche_bd"];

  /// @todo this part of the code needs still to be implemented dimension-independent
  for (System_per_grid& s: this->systems)
  {
    const FESpace* v_space = s.velocity_space.get();

    if (n_windkessel_bd)
    {
      // Windkessel BC
      std::vector<size_t> windkessel_id = e_db["windkessel_id"];

      double theta = e_db["windkessel_theta"];

      for (unsigned int k = 0; k < windkessel_id.size(); k++)
      {
        double q_out;

        if (first_nonlinear_step)
        {
          if (theta == 1.0)
          {
            q_out = s.extrapolate_u.compute_flux(windkessel_id[k]);
          }
          else
          {
            BlockVector extrapolate_s(s.solution);
            extrapolate_s.scale(1.0 - theta);
            extrapolate_s.add_scaled(s.extrapolate_sol, theta);
            extrapolate_s.copy_nonactive(s.solution);

            FEVectFunct extrapolate_u(s.velocity_space, "u",
                                      extrapolate_s.get_entries(), d);

            q_out = extrapolate_u.compute_flux(windkessel_id[k]);
          }
        }
        else
        {
          if (theta == 1.0)
          {
            q_out = s.u.compute_flux(windkessel_id[k]);
          }
          else
          {
            BlockVector interpolate_s(s.solution);
            interpolate_s.scale(theta);
            interpolate_s.add_scaled(s.solution_m1, 1.0 - theta);
            interpolate_s.copy_nonactive(s.solution);

            FEVectFunct interpolate_u(s.velocity_space, "u",
                                      interpolate_s.get_entries(), d);

            q_out = interpolate_u.compute_flux(windkessel_id[k]);
          }
        }

        // solve the model for pressure
        double p_out = wk_boundary[k].solve(q_out, TDatabase::TimeDB->TIMESTEPLENGTH);
        wk_p_out[k] = p_out;

        Output::root_info<2>("TNSE", "Windkessel BC on ", windkessel_id[k],
          ", Rp: ", wk_boundary[k].get_proximal_resistance(),
          ", Rd: ", wk_boundary[k].get_distal_resistance(),
          ", C: ", wk_boundary[k].get_capacitance(),
          ", Q: ", q_out, " p: ", p_out);

#ifdef __2D__
        BoundaryAssembling2D::rhs_g_v_n(s.rhs, v_space, nullptr,
          windkessel_id[k], -p_out);
#else
        std::vector<TBoundFace*> boundaryFaceList;

        v_space->GetCollection()
          ->get_face_list_on_component(windkessel_id[k], boundaryFaceList);

        BoundaryAssembling3D ba;
        ba.rhs_g_v_n(s.rhs, v_space, nullptr, boundaryFaceList,
          windkessel_id[k], -p_out);
#endif
      }
    } /* if(n_windkessel_bd) */

    if (n_neumann_bd)
    {
      // Neumann BC
      std::vector<size_t> neumann_id = e_db["neumann_id"];
      std::vector<double> neumann_value = e_db["neumann_value"];

      for (unsigned int k = 0; k < neumann_id.size(); k++)
      {
        Output::root_info<2>("TNSE", "Neumann BC on ", neumann_id[k],
          ", value = ", neumann_value[k]);

#ifdef __2D__
        BoundaryAssembling2D::rhs_g_v_n(s.rhs, v_space, nullptr,
          neumann_id[k], -neumann_value[k]);
#else
        std::vector<TBoundFace*> boundaryFaceList;

        v_space->GetCollection()
          ->get_face_list_on_component(neumann_id[k], boundaryFaceList);

        BoundaryAssembling3D ba;
        ba.rhs_g_v_n(s.rhs, v_space, nullptr, boundaryFaceList,
          neumann_id[k], -neumann_value[k]);
#endif
      }
    } /* if(n_neumann_bd) */

    if (n_directional_bd > 0)
    {
#if __2D__
      Output::root_warn("TNSE", "Directional Neumann/do-nothing BC are not "
        "supported in 2D!");
#elif __3D__
      if (!linearized() && !is_rhs_and_mass_matrix_nonlinear)
      {
        Output::root_warn("TNSE", "Directional Neumann/do-nothing BC will "
          "introduce explicitness to the simulation!");
      }

      double theta = e_db["directional_do_nothing_theta"];
      std::vector<size_t> directional_id = e_db["directional_do_nothing_id"];
      std::vector<double> directional_value = e_db["directional_do_nothing_value"];

      bool is_smooth = e_db["directional_do_nothing_type"].is("smooth");
      double delta = e_db["directional_do_nothing_delta"];
      double nu_D0 = e_db["directional_do_nothing_D0"];
      nu_D0 *= example.get_nu();

      n_directional_bd = std::min(n_directional_bd, (int)directional_id.size());

      if (n_directional_bd > 0)
      {
        BlockVector extrapolate_s;

        if (first_nonlinear_step)
        {
          extrapolate_s = s.solution_m1;

          if (theta != 0.0)
          {
            double tau_ratio = theta * TDatabase::TimeDB->TIMESTEPLENGTH / s.tau_m1;
            extrapolate_s.scale(1.0 + tau_ratio);
            extrapolate_s.add_scaled(s.solution_m2, -tau_ratio);
          }
        }
        else
        {
          if (theta != 1.0)
          {
            extrapolate_s = s.solution_m1;
            extrapolate_s.scale(1.0 - theta);
            extrapolate_s.add_scaled(s.solution, theta);
          }
          else
          {
            extrapolate_s = s.solution;
          }
        }

        extrapolate_s.copy_nonactive(s.solution);

#ifdef _MPI
        auto comms = s.matrix.get_communicators();
        for (size_t bl = 0; bl < comms.size(); ++bl)
        {
          comms[bl]->consistency_update(extrapolate_s.block(bl), 3);
        }
#endif

        BlockVector extrapolate_dsdt;
        FEVectFunct extrapolate_dudt;
        FEVectFunct extrapolate_u(s.velocity_space, "u",
                                  extrapolate_s.get_entries(), d);

        if (is_smooth && nu_D0 > 0.0)
        {
          if (first_nonlinear_step)
          {
            extrapolate_dsdt = s.solution_m1;
            extrapolate_dsdt -= s.solution_m2;
            extrapolate_dsdt *= 1.0 / s.tau_m1;
          }
          else
          {
            extrapolate_dsdt = s.solution;
            extrapolate_dsdt -= s.solution_m1;
            extrapolate_dsdt *= 1.0 / TDatabase::TimeDB->TIMESTEPLENGTH;
          }

          extrapolate_dudt = FEVectFunct(s.velocity_space, "du_dt",
                                         extrapolate_dsdt.get_entries(), d);
        }

        for (int k = 0; k < n_directional_bd; k++)
        {
          double value = 1.0;
          if (k < (int)directional_value.size())
          {
            value = directional_value[k];
          }

          Output::root_info<2>("TNSE", "Directional BC on ", directional_id[k],
            ", value = ", value);

          std::vector<TBoundFace*> boundaryFaceList;

          v_space->GetCollection()
            ->get_face_list_on_component(directional_id[k], boundaryFaceList);

          BoundaryAssembling3D ba;

          if (!is_smooth)
          {
            ba.rhs_directional_do_nothing(s.rhs, v_space, nullptr,
              extrapolate_u, boundaryFaceList, directional_id[k], value);
          }
          else if (nu_D0 > 0.0)
          {
            ba.rhs_directional_do_nothing_smoothstep(s.rhs, v_space,
              extrapolate_u, extrapolate_dudt,
              boundaryFaceList, directional_id[k],
              nu_D0, value, delta);
          }
          else
          {
            ba.rhs_directional_do_nothing_smoothstep(s.rhs, v_space,
              extrapolate_u, boundaryFaceList, directional_id[k],
              value, delta);
          }
        }
      }
#endif
    }

    if (n_nitsche_bd)
    {
      // Nitsche penalty for weak essential BC
      std::vector<size_t> nitsche_id = e_db["nitsche_id"];
      std::vector<double> nitsche_penalty = e_db["nitsche_penalty"];
      double effective_viscosity = this->example.get_nu();

#ifdef __2D__
      for (unsigned int k = 0; k < nitsche_id.size(); k++)
      {
        const FESpace* p_space = s.pressure_space.get();

        int sym_u = e_db["symmetric_nitsche_u"];
        int sym_p = e_db["symmetric_nitsche_p"];

        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];

        Output::root_info<2>(" Nitsche BC on ", nitsche_id[k],
          ", nitsche penalty: ", nitsche_penalty[k]);

        BoundaryAssembling2D::nitsche_bc(s.matrix, s.rhs,
                                         v_space, p_space,
                                         this->example.get_bd(0),
                                         this->example.get_bd(1),
                                         nitsche_id[k], nitsche_penalty[k],
                                         effective_viscosity,
                                         sigma, L_0,
                                         sym_u, sym_p);
      }

      double corner_stab = e_db["corner_stab"];
      if (corner_stab)
      {
        Output::print<2>(" Corner stabilization is applied, corner_stab = ",
                         corner_stab);
        double sigma = this->example.get_inverse_permeability();
        double L_0 = db["L_0"];
        corner_stab = corner_stab * (effective_viscosity + sigma * L_0 * L_0);

        BoundaryAssembling2D::matrix_and_rhs_corner_stabilization(
                                                       s.matrix, s.rhs,
                                                       v_space,
                                                       this->example.get_bd(0),
                                                       this->example.get_bd(1),
                                                       nitsche_id,
                                                       corner_stab);
      } /* if(corner_stab) */
#else
      const FESpace* p_space = s.pressure_space.get();

      std::vector<TBoundFace*> boundaryFaceList;

      int sym_u = e_db["symmetric_nitsche_u"];
      int sym_p = e_db["symmetric_nitsche_p"];

      for (size_t k = 0; k < nitsche_id.size(); k++)
      {
        Output::root_info<2>("TNSE", "Nitsche BC on ", nitsche_id[k],
          ", nitsche penalty: ", nitsche_penalty[k]);

        boundaryFaceList.clear();
        v_space->GetCollection()
          ->get_face_list_on_component(nitsche_id[k], boundaryFaceList);

        BoundaryAssembling3D ba;
        ba.nitsche_bc(s.matrix, s.rhs,
                      v_space, p_space,
                      nullptr, nullptr, nullptr, // from Example
                      nullptr, // given FE function
                      boundaryFaceList,
                      nitsche_id[k], nitsche_penalty[k],
                      effective_viscosity,
                      sym_u, sym_p);
      }
#endif
    } /* if(n_nitsche_bd) */
  } /* for(System_per_grid& s : this->systems) */
}

template<int d>
void TimeNavierStokes<d>::prepare_vms(const TDomain& domain)
{
  auto collections = domain.get_grid_collections();

  // projection-based VMS, only on the finest level

  // VMS projection order

  int projection_order = db["vms_projection_space_order"];

  if (projection_order < 0)
  {
    projection_order = 0;
  }

  // projection space
  projection_space_ =
     std::make_shared<FESpace>(collections.front(), "L", example.get_bc(2),
                               DiscP_PSpace, projection_order);

  int ndofP = projection_space_->get_n_dof();

  // create vector for vms projection, used if small resolved scales are needed
  this->vms_small_resolved_scales.resize(6 * ndofP);

  // finite element vector function for vms projection
  this->vms_small_resolved_scales_fefct =
     std::make_shared<FEVectFunct>(projection_space_, "v",
                                   &vms_small_resolved_scales[0], 6);

  // create the label space, used in the adaptive method
  int n_cells = collections.front()->GetN_Cells();

  // initialize the piecewise constant vector
  if (projection_order == 0)
  {
    this->label_for_local_projection.resize(n_cells, 0);
  }
  else if (projection_order == 1)
  {
    this->label_for_local_projection.resize(n_cells, 1);
  }
  else
  {
    ErrThrow("local projection space is not defined");
  }

  // create fefunction for the labels such that they can be passed to the assembling routines
  label_for_local_projection_space_ =
    std::make_shared<FESpace>(collections.front(), "label_for_local_projection",
                              example.get_bc(2), DiscP_PSpace, 0);

  // finite element function for local projection space
  this->label_for_local_projection_fefct =
    std::make_shared<FEFunction>(label_for_local_projection_space_,
                                 "vms_local_projection_space_fefct",
                                 &label_for_local_projection[0]);

  // matrices for the vms method needs to assemble only on the
  // finest grid. On the coarsest grids, the SMAGORINSKY model
  // is used and for that matrices are available on all grids:

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1: case 2:
      ErrThrow("VMS projection cannot be supported for NSTYPE  ",
               TDatabase::ParamDB->NSTYPE);
      break;

    case 3: case 4: case 14:
      auto velocity_space = this->get_velocity_space();

      // matrices G_tilde
      matrices_for_turb_mod.push_back( std::make_shared<FEMatrix>(velocity_space, projection_space_));
      matrices_for_turb_mod.push_back( std::make_shared<FEMatrix>(velocity_space, projection_space_));
      matrices_for_turb_mod.push_back( std::make_shared<FEMatrix>(velocity_space, projection_space_));

      // matrices G
      matrices_for_turb_mod.push_back(std::make_shared<FEMatrix>(projection_space_, velocity_space));
      matrices_for_turb_mod.push_back(std::make_shared<FEMatrix>(projection_space_, velocity_space));
      matrices_for_turb_mod.push_back(std::make_shared<FEMatrix>(projection_space_, velocity_space));

      // mass matrix 
      matrices_for_turb_mod.push_back(std::make_shared<FEMatrix>(projection_space_, projection_space_));
      break;
  }
}

template<int d>
void TimeNavierStokes<d>::get_time_derivative(BlockVector& dt)
{
  System_per_grid& s = systems.front();

  BlockFEMatrix tmp_matrix = s.matrix;
  BlockVector tmp_rhs = s.rhs;

  call_assembling_routine(s, LocalAssembling_type::TimeNavierStokesAll,
    false, false, true);

  s.matrix.apply_transpose(s.solution, dt);
  dt -= s.rhs;

  s.matrix = tmp_matrix;
  s.rhs = tmp_rhs;

  dt.ResetBoundary();
}

#ifdef __3D__
template class TimeNavierStokes<3>;
#else
template class TimeNavierStokes<2>;
#endif
