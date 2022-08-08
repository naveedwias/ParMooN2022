/**
 * @brief A test program for the solving of NSE3D problems.
 *
 * This test program is intended to check whether different solvers for NSE3D
 * are able to solve three artificial Navier--Stokes problems with known
 * analytic solution, where the solutions lie in the ansatz spaces.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 *
 * The program tests a selection of combinations of
 *  - the polynomial examples -1 to -3,
 *  - the stable finite element pairs P2/P1, P3/P2, Q2/Q1, Q2/P1_disc, Q3/Q2
 *  - the NSTYPES 1,2,3,4
 * on the default unit cube geometry.
 * We're only testing examples whose analytic solution is in the ansatz space,
 * thus we expect very small errors every time.
 * Fixed are discretization_type 1 (galerkin), LAPLACETYPE 0,
 * 'nse_nonlinear_form' 0.
 * Note that for these discretizations it is futile to choose any other NSTYPE
 * than 1 (see e.g. MooNMD documentation p. 35), but we vary them anyway, just
 * for the sake of testing the different types.
 *
 * So far the test is adapted to:
 *  - testing the umfpack solver when compiled SEQUENTIAL
 *  - testing lsc preconditioned fgmres SEQUENTIAL
 *  - testing multigrid preconditioned fgmres SEQUENTIAL
 *  - testing the mumps solver when compiled MPI
 *
 * The MPI Mumps test contains one example for the three combinations
 * 1 - 2 (hexa), 2 - 12 (hexa), 4 - 2 (hexa) of nstype and velocity space.
 * The tests for 3rd order elements cannot be run, because these elements are
 * not fitted for mpi yet.
 *
 * @author Clemens Bartsch (heavily inspired by Najib's NSE2D Test program)
 *
 * @date 2016/04/04
 */

#include <cmath>
#include "NavierStokes.h"

#include <Database.h>
#include <Multigrid.h>
#include <iomanip>
#include "LocalAssembling.h"
#include "Saddle_point_preconditioner.h"
#include "Utilities.h"
#include "ParMooN.h"

#ifdef _MPI
#include <mpi.h>
#endif

void compare(const NavierStokes<3>& nse3d, std::array<double, int(5)> errors,
             double tol)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  std::array<double, 8> computed_errors = nse3d.get_errors();

  // check the L2-error of the velcoity
  if( my_rank == 0 && (!utilities::are_equal(computed_errors[0], errors[0], tol) ||
      computed_errors[0] != computed_errors[0])) //check for nan!
  {
    ErrThrow("L2 norm of velocity: ", std::setprecision(14), computed_errors[0], "  ", errors[0], "  ", std::abs(computed_errors[0] - errors[0]));
  }
  // check the H1-error of the velcoity
  if( my_rank == 0 && (!utilities::are_equal(computed_errors[1], errors[1], tol) ))
  {
    ErrThrow("H1 norm of velocity: ", std::setprecision(14), computed_errors[1], "  ", errors[1], "  ", std::abs(computed_errors[1] - errors[1]));
  }
  // check the L2-error of the pressure
  if( my_rank == 0 && (!utilities::are_equal(computed_errors[3], errors[3], tol)))
  {
    ErrThrow("L2 norm of pressure: ", std::setprecision(14), computed_errors[3], "  ", errors[3], "  ", std::abs(computed_errors[3] - errors[3]));
  }
  // check the H1-error of the pressure
  if(my_rank == 0 && (!utilities::are_equal(computed_errors[4], errors[4], tol) ))
  {
    ErrThrow("H1 norm of pressure: ", std::setprecision(14), computed_errors[4], "  ", errors[4], "  ", std::abs(computed_errors[4] - errors[4]));
  }
}
void check(ParameterDatabase& db, const TDomain& domain, int velocity_order,
           int pressure_order, int nstype, std::array<double, int(5)> errors,
           double tol)
{
  int my_rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  if (my_rank ==0)
  {
    Output::print("******* Check ", db["example"], " ", velocity_order, " ",
                  pressure_order, " ", nstype, " *******");
  }

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;

  //Perform usual checks on the parameter consistency
  check_parameters_consistency_NSE(db);

  // Construct the nse3d problem object.
  NavierStokes<3> nse3d(domain, db);

  nse3d.assemble_linear_terms();

  // check stopping criterion
  nse3d.stop_it(0);
  for(unsigned int k=1;; k++)
  {
    Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
                     nse3d.get_residuals());
    nse3d.solve();

    // checking the first nonlinear iteration
    nse3d.assemble_nonlinear_term();;
    if(nse3d.stop_it(k))
      break;
  }
  nse3d.output();
  // now compare the errors
  compare(nse3d, errors, tol);
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  db["solver_type"] = std::string("iterative");
  db["direct_solver_type"] = std::string("umfpack");
  db["iterative_solver_type"] = std::string("fgmres");
  db["preconditioner"] = std::string("no_preconditioner");
  db["residual_tolerance"] = 1.0e-13;
  
  if (solver_name.compare("lsc") == 0)
  {
    db["preconditioner"] = "least_squares_commutator";
    db["nonlinloop_epsilon"] = 1e-12;
    // If you want to run tests with an iterative solver for the
    // velocity block, comment out these lines of code here.
    ParameterDatabase velocity_solver_db = Solver<>::default_solver_database();
    velocity_solver_db.set_name(Saddle_point_preconditioner::database_name_velocity_solver);
    velocity_solver_db["solver_type"] = "iterative";
    velocity_solver_db["iterative_solver_type"] = "bi_cgstab";
    velocity_solver_db["preconditioner"] = "ssor";
    velocity_solver_db["sor_omega"] = 1.0;
    velocity_solver_db["max_n_iterations"] = 1000;
    velocity_solver_db["residual_tolerance"] = 1.0e-10;
    velocity_solver_db["residual_reduction"] = 1.0e-10;
    velocity_solver_db["damping_factor"] = 1.0;
    db.add_nested_database(velocity_solver_db);
    ParameterDatabase pressure_solver_db = velocity_solver_db; // copy
    pressure_solver_db.set_name(Saddle_point_preconditioner::database_name_pressure_solver);
    // choose a different solver compared with the velocity solver
    pressure_solver_db["solver_type"] = "direct";
    db.add_nested_database(pressure_solver_db);
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    db.merge(Multigrid::default_multigrid_database());
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 1;
    //control nonlinear loop
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 5;
    // New multigrid parameters
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "patch_vanka_store";
    db["multigrid_smoother_coarse"] = "nodal_vanka_store";
    db["multigrid_correction_damp_factor"] = 1.0;
    db["multigrid_n_pre_smooth"] = 1;
    db["multigrid_n_post_smooth"] = 1;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["multigrid_vanka_damp_factor"]=1.0;

  }
  else if (solver_name.compare("cell_vanka_jacobi") == 0)
  {
    db.merge(Multigrid::default_multigrid_database());
    db["iterative_solver_type"] = std::string("fgmres");
    db["preconditioner"] = "multigrid";
    db["residual_tolerance"] = 1.0e-8;
    db["max_n_iterations"] = 10;
    db["min_n_iterations"] = 3;

    db["refinement_n_initial_steps"] = 2;

    //control nonlinear loop
    db["nonlinloop_epsilon"] = 1e-12;
    db["nonlinloop_maxit"] = 15;
    // New multigrid parameters
    db["multigrid_n_levels"] = 2;
    db["multigrid_cycle_type"] = "V";
    db["multigrid_smoother"] = "cell_vanka_jacobi";
    db["multigrid_smoother_coarse"] = "cell_vanka_jacobi";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 1;
    db["multigrid_n_post_smooth"] = 1;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 10;
    db["multigrid_vanka_damp_factor"]=0.8;


  }
  else if(solver_name.compare("petsc") == 0)
  {
    db["solver_type"] = "petsc";
  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "umfpack";
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 6;
  }
  else if(solver_name.compare("pardiso") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "pardiso";
    ErrThrow("pardiso not yet set!");
  }
#else
  else if (solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "mumps";
    db["nonlinloop_epsilon"] = 1e-15;
    db["nonlinloop_maxit"] = 6;
  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for NSE3D problem!");
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?

  if(solver_name.compare("cell_vanka_jacobi") == 0)
    return 1e-8;
  if(solver_name.compare("petsc") == 0)
    return 1e-9;
  if(solver_name.compare("lsc") == 0)
    return 1e-9;
#ifndef _MPI
  if(solver_name.compare("umfpack") == 0)
    return 1e-9;
  if(solver_name.compare("multigrid") == 0)
    return 1e-8;
#else
  if(solver_name.compare("mumps") == 0)
    return 1e-9 ;
#endif
    throw std::runtime_error("Unknown solver for NSE3D problem!");

  return 0;
}


int main(int argc, char* argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  
  int my_rank = 0;
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  Output::setVerbosity(3);

  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(LocalAssembling3D::default_local_assembling_database());

  db["problem_type"].set<size_t>(5);
  
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 3);

  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5; // flow problem type

  db["space_discretization_type"] = "galerkin"; //Galerkin discretization, nothing else implemented
  db["nse_nonlinear_form"] = "convective";
  TDatabase::ParamDB->LAPLACETYPE = 0;

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  set_solver_globals(std::string(argv[1]), db);

  double tol = get_tolerance(std::string(argv[1]));

  //===========================================================
  if(my_rank==0)
    Output::print<1>(">>>>> Starting computations with solver: <<<<<"
        , std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(5)> errors;
  errors = {{0.0, 0.0, 0.0, 0.0, 0.0}};

  //============= Tests on hexa grid ==========================
  if(my_rank==0)
    Output::print<1>(">>>>> Hexahedra grid. <<<<<");
  //===========================================================
  {
    //do the domain thingy
    db.add("boundary_file", "Default_UnitCube", "");
    db.add("geo_file", "Default_UnitCube_Hexa", "",
	   {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
    TDomain domain_hex(db);

    // Intial refinement and grabbing of grids for multigrid.
    domain_hex.refine_and_get_hierarchy_of_collections(db);
    db.merge(Example_NSE3D::default_example_database());

    if(std::string(argv[1]) == std::string("cell_vanka_jacobi"))
    {//the cell vanka case - test only on Q2/P1-disc element (MPI parallelized & disc press)
      db["example"] = -1;
      size_t nstype = 4; //nstype should not matter much here
      check(db, domain_hex, 22, -4711, nstype, errors, tol);
#ifdef _MPI
      MPI_Finalize();
#endif
      return 0;
    }

    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q2/Q1 element on hexahedral grid. <<<<<");
      db["example"] = -3;
      size_t nstype = 4;
      Parameter p("petsc_arguments",
         std::string("-pc_type lu -pc_factor_mat_solver_package umfpack "
                     "-ksp_monitor"), "");
      db["petsc_arguments"].impose(p);
      check(db, domain_hex, 2, -4711, nstype, errors, tol);
    }
    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q2/P1^disc element on hexahedral grid. <<<<<");
      db["example"] = -3;
      size_t nstype = 2;
      check(db, domain_hex, 22, -4711, nstype, errors, tol);
    }
#ifndef _MPI//only for seq, 3rd order elements are not yet adapted for parallel
    {
      if(my_rank==0)
        Output::print<1>("\n>>>>> Q3/Q2 element on hexahedral grid. <<<<<");
      db["example"] = -4;
      size_t nstype = 3;
      check(db, domain_hex, 3, -4711, nstype, errors, tol);
    }
#endif
  }

  //============= Tests on tetra grid =========================
  if(my_rank==0)
    Output::print<1>(">>>>> Tetrahedral grid. <<<<<");
  //===========================================================
  {
    //do the domain thingy
    db["geo_file"]= "Default_UnitCube_Tetra";
    TDomain domain_tet(db);
    // Intial refinement and grabbing of grids for multigrid.
    domain_tet.refine_and_get_hierarchy_of_collections(db);

    {
      db["example"] = -3;
      size_t nstype = 4;
      if(my_rank==0)
        Output::print<1>("\n>>>>> P2/P1 element on tetrahedral grid. <<<<<");
      check(db, domain_tet, 2,-4711, nstype, errors, tol);
    }

    {
      db["example"] = -2;
      size_t nstype = 3;
      if(my_rank==0)
        Output::print<1>("\n>>>>> CR/P0 element on tetrahedral grid. <<<<<");
      check(db, domain_tet, -1,-4711, nstype, errors, tol);
    }

#ifndef _MPI
    {
      if(std::string(argv[1]) == std::string("multigrid"))
        return 0; //no convergence for multigrid
      if(my_rank==0)
        Output::print<1>("\n>>>>> P3/P2 element on tetrahedral grid. <<<<<");
      db["example"] = -4;
      size_t nstype = 4;
      check(db, domain_tet, 3,-4711, nstype, errors, tol);
    }
#endif
  }
  parmoon::parmoon_finalize();
}
