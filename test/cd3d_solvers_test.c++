/**
 * @brief A test program for the solving of CD3D problems.
 *
 * This test program is intended to check whether different solvers for CD3D
 * are able to solve three artificial CDR problems with analytic solution,
 * where the solution lies in the ansatz space.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 *
 * @author Clemens Bartsch (heavily inspired by Najib's NSE2D Test program)
 *
 * @date 2016/03/31
 */

#include <cmath>

#include "ConvectionDiffusion.h"
#include "LocalAssembling.h"
#include <Solver.h>
#include <Multigrid.h>
#include "ParMooN.h"

#include <Database.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
#endif

void compare(const ConvectionDiffusion<3>& cd3d, std::array<double, int(2)> errors, double tol)
{
  // check the L2-error
  if(std::abs(cd3d.get_L2_error() - errors.at(0)) > tol)
  {
    ErrThrow("L2 norm: ", cd3d.get_L2_error(), "  ", errors.at(0));
  }
  // check the H1-error
  if(std::abs(cd3d.get_H1_semi_error() - errors.at(1)) > tol)
  {
    ErrThrow("H1-semi norm: ", cd3d.get_H1_semi_error(), "  ", errors.at(1));
  }
}

void check(ParameterDatabase& db, int ansatz_order,
           std::array<double, int(2)> errors, double tol)
{
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
  // fresh domain object
  TDomain domain(db);

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  // Choose and construct example.
  Example_CD3D example_obj(db);

  // Construct the cd3d problem object.
  ConvectionDiffusion<3> cd3d(domain, db, example_obj);

  cd3d.assemble();
  cd3d.solve();
  cd3d.output();

  compare(cd3d, errors, tol);
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{

  if (solver_name.compare("jacobi") == 0)
  {
    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "richardson";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1.0e-10;
    db["residual_reduction"] =  0.0;
    db["max_n_iterations"] =  1000;
    db["min_n_iterations"] =  5;
  }
  else if (solver_name.compare("multigrid") == 0)
  {
    db.merge(Multigrid::default_multigrid_database() ,true);

    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "fgmres";
    db["preconditioner"] = "multigrid";
    db["refinement_n_initial_steps"] = 2;
    db["multigrid_n_levels"] = 2;
    db["max_n_iterations"] =  100;
    db["residual_tolerance"] = 1.0e-15;
    db["residual_reduction"] =  0.0;
    // Multigrid parameters
    db["multigrid_cycle_type"] = "W";
    db["multigrid_smoother"] = "sor";
    db["multigrid_smoother_coarse"] = "direct_solve";
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 3;
    db["multigrid_n_post_smooth"] = 3;

    Output::setVerbosity(3);

  }
  else if(solver_name.compare("petsc") == 0)
  {
    db["solver_type"] = "petsc";
    db["max_n_iterations"] = 1000;
    db["residual_tolerance"] = 1.0e-13;
    db["residual_reduction"] =  0.0;
  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
  }
#endif
#ifdef _MPI
  else if (solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
  }
#endif
  else
  {
    ErrThrow("Unknown solver for CD3D problem! ", solver_name);
  }

}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?
  if (solver_name.compare("jacobi") == 0)
    return 1e-8;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-12;
  else if(solver_name.compare("petsc") == 0)
    return 1e-12;
  
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#else
  else if (solver_name.compare("mumps") == 0)
    return 1e-13;
#endif

  else
    ErrThrow("Unknown solver for CD3D problem!");
  return 0;
}


int main(int /*argc*/, char* argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  int my_rank = 0;
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  db.merge(Solver<>::default_solver_database() ,true);
  db.merge(Example3D::default_example_database());
  db.merge(LocalAssembling3D::default_local_assembling_database());
  db["problem_type"] = 1;
  db.add("refinement_n_initial_steps",(size_t) 2,"",(size_t) 0, (size_t) 3);
  db.add("multigrid_n_levels", (size_t) 0, "",(size_t) 0, (size_t) 3);
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "",
	 {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  
  db["space_discretization_type"] = "galerkin"; //Galerkin discretization, nothing else implemented

  TDatabase::ParamDB->Par_P0 = 0; // process responsible for the output
  TDatabase::ParamDB->Par_P3 = 1; // use mesh partitioning with halo cells

  set_solver_globals(std::string(argv[1]), db);

  double tol = get_tolerance(std::string(argv[1]));

  //===========================================================
  if(my_rank==0)
    Output::print<1>("Starting computations with solver: ",
                     std::string(argv[1]), ".");
  //===========================================================
  std::array<double, int(2)> errors;
  errors = {{0.0, 0.0}};
  
  if(my_rank==0)
    Output::print<1>("Hexahedra grid.");
  db["geo_file"] = "Default_UnitCube_Hexa";
  
  db["example"] = -1; // Example -1: constant solution, order 1 elements
  check(db, 1, errors, tol);
  
  db["example"] = -2; // Example -2: linear solution, order 1 elements
  check(db, 1, errors, tol);
  
  db["example"] = -3; // Example -3: quadratic solution, order 2 elements
  check(db, 2, errors, tol);

    if(my_rank==0)
      Output::print<1>("Tetrahedra grid.");
    db["geo_file"] = "Default_UnitCube_Tetra";
    
    db["example"] = -1; // Example -1: constant solution, order 1 elements
    check(db, 1, errors, tol);
    
    db["example"] = -2; // Example -2: linear solution, order 1 elements
    check(db, 1, errors, tol);
    
    db["example"] = -3; // Example -3: quadratic solution, order 2 elements
    check(db, 2, errors, tol);

  parmoon::parmoon_finalize();
}
