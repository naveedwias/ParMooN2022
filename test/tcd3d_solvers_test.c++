/**
 * @brief A test program for the solution of Time dependent convection diffusion reaction equations
 * 
 * This test program is intended to check whether different solvers for CD3D
 * are able to solve three artificial CDR problems with analytic solution,
 * where the solution lies in the ansatz space.
 * It can be easily adapted to include more solvers, the control which solver
 * to use should come from the outside.
 * 
 * @author 
 * @history 14.05.2016
 * 
 */

#include <cmath>

#include "TimeConvectionDiffusion.h"
#include <Database.h>
#include <Domain.h>
#include <TimeDiscRout.h>
#include <TimeDiscretizations.h>
#include <Multigrid.h>
#include "ParMooN.h"

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
#endif

// function to compare the errors
void compare(const TimeConvectionDiffusion<3>& tcd3d,
             std::array<double, int(3)>errors, double tol)
{
  std::array<double, int(3)> computed_errors;
  computed_errors = tcd3d.get_errors();
  //check the L2-error
  if(std::abs(computed_errors.at(0)-errors.at(0)) > tol )
  {
    ErrThrow("L2(0,t,L2): " , computed_errors.at(0), "  ", errors.at(0));
  }
  if(std::abs(computed_errors.at(1) - errors.at(1)) > tol )
  {
    ErrThrow("L2(0,t,H1): " , computed_errors.at(1), "  ", errors.at(1));
  }
}

void check(ParameterDatabase& db, int ansatz_order, int time_disc,
           std::array<double, int(3)> errors, double tol)
{
  TDatabase::ParamDB->ANSATZ_ORDER = ansatz_order;
  if(time_disc == 1)
    db["time_discretization"] = "backward_euler";
  if(time_disc == 2)
    db["time_discretization"] = "crank_nicolson";
  TDatabase::TimeDB->TIMESTEPLENGTH = 0.1;
  db["time_step_length"] = 0.1;
  db["time_end"] = 1.;
  db["time_start"] = 0;
  TDatabase::TimeDB->CURRENTTIME = db["time_start"];
  
  TDomain domain(db);

  // Intial refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  
  // set the time discretization 
  SetTimeDiscParameters(0);

  // example object
  Example_TimeCD3D example_obj(db);
  TimeConvectionDiffusion<3> tcd3d(domain, db, example_obj);
  TimeDiscretization& tss = tcd3d.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.current_time_ = db["time_start"];
  // assemble the matrices and right hand side at the start time
  tcd3d.assemble_initial_time();
  
  while(tss.current_time_ < tss.get_end_time()-1e-10)
  {
    tss.current_step_++;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);

    double tau = tss.get_step_length();
    tss.current_time_ += tau;
    TDatabase::TimeDB->CURRENTTIME = tss.current_time_;
    
    Output::print<1>("\nCURRENT TIME: ", tss.current_time_);
    tcd3d.assemble();
    tcd3d.solve();    
    tcd3d.output();
    compare(tcd3d, errors, tol);
  }
}

// Choose the solver according to the input string and set global database
// entries accordingly.
void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  if(solver_name.compare("jacobi")==0)
  {
    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "richardson";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1.0e-13;
    db["residual_reduction"] =  0.0;
    db["max_n_iterations"] =  1000;
    db["min_n_iterations"] =  5;
    // Output::setVerbosity(2);
  }
  else if(solver_name.compare("multigrid") ==0)
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

    Output::setVerbosity(2);
  }
  else if(solver_name.compare("petsc") ==0)
  {
    db["solver_type"] = "petsc";
    db["max_n_iterations"] = 1000;
    db["residual_tolerance"] = 1.0e-12;
    db["residual_reduction"] =  0.0;
  }
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
  }
#endif
#ifdef _MPI
  else if(solver_name.compare("mumps") == 0)
  {
    db["solver_type"] = "direct";
  }
#endif
  else
  {
    throw std::runtime_error("Unknown solver for Time_CD3D problem!");
  }
}

double get_tolerance(std::string solver_name)
{//solver dependent tolerance?
  if (solver_name.compare("jacobi") == 0)
    return 1e-10;

  else if (solver_name.compare("multigrid") == 0)
    return 1e-10;
  else if(solver_name.compare("petsc") == 0)
    return 1e-10;
#ifndef _MPI
  else if(solver_name.compare("umfpack") == 0)
    return 1e-13 ;
#else
  else if (solver_name.compare("mumps") == 0)
    return 1e-13;
#endif
  else
    throw std::runtime_error("Unknown solver for Time_CD3D problem!");
}

int main(int /*argc*/, char* argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  int my_rank = 0;
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

  db.merge(Solver<>::default_solver_database(), true);
  db.merge(Multigrid::default_multigrid_database());
  db.merge(Example3D::default_example_database());
  db.merge(LocalAssembling3D::default_local_assembling_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  db["problem_type"] = 1;
  db.add("refinement_n_initial_steps",(size_t) 2,"",(size_t) 0, (size_t) 2);
  db["multigrid_n_levels"] = 1;
  // db.add("n_multigrid_levels", (size_t) 0, "",(size_t) 0, (size_t) 2);
  // db.add("solver_type", std::string("direct"), "", {"direct", "iterative"});
  // db.add("preconditioner", std::string("multigrid"), "",{"jacobi", "multigrid"});
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Hexa", "", 
         {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  // db["output_write_vtk"] = false;

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
  std::array<double, int(3)> errors;
  errors = {{0.0, 0.0, 0.0}};
  
  if(my_rank==0)
    Output::print<1>("Hexahedra grid.");
  db["geo_file"] = "Default_UnitCube_Hexa";
    
  db["example"] = -2; // Example -4: linear space and time 
  check(db, 1, 1, errors, tol); // time discretization is also included in the function
  
  db["example"] = -1; // Example -4: quadratic space time
  check(db, 2, 2, errors, tol); // time discretization is also included in the function
  
    if(my_rank==0)
      Output::print<1>("Tetrahedra grid.");
    db["geo_file"] = "Default_UnitCube_Tetra";
    
    db["example"] = -2; // Example -4: linear space and time 
    check(db, 1, 1, errors, tol); // time discretization is also included in the function
    
    db["example"] = -1; // Example -4: quadratic space time
    check(db, 2, 2, errors, tol); // time discretization is also included in the function
  
  parmoon::parmoon_finalize();
}
