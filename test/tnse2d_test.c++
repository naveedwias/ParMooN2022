/**
 * @brief A test program to test Time-dependent Navier--Stokes
 * 
 * Residuals are checked at different time steps
 *
 * @history Reworked 22.10.2018 by Najib -
 * added solvers as arguments + activated lsc test
 */
#include <Domain.h>
#include <Database.h>
#include "TimeNavierStokes.h"
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <FEFunction2D.h>
#include <Multigrid.h> // newly implemented by clemens
#include "ParMooN.h"

void compare(TimeNavierStokes<2>& tnse2d, std::array<double, 4> errors,
             double tol)
{
  auto computed_errors = tnse2d.get_errors();

  // check the L2-error of the velocity
  if( std::abs(computed_errors[0]-errors[0]) > tol )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the H1-error of the velocity
  if( std::abs(computed_errors[2] - errors[1]) > tol )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[2], "  ", errors[1]);
  }
  // check the L2-error of the pressure
  if( std::abs(computed_errors[3] - errors[2]) > tol)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[3], "  ", errors[2]);
  }
  // check the H1-error of the pressure
  if(std::abs(computed_errors[4] - errors[3]) > tol )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[4], "  ", errors[3]);
  }
}

void compute(TDomain& domain, ParameterDatabase& db,
             std::array<std::array<double, 4>, 4> errors, double tol,
             int time_disc)
{
  TDatabase::TimeDB->TIME_DISC=time_disc;
  TDatabase::TimeDB->CURRENTTIME = db["time_start"];

  TimeNavierStokes<2> tnse2d(domain, db);

  tnse2d.get_time_stepping_scheme().current_step_ = 0;
  tnse2d.get_time_stepping_scheme().current_time_ = db["time_start"];
  tnse2d.get_time_stepping_scheme().set_time_disc_parameters();

  tnse2d.assemble_initial_time();
  tnse2d.output();
  double end_time = db["time_end"];
  while(TDatabase::TimeDB->CURRENTTIME < end_time-1e-10)
  {
    tnse2d.get_time_stepping_scheme().current_step_++;
    tnse2d.get_time_stepping_scheme().set_time_disc_parameters();
    double tau = tnse2d.get_time_stepping_scheme().get_step_length();
    TDatabase::TimeDB->CURRENTTIME += tau;    
    tnse2d.get_time_stepping_scheme().current_time_ += tau;
    
    Output::print("\nCURRENT TIME: ", TDatabase::TimeDB->CURRENTTIME);    
    // prepare the right hand side vector
    tnse2d.assemble_matrices_rhs(0);
    // nonlinear iteration
    for(unsigned int i=0;; i++)
    {
      if(tnse2d.stop_it(i))
        break;
      tnse2d.solve();
      tnse2d.assemble_matrices_rhs(i+1);
    }
    Output::print(" current step : ",
                  tnse2d.get_time_stepping_scheme().current_step_);
    // post processing: error computations
    tnse2d.output();
    std::cout << std::endl;
    // check the errors
    if(tnse2d.get_time_stepping_scheme().current_step_ == 1)
      compare(tnse2d, errors[0], tol);
    else if(tnse2d.get_time_stepping_scheme().current_step_ == 2)
      compare(tnse2d, errors[1], tol);
    else if(tnse2d.get_time_stepping_scheme().current_step_ == 20)
      compare(tnse2d, errors[2], tol);
  }
}

void check(TDomain& domain, ParameterDatabase& db, int velocity_order, int,
           int nstype, int laplace_type, int nonlinear_form, int time_disc,
           std::array<std::array<double, int(4)>,4> errors, double tol)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  if(laplace_type == 0)
    db["laplace_type_deformation"] = false;
  else if(laplace_type == 1)
    db["laplace_type_deformation"] = true;
  else
    ErrThrow("laplace_type can only be 0 or 1. ", laplace_type);
  if(nonlinear_form == 0)
    db["nse_nonlinear_form"] = "convective";
  else if(nonlinear_form == 1)
    db["nse_nonlinear_form"] = "skew_symmetric";
  else if(nonlinear_form == 2)
    db["nse_nonlinear_form"] = "rotational";
  else if(nonlinear_form == 3)
    db["nse_nonlinear_form"] = "emac";
  else if(nonlinear_form == 4)
    db["nse_nonlinear_form"] = "divergence";
  else
    ErrThrow("other nonlinear forms are not yet tested");

  // rough check
  if(time_disc == 1)
    db["time_discretization"] = "backward_euler";
  else if(time_disc == 2)
    db["time_discretization"] = "crank_nicolson";
  else if(time_disc == 3){
    db["time_discretization"]= "bdf_two";
  }

  db["time_start"]= 0.;
  db["time_end"]= 1.;
  db["time_step_length"]= 0.05;

  Output::print("============Time Discretization=============== ", db["time_discretization"]);
  compute(domain,db,errors,tol,time_disc);
}

void set_solver_globals(std::string solver_name, ParameterDatabase& db)
{
  db["solver_type"] = std::string("iterative");
  db["direct_solver_type"] = std::string("umfpack");
  db["iterative_solver_type"] = std::string("fgmres");
  db["preconditioner"] = std::string("no_preconditioner");
  db["residual_tolerance"] = 1.0e-12;

  if (solver_name.compare("lsc") == 0)
  {
    db["preconditioner"] = "least_squares_commutator";
    db["nonlinloop_epsilon"] = 1e-10;
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
    db["multigrid_correction_damp_factor"] = 0.8;
    db["multigrid_n_pre_smooth"] = 2;
    db["multigrid_n_post_smooth"] = 2;
    db["multigrid_coarse_residual"] = 1.0e-1;
    db["multigrid_coarse_max_n_iterations"] = 5;
    db["multigrid_vanka_damp_factor"]=0.7;

  }
  else if(solver_name.compare("umfpack") == 0)
  {
    db["solver_type"] = "direct";
    db["direct_solver_type"] = "umfpack";
    db["nonlinloop_epsilon"] = 1e-10;
    db["nonlinloop_maxit"] = 5;
  }
  else
  {
    throw std::runtime_error("Unknown solver for Time_NSE2D problem!");
  }
}

double get_tolerance(std::string solver_name)
{
  double val;
  if(solver_name.compare("umfpack") == 0)
    val = 1e-7;
  if(solver_name.compare("lsc") == 0)
    val = 1e-7;
  if(solver_name.compare("multigrid") == 0)
    val = 1e-7;
  return val;
}

void set_errors(std::array<std::array<double, int(4)>,4>& errors)
{
  errors[0] = {{0.0, 0.0, 0.0, 0.0}};
  errors[1] = {{0.0, 0.0, 0.0, 0.0}};
  errors[2] = {{0.0, 0.0, 0.0, 0.0}};
}

int main(int, char* argv[])
{
  parmoon::parmoon_initialize();
  bool testall = false;
  if (argv[2])
    testall = (std::string(argv[2]).compare("testall") == 0);
  ParameterDatabase db = TimeNavierStokes<2>::default_tnse_database();
  db.merge(Solver<>::default_solver_database());
  db.merge(Example2D::default_example_database());
  db.merge(TimeDiscretization::default_TimeDiscretization_database());
  
  db["problem_type"].set<size_t>(6);
  db["nonlinloop_slowfactor"]=1.; 
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "", {"UnitSquare","TwoTriangles"});
  db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2);
  
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 6; // flow problem type
  db["space_discretization_type"] = "galerkin"; //Galerkin discretization, nothing else implemented
  
  set_solver_globals(std::string(argv[1]), db);
  double tol = get_tolerance(std::string(argv[1]));
  std::array<std::array<double, int(4)>,4> errors;    
  //=======================================================================
  //============= Test 1 : Quads ==========================================
  //=======================================================================
  {    
    TDomain domain_quads(db);
    domain_quads.refine_and_get_hierarchy_of_collections(db);
    db["example"] = 1;
    int lapacetype = 0; 
    int nonlineartype = 0;
    int time_disc = 1;
    Output::print<1>("Testing the Q2/P1-disc elements");
    set_errors(errors);
    check(domain_quads, db, 12, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    check(domain_quads, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      nonlineartype = 4; //divergence form
      check(domain_quads, db, 12, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      nonlineartype = 0; // convective form
      
      Output::print<1>("Testing the Q3/P2 elements");
      check(domain_quads, db, 13, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 1, TimeDiscretization 1=====================================    
    lapacetype = 1;
    Output::print<1>("Testing the Q2/P1-disc elements");
    set_errors(errors);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    if(testall)
    {
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      nonlineartype = 4; //divergence form    
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);    
      nonlineartype = 0; // convective form 
      
      Output::print<1>("Testing the Q3/P2-disc elements");
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 0, TimeDiscretization 2=====================================    
    time_disc = 2;
    lapacetype = 0;
    check(domain_quads, db, 12, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    check(domain_quads, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      nonlineartype = 4; //divergence form
      check(domain_quads, db, 12, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      nonlineartype = 0; // convective form
      
      Output::print<1>("Testing the Q3/P2-disc elements");
      check(domain_quads, db, 13, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 1, TimeDiscretization 2=====================================    
    lapacetype = 1; 
    Output::print<1>("Testing the Q2/P1-disc elements");
    set_errors(errors);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    if(testall)
    {
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      nonlineartype = 4; //divergence form
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);    
      nonlineartype = 0; // convective form
      
      Output::print<1>("Testing the Q3/P2-disc elements");
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    
    Output::print<1>("Testing imex scheme");
    // imex scheme 
    db["time_discretization_nonlinear_term"] = "imex";
    db["extrapolation_type"] = "linear";
    lapacetype = 0; 
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    lapacetype = 1; 
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing fully_explicit scheme");
    db["time_discretization_nonlinear_term"] = "fully_explicit";
    lapacetype = 0; 
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    lapacetype = 1; 
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    db["time_discretization_nonlinear_term"] = "fully_implicit";
    
    //============= lapacetype 0, TimeDiscretization 3=====================================    
    time_disc = 3;
    lapacetype = 0;
    Output::print<1>("Testing the Q2/P1-disc elements");
    check(domain_quads, db, 12, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    check(domain_quads, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_quads, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the Q3/P2-disc elements");
      check(domain_quads, db, 13, -4711, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 1, TimeDiscretization 2=====================================    
    lapacetype = 1; 
    Output::print<1>("Testing the Q2/P1-disc elements");
    set_errors(errors);
    if(testall)
    {
      check(domain_quads, db, 12, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    Output::print<1>("Testing the Q2/Q1 elements");
    if(testall)
    {
      check(domain_quads, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_quads, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the Q3/P2-disc elements");
      check(domain_quads, db, 13, -4711, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 13, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q3/Q2 elements");
      check(domain_quads, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the Q4/Q3 elements");
      check(domain_quads, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_quads, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
  }
  //=======================================================================
  //============= Test 1 : TwoTriangles====================================
  //=======================================================================
  {
    ParameterDatabase db = TimeNavierStokes<2>::default_tnse_database();
    db.merge(Solver<>::default_solver_database());
    db.merge(Example2D::default_example_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "TwoTriangles", "", {"UnitSquare","TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 2,"", (size_t) 0, (size_t) 2);    
    
    TDomain domain_tri(db);
    domain_tri.refine_and_get_hierarchy_of_collections(db);
    
    db["example"] = 1;
    int lapacetype = 0; int nonlineartype = 0; int time_disc = 1;
    
    set_errors(errors);    
    Output::print<1>("Testing the P2/P1 elements");
    check(domain_tri, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_tri, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      nonlineartype = 4;
      check(domain_tri, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      nonlineartype = 0;
      
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
      Output::print<1>("Testing the P4/P3-disc elements");
      check(domain_tri, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 1, TimeDiscretization 1=====================================    
    lapacetype = 1;     
    set_errors(errors);    
    Output::print<1>("Testing the P2/P1 elements");
    if(testall)
    {
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
      Output::print<1>("Testing the P4/P3 elements");
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    //============= lapacetype 0, TimeDiscretization 2=====================================    
    time_disc = 2;
    lapacetype = 0;    
    Output::print<1>("Testing the P2/P1 elements");
    check(domain_tri, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_tri, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
      Output::print<1>("Testing the P4/P3-disc elements");
      check(domain_tri, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    
    //============= lapacetype 1, TimeDiscretization 2=====================================    
    lapacetype = 1; 
    set_errors(errors);    
    Output::print<1>("Testing the P2/P1 elements");
    if(testall)
    {
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
      Output::print<1>("Testing the P4/P3 elements");
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    
    //============= lapacetype 0, TimeDiscretization 3=====================================    
    time_disc = 3;
    lapacetype = 0;    
    Output::print<1>("Testing the P2/P1 elements");
    check(domain_tri, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      check(domain_tri, db, 2, 1, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    if(testall)
    {
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the P4/P3 elements");
      check(domain_tri, db, 4, 3, 1, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 2, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    
    //============= lapacetype 1, TimeDiscretization 3=====================================    
    lapacetype = 1; 
    set_errors(errors);
    Output::print<1>("Testing the P2/P1 elements");
    if(testall)
    {
      check(domain_tri, db, 2, 1, 3, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    if(testall)
    {
      Output::print<1>("Testing the P3/P2 elements");
      check(domain_tri, db, 3, 2, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 3, 2, 4, lapacetype, nonlineartype, time_disc, errors, tol);
      
      Output::print<1>("Testing the P4/P3 elements");
      check(domain_tri, db, 4, 3, 3, lapacetype, nonlineartype, time_disc, errors, tol);
      check(domain_tri, db, 4, 3, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    }
    // just one test with IMEX scheme
    Output::print<1>("Testing P2/P1 elements with imex, bdf2");
    db["time_discretization_nonlinear_term"] = "imex";
    db["extrapolation_type"] = "linear";
    lapacetype = 0; 
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    db["time_discretization_nonlinear_term"] = "fully_explicit";
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    db["time_discretization_nonlinear_term"] = "fully_implicit";
  }
  //=======================================================================
  //============= Old Test : Quads ==========================================
  //=======================================================================
  // test for quads
  {
    ParameterDatabase db = TimeNavierStokes<2>::default_tnse_database();
    db.merge(Solver<>::default_solver_database());
    db.merge(Example2D::default_example_database());
    db.merge(TimeDiscretization::default_TimeDiscretization_database());
    
    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "UnitSquare", "", {"UnitSquare","TwoTriangles"});
    db.add("refinement_n_initial_steps", (size_t) 1,"", (size_t) 0, (size_t) 2); 
    db.merge(ParameterDatabase::default_output_database());
    db["example"] = 0;
    
    //  declaration of databases
    TDomain domain(db);    
    domain.refine_and_get_hierarchy_of_collections(db);

    //=============================================================================
    // CRANK-NICOLSON TIME STEPPING SCHEME
    int time_disc=2; int laplace_type; int nl_form=0;
    //=============================================================================
    Output::print<1>("Testing the Q2/P1-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001112092818, 0.01455279725, 0.02027240218, 0.09770369527}};
    errors[1] = {{0.002261228092, 0.02900169754, 0.02290308963, 0.1548609258}};
    errors[2] = {{0.019173655437058, 0.24443064635394, 0.1019439334998, 1.2138692181926}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    if(testall)
    {
      check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
      check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.001100290828, 0.01461603719, 0.02009574353, 0.09725227236}};
    errors[1] = {{0.002224014719, 0.02909117943, 0.0220984703, 0.1527635601}};
    errors[2] = {{0.018817584454803, 0.24503529699343, 0.087813819645541, 1.1853696236832}};
    laplace_type = 1;
    if(testall)
    {
      check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    //=============================================================================
    Output::print<1>("Testing the Q3/P2-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{9.921682115e-05, 0.001899852996, 0.01927905366, 0.07644264853}};
    errors[1] = {{0.0001999247874, 0.003792810533, 0.01926073977, 0.07440290132}};
    errors[2] = {{0.001688312120429, 0.031969723353629, 0.015089809171047, 0.24683813828767}};
    laplace_type=0;
    check(domain, db, 13, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    if(testall)
    {
      check(domain, db, 13, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
      check(domain, db, 13, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 13, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.00010012275994756, 0.0019037255002744, 0.019279007255799, 0.076445700453971}};
    errors[1] = {{0.00020164329006377, 0.0038049726341925, 0.019260939678099, 0.074427427020247}};
    errors[2] = {{0.0017012192164184, 0.032072166618082, 0.015157980380392, 0.24774350544933}};
    laplace_type=1; 
    if(testall)
    {
      check(domain, db, 13, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 13, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    
    //=============================================================================
    // BACKWARD-EULER TIME STEPPING SCHEME
    time_disc = 1;
    //=============================================================================
    //=============================================================================
    Output::print<1>("Testing the Q2/P1-disc elements");
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.0011233147920272, 0.014527988618419, 0.0081517719198929, 0.074300286739799}};
    errors[1] = {{0.0022578528432466, 0.029004021823278, 0.013818200240636,0.14595696769918}};
    errors[2] = {{0.019188837321532, 0.24443049572246, 0.10004389887316, 1.2161726946242}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    if(testall)
    {
      check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
      check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    Output::print<1>("LAPLACETYPE: ", 1, " NSTYPE's: ", 3, ", " , 4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.0011076940040252, 0.014582185030344, 0.0069226734930213, 0.072278525899217}};
    errors[1] = {{0.0022218289200176, 0.029097247028845, 0.011676670594301, 0.14224315633584}};
    errors[2] = {{0.018829987146181, 0.24501329407406, 0.086168813250833, 1.1875928535871}};
    laplace_type = 1;
    if(testall)
    {
      check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
    
    //=============================================================================
    // BDF2 STEPPING SCHEME
    time_disc = 3;
    //=============================================================================
    Output::print<1>("LAPLACETYPE: ", 0, " NSTYPE's: ",1,", ",2,", ",3,", ",4, 
    " TIME_DISC: ", time_disc );
    errors[0] = {{0.0011233147920272, 0.014527988618419, 0.0081517719198929, 0.074300286739799}};
    errors[1] = {{0.0022562146863356, 0.029005332276073, 0.01401452192213, 0.14624345735295}};
    errors[2] = {{0.019173925569774, 0.24443063251283, 0.10235449879806, 1.2198025868096}};
    laplace_type = 0;
    check(domain, db, 12, -4711, 1, laplace_type, nl_form, time_disc, errors, tol);
    if(testall)
    {
      check(domain, db, 12, -4711, 2, laplace_type, nl_form, time_disc, errors, tol);
      check(domain, db, 12, -4711, 3, laplace_type, nl_form, time_disc, errors, tol);
    }
    check(domain, db, 12, -4711, 4, laplace_type, nl_form, time_disc, errors, tol);
  }
  parmoon::parmoon_finalize();
}
