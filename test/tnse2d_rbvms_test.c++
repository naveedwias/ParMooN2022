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
  else if(nonlinear_form ==1 )
    db["nse_nonlinear_form"] = "skew_symmetric";
  else if(nonlinear_form ==2 )
    db["nse_nonlinear_form"] = "rotational";
  else if(nonlinear_form ==3 )
    db["nse_nonlinear_form"] = "emac";
  else
    ErrThrow("Nonlinear form ", nonlinear_form," Not supported");

  // rough check
  if(time_disc != 3)
    ErrThrow("SUPG is only implemented for BDF2");

  db["time_discretization"]= "bdf_two";

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
  db["space_discretization_type"] = "supg";
  db["reynolds_number"] = 1e5;
  db["supg_delta0"] = 0.2;
  db["graddiv_stab"] = 0.2;
  
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
    int time_disc = 3;
    Output::print<1>("Testing the Q2/P1-disc elements");
    set_errors(errors);
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    //db["supg_delta0"] = 0.2;
    //db["graddiv_stab"] = 100;
    //check(domain_quads, db, 12, -4711, 4, lapacetype, 1, time_disc, errors, tol);
    db["example"] = 0;
    db["supg_delta0"] = 0.25;
    db["graddiv_stab"] = 0.25;
    errors[0] = {{0.0011033358362044, 0.014606538172967, 0.0073933510910158, 0.072530590832389}};
    errors[1] = {{0.00220873442159, 0.02915516583663, 0.014637455246588, 0.14708827476825}};
    errors[2] = {{0.0186511945546, 0.24656859592068, 0.44963304917506, 3.1357628660234}};   
    check(domain_quads, db, 12, -4711, 4, lapacetype, 0, time_disc, errors, tol);
    
    errors[0] = {{0.0011033358564588, 0.014606544096216, 0.0073925262970495, 0.072528253954498}};
    errors[1] = {{0.0022087339357477, 0.029155207879098, 0.014630835658039, 0.14706156149486}};
    errors[2] = {{0.018649712723659, 0.24659836364582, 0.44862454763272, 3.1273505315397}};
    nonlineartype = 1;
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    errors[0] = {{0.0011033619477922, 0.014606944442135, 0.0073962910146037, 0.072570393087597}};
    errors[1] = {{0.0022089635054129, 0.029158642738862, 0.014659661019378, 0.14726985966839}};
    errors[2] = {{0.018810439118855, 0.24906925801408, 0.45462757624864, 3.1500510693073}};
    nonlineartype = 2;
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    
    errors[0] = {{0.0011033642093273, 0.014606764391304, 0.0073904215892819, 0.072494614190416}};
    errors[1] = {{0.0022089366605954, 0.029156722697324, 0.014615470347597, 0.14691996338383}};
    errors[2] = {{0.01877597319626, 0.24739039074688, 0.44530462575901, 3.1221756618101}};
    nonlineartype=3;
    check(domain_quads, db, 12, -4711, 4, lapacetype, nonlineartype, time_disc, errors, tol);
    //Mostly the errors are same upto 6 to 7 decimel places
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
    int lapacetype = 0; int nonlineartype = 0; int time_disc = 3;
    
    set_errors(errors);    
    Output::print<1>("Testing the P2/P1 elements");
    check(domain_tri, db, 2, 1, 1, lapacetype, nonlineartype, time_disc, errors, tol);
    db["example"] = 0;
    errors[0] = {{0.00028488091109697, 0.0091858483384685, 0.0015346590648524, 0.03596886879566}};
    errors[1] = {{0.00057003816214114, 0.018347933250905, 0.0031449231788876, 0.071788241995252}};
    errors[2] = {{0.0048158023824515, 0.15464093290169, 0.027392896688359, 0.60529868722018}};   
    check(domain_tri, db, 2, 1, 4, lapacetype, nonlineartype, time_disc, errors, tol);
  }
  parmoon::parmoon_finalize();
}
