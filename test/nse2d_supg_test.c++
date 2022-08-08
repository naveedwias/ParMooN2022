/**
 * @brief A test program for the solving of NSE2D problems.
 *
 * This serves as a test for the solving of NSE2D problems. It is intended to
 * perform NSE2D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes wrong)
 * the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then this tests
 * shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 *
 * @author Naveed, Ulrich, Clemens
 *
 */
#include <cmath>

#include <Domain.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include "NavierStokes.h"
#include <Example_NSE2D.h>
#include <Multigrid.h>
#include <Chrono.h>
#include <algorithm>
#include "LocalAssembling.h"
#include "ParMooN.h"

#include <ParMooN_repository_info.h>

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/data/mesh/";

double accuracy = 1e-6;

void compare(const NavierStokes<2>& nse2d, std::array<double, int(5)> errors)
{
  auto computed_errors = nse2d.get_errors();
  
  // check the L2-error of the velcoity
  if( std::abs(computed_errors[0]-errors[0]) > accuracy )
  {
    ErrThrow("L2 norm of velocity: ", computed_errors[0], "  ", errors[0]);
  }
  // check the L2-error of the divergence of the velocity
  if( std::abs(computed_errors[1] - errors[1]) > accuracy )
  {
    ErrThrow("L2 norm of divergence: ", std::setprecision(14), computed_errors[1], "  ", errors[1]);
  }
  // check the H1-error of the velcoity
  if( std::abs(computed_errors[2] - errors[2]) > accuracy )
  {
    ErrThrow("H1 norm of velocity: ", computed_errors[2], "  ", errors[2]);
  }
  // check the L2-error of the pressure
  if( std::abs(computed_errors[3] - errors[3]) > accuracy)
  {
    ErrThrow("L2 norm of pressure: ", computed_errors[3], "  ", errors[3]);
  }
  // check the H1-error of the pressure
  if(std::abs(computed_errors[4] - errors[4]) > accuracy )
  {
    ErrThrow("H1 norm of pressure: ", computed_errors[4], "  ", errors[4]);
  }
}

void compute(TDomain &domain, ParameterDatabase& db,
             std::array<double, int(5)> errors)
{
  NavierStokes<2> nse2d(domain, db);
  nse2d.assemble_linear_terms();
  // check stopping criterion
  nse2d.stop_it(0);
  for(unsigned int k=1;; k++)
  {
    Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
                     nse2d.get_residuals());
    nse2d.solve();
    int lintype = TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR;
    if(lintype)
      break;
    // checking the first nonlinear iteration    
    nse2d.assemble_nonlinear_term();;
    if(nse2d.stop_it(k))
      break;
  }
  nse2d.output();
  // compare now the errors
  compare(nse2d, errors);
}

void check(TDomain &domain, ParameterDatabase db,
           int velocity_order, int pressure_order, 
           int nstype, int laplace_type,
           std::string nonlinear_form,
           std::array<double, int(5)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", nstype=", nstype, ", laplace_type=", laplace_type, 
                ", and nonlinear_form=", nonlinear_form);
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(Multigrid::default_multigrid_database());
  db["problem_type"] = 3;
  db["solver_type"] = "direct";
  db["iterative_solver_type"] = "fgmres";
  db["residual_tolerance"] = 1.e-12;
  db["preconditioner"] = "least_squares_commutator";
  
  db["nonlinloop_maxit"] = 50;
  db["nonlinloop_epsilon"] = 1e-10;

  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = pressure_order;
  TDatabase::ParamDB->NSTYPE = nstype;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  db["nse_nonlinear_form"] = nonlinear_form;
  
  Chrono timer;
  compute(domain, db, errors);  
}


// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
  /** Program 1
   *  test for SV element 
   */
  ParameterDatabase db = parmoon::parmoon_initialize();
  {
    db.merge(ParameterDatabase::default_nonlinit_database());
    db.merge(ParameterDatabase::default_output_database());
    db.merge(Example2D::default_example_database());
    db.merge(LocalAssembling2D::default_local_assembling_database());

    db["problem_type"]=3;
    db["example"] = 17;

    db.add("refinement_n_initial_steps", (size_t) 5,"");
    db.add("refinement_final_step_barycentric", true,"");

    db["nonlinloop_epsilon"] = 1e-10;

    db.add("boundary_file", "Default_UnitSquare", "");
    db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});

    // default construct a domain object
    TDomain domain(db);

    db["reynolds_number"] = 1e6;
    TDatabase::ParamDB->FLOW_PROBLEM_TYPE=3;
    db["space_discretization_type"] = "galerkin";
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = true;
    
    TDatabase::ParamDB->P1 = 0;
    TDatabase::ParamDB->NSTYPE = 4;
    // refine grid
    domain.refine_and_get_hierarchy_of_collections(db);
    std::array<double, int(5)> errors;

    //=========================================================================
    Output::print<1>("\nTesting the SV-2 elements");
    errors = {{ 0., 0., 0., 0.0045514002205762, 0.84192901372277 }};
    check(domain, db, 2, -110, 4, 0, "convective", errors);
  }// end program 1
  parmoon::parmoon_finalize();
}
