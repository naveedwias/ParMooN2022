/**
 * @brief A test program for the solving of CD2D problems.
 *
 * This serves as a test for the solving of CD2D problems. It is intended to
 * perform CD2D calculations with different examples in different setups to test
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
 * @date 2015/11/11
 * @author Clemens Bartsch
 *
 */
#include <cmath>

#include <AlgebraicFluxCorrection.h>
#include <Domain.h>
#include <Database.h>
#include "FEDatabase.h"
#include "ConvectionDiffusion.h"
#include <Multigrid.h>
#include "all_defines_external_libraries.h"
#include "ParMooN.h"

#include <Chrono.h>
#include "LocalAssembling.h"
#include <MainUtilities.h> //for error measuring

// compare the computed errors in the Darcy2D object with the given ones in 
// the array
void compareErrors(const ConvectionDiffusion<2>& cd2d,
                   std::array<double, 4> errors)
{
  const double eps = 1e-11;
  
  // check the errors
  if( std::abs(cd2d.get_L2_error() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             cd2d.get_L2_error() - errors[0]);
  }
  if( std::abs(cd2d.get_H1_semi_error() - errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             cd2d.get_H1_semi_error() - errors[1]);
  }
  if( std::abs(cd2d.get_SD_error() - errors[2]) > eps )
  {
    ErrThrow("Program 1: sd error not correct.",
             cd2d.get_SD_error() - errors[2]);
  }
  if( std::abs(cd2d.get_L_inf_error() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L_inf error not correct.",
             cd2d.get_L_inf_error() - errors[3]);
  }
}

// Here the actual computations take place
void check_cd2d(TDomain & domain, const ParameterDatabase& db, int element_code,
                std::array<double, 4> errors)
{
  TDatabase::ParamDB->ANSATZ_ORDER = element_code;
  
  ConvectionDiffusion<2> cd2d(domain, db);
  cd2d.assemble();
  cd2d.solve();
  cd2d.output();
  // compare computed with given errors
  compareErrors(cd2d, errors); // throws upon a difference
}

void tests_on_quads(ParameterDatabase& db)
{
  db["geo_file"] = "UnitSquare";
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);

  std::array<double, 4> errors;
    
  Output::print("\nstarting with Q1");
  errors = {{ 0.030136750989447, 0.50136473756928, 0.50136473756928, 
              0.086620472046938 }};
  check_cd2d(domain, db, 1, errors);
  
  Output::print("\nstarting with Q2");
  errors = {{ 0.0019271846609828, 0.050971302154831, 0.050971302154831, 
              0.0033739148851507 }};
  check_cd2d(domain, db, 2, errors);
  
  Output::print("\nstarting with Q3");
  errors = {{ 8.8160702476897e-05, 0.0033762284959484, 0.0033762284959484,
              0.00033253078309192 }};
  check_cd2d(domain, db, 3, errors);
  
  Output::print("\nstarting with Q4");
  errors = {{ 3.3497247653192e-06, 0.00016700003852901, 0.00016700003852901,
              7.4417876531063e-06 }};
  check_cd2d(domain, db, 4, errors);
  
  if(db["solver_type"].is("iterative") && db["preconditioner"].is("multigrid"))
  {
    // the higher order elements do not/very slowly converge.
    // @todo find out why multigrid does not allow high order elements for CD2D
    return;
  }
  Output::print("\nstarting with Q5");
  errors = {{ 1.0746418844078e-07, 6.5922012588569e-06, 6.5922012588569e-06,
              4.0247144117433e-07 }};
  check_cd2d(domain, db, 5, errors);
  
  Output::print("\nstarting with Q6");
  errors = {{ 2.9751487986477e-09, 2.1654097614847e-07, 2.1654097614847e-07,
              6.0485347841421e-09 }};
  check_cd2d(domain, db, 6, errors);
  
  if(db["solver_type"].is("iterative"))
  {
    // The gmres residual is changing during restart.
    // Maybe the basis functions of these elements are not implemented with
    // high enough precision.
    return;
  }
  Output::print("\nstarting with Q7");
  errors = {{ 7.2354841489073e-11, 6.0911844559287e-09, 6.0911844559287e-09,
              2.7735080898594e-10 }};
  check_cd2d(domain, db, 7, errors);
  
  Output::print("\nstarting with Q8");
  errors = {{ 2.0882090117857e-12, 1.4998643575721e-10, 1.4998643575721e-10,
              5.4845017416483e-12 }};
  check_cd2d(domain, db, 8, errors);
}

void tests_on_triangles(ParameterDatabase& db, bool other_diagonal)
{
  db["geo_file"] = other_diagonal ? "TwoTriangles2" : "TwoTriangles";
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  domain.refine_and_get_hierarchy_of_collections(db);
  
  std::array<double, 4> errors;

  Output::print("\nstarting with P1");
  errors = {{ 0.068969534539772, 0.83003777179706, 0.83003777179706,
              0.16974869712877 }};
  check_cd2d(domain, db, 1, errors);
  
  Output::print("\nstarting with P2");
  errors = {{ 0.00421775248811, 0.12664905657584, 0.12664905657584,
              0.013686550806003 }};
  check_cd2d(domain, db, 2, errors);
  
  Output::print("\nstarting with P3");
  errors = {{ 0.0003326014029943, 0.013049293861338, 0.013049293861338,
              0.0013647682269053 }};
  check_cd2d(domain, db, 3, errors);
  
  Output::print("\nstarting with P4");
  errors = {{ 2.4012013491882e-05, 0.0011161787258859, 0.0011161787258859,
              8.9920045308348e-05 }};
  check_cd2d(domain, db, 4, errors);
  
  if(db["solver_type"].is("iterative") && db["preconditioner"].is("multigrid"))
  {
    // the higher order elements do not/very slowly converge.
    // @todo find out why multigrid does not allow high order elements for CD2D
    return;
  }
  
  Output::print("\nstarting with P5");
  errors = {{ 1.431573900717e-06, 7.8967826027302e-05, 7.8967826027302e-05,
              5.8590398090697e-06 }};
  check_cd2d(domain, db, 5, errors);
  
  Output::print("\nstarting with P6");
  errors = {{ 7.4027814432244e-08, 4.7847925721341e-06, 4.7847925721341e-06,
              3.2946208161633e-07 }};
  check_cd2d(domain, db, 6, errors);
}


// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
  auto db = parmoon::parmoon_initialize();
  Output::setVerbosity(2);
  
  size_t nRefinements = 2;
  
  //ParameterDatabase db = ParameterDatabase::parmoon_default_database();
  db.merge(Solver<>::default_solver_database());
  db.merge(Multigrid::default_multigrid_database());
  db.merge(Example2D::default_example_database());
  db.merge(ParameterDatabase::default_output_database());
  db["problem_type"] = 0; // problem type is not needed
  db["example"] = 0; // known sin-cos solution
  db["residual_tolerance"] = 1.0e-13;
  
  db["output_compute_errors"] = true;
  
  // from 'TDomain::default_domain_parameters()'
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "UnitSquare", "", {"UnitSquare", "TwoTriangles",
                                        "TwoTriangles2"});
  db.add("refinement_n_initial_steps", (size_t) nRefinements,"");
  
  Chrono time; // measure the time spent for each solver
  
  Output::print("\n\n ----------- direct solver -----------\n");
  db["solver_type"] = "direct";
  tests_on_quads(db);
  tests_on_triangles(db, false);
  tests_on_triangles(db, true);
  time.restart_and_print("all tests, direct solver");
  
#ifdef PARMOON_WITH_PETSC
  Output::print("\n\n ----------- petsc solver -----------\n");
  db["solver_type"] = "petsc";
  tests_on_quads(db);
  tests_on_triangles(db, false);
  tests_on_triangles(db, true);
  time.restart_and_print("all tests, petsc solver");
#endif // PARMOON_WITH_PETSC  
  
  
  Output::print("\n\n --------- fgmres+ssor solver ---------\n");
  db["solver_type"] = "iterative";
  db["iterative_solver_type"] = "fgmres";
  db["max_n_iterations"] = 1000;
  db["preconditioner"] = "ssor";
  tests_on_quads(db);
  tests_on_triangles(db, false);
  tests_on_triangles(db, true);
  time.restart_and_print("all tests, fgmres solver and ssor preconditioner");
  
  
  
  Output::print("\n\n --------- fgmres+multigrid solver ---------\n");
  db["iterative_solver_type"] = "fgmres";
  db["preconditioner"] = "multigrid";
  db["multigrid_n_levels"] = 2;
  db["multigrid_cycle_type"] = "W";
  db["multigrid_smoother"] = "jacobi";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_n_pre_smooth"] = 2;
  db["multigrid_n_post_smooth"] = 2;
  tests_on_quads(db);
  tests_on_triangles(db, false);
  tests_on_triangles(db, true);
  time.restart_and_print("all tests, fgmres solver and multigrid "
                         "preconditioner");
  
  parmoon::parmoon_finalize();
}
