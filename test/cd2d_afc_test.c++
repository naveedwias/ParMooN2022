/**
 * @brief A test program for the solving of Algebraic Flux Correction (AFC) in CD2D problems.
 *
 * This serves as a test for the solving of AFC schemes for CD2D problems. It is intended to
 * perform CD2D calculations with different examples in different setups to test
 * a wide variety of ParMooN core functionality.
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
 * @date 2019/03/13
 * @author Abhinav Jha
 *
 */
#include <cmath>
#include <AlgebraicFluxCorrection.h>
#include <Domain.h>
#include <Database.h>
#include <ConvectionDiffusion_AFC.h>
#include <Multigrid.h>
#include "ParMooN.h"

#include <LocalAssembling.h>
#include <MainUtilities.h> //for error measuring

// compare the computed errors in the CD2D object with the given ones in 
// the array
void compareErrors(const ConvectionDiffusion_AFC<2>& cd2d, std::array<double, 4> errors)
{
  const double eps = 1e-8;
  
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
void check_cd2d(TDomain & domain, ParameterDatabase& db, std::array<double,4> errors)
{  
  ConvectionDiffusion_AFC<2> cd2d(domain, db);
  cd2d.assemble(0);
  cd2d.solve(0);
  //maximum number of iterations for non lienar loop
  unsigned int max_it=1000;
  for(unsigned int k=1;;k++)
  {
    bool converged;
    converged = cd2d.solve(k);
    if ((converged)||(k>= max_it))
      break;
  }
  cd2d.output();
  //comparison of errors
  compareErrors(cd2d, errors); 
}

// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  Output::setVerbosity(2);
  Output::print("\ntesting with algebraic flux correction");
  
  db.merge(Example2D::default_example_database());
  db["example"] = 3; //Sharp Boundary Layer Example part of Example Database
  
  //only direct solver for CD2D class
  db.add("solver_type", std::string("direct"), "");
  db.add("refinement_n_initial_steps", (size_t) 3,"");
  db.add("multigrid_n_levels", (size_t) 0, "");
  
  //addition for AFC database
  db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
  db["algebraic_flux_correction"].set("afc");
  
  // default construct a domain object
  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "TwoTriangles", "", {"UnitSquare", "TwoTriangles"});
  TDomain domain(db);
  
  //AFC only applicable to P1 elements
  TDatabase::ParamDB->ANSATZ_ORDER = 1; //P1 elements
  //db["space_discretization_type"] = "galerkin"; //Galerkin Desicreitzation
  
  //maximum number of iterations for non linear loop
  db["afc_nonlinloop_maxit"]=1000;
 
  // refine grid up to the coarsest level
  domain.refine_and_get_hierarchy_of_collections(db);
  
  //Computation for P1 elements
  //P1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  std::array<double,4> errors = {{0.53300159222151, 2.2679598767339,
    0.40567403714271,0.99999999999964 }};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.53300159151439, 2.2679598717208,
    0.40567403673183, 0.99999999999294}};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.533001591476, 2.2679598714175,
     0.40567403666082,0.99999999999266}};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.54553405607575, 2.5382226427268, 0.48949470018912, 1}};
  check_cd2d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{ 0.5455340557437, 2.5382226441468, 0.48949469803866, 1}};
  check_cd2d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.54553405576562, 2.5382226446816, 0.48949469801065, 1}};
  check_cd2d(domain, db, errors);
  //=========================================================================
  
  //P1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n ---------  P1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.53300159242554, 2.2679598769077, 0.40567403742442, 1.0000000009584 }};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.53300159208051, 2.2679598729401, 0.4056740371515, 1.0000000000071}};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.53300159216948, 2.2679598744846, 0.40567403724517,1.0000000000082}};
  check_cd2d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.54553405581158, 2.5382226459933, 0.48949469868059, 1.0000000001146}};
  check_cd2d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.54553405625277, 2.5382226482942, 0.48949469963447, 1.0000000000222}};
  check_cd2d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.54553405606112, 2.5382226468272, 0.48949469950613, 1.0000000000052}};
  check_cd2d(domain, db, errors);
  //=========================================================================

  //Computations for Q1 elements
  //Note: BJK17 limiter not applicable for Q1 elements.
  
  db["geo_file"].set<>("UnitSquare");
  TDomain domain_quad(db);
  // refine grid up to the coarsest level
  domain_quad.refine_and_get_hierarchy_of_collections(db);
  
  //Q1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.52464493314744, 2.2075436785016,
    0.38778102739911, 0.99998511032586}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.52464493317126, 2.2075436776853,
    0.38778102729849, 0.99998511033174}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.5246449330415, 2.2075436769812,
     0.38778102725893, 0.99998511032616}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  
  //Q1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.52464493309035, 2.2075436770111, 0.38778102708239, 0.99998511033224}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.52464493324716, 2.2075436781711, 0.38778102729767, 0.9999851103142}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.52464493319966, 2.2075436785091, 0.38778102729507, 0.99998511031822}};
  check_cd2d(domain_quad, db, errors);
  //=========================================================================
  parmoon::parmoon_finalize();
}
