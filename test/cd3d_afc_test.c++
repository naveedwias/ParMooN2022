/**
 * @brief A test program for the solving of Algebraic Flux Correction (AFC) in CD3D problems.
 *
 * This serves as a test for the solving of AFC schemes for CD3D problems. It is intended to
 * perform CD3D calculations with different examples in different setups to test
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

#include <MainUtilities.h> //for error measuring

// compare the computed errors in the CD3D object with the given ones in 
// the array
void compareErrors(const ConvectionDiffusion_AFC<3>& cd3d, std::array<double, 2> errors)
{
  const double eps = 1e-11;
  
  // check the errors
  if( std::abs(cd3d.get_L2_error() - errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 error not correct. ",
             cd3d.get_L2_error() - errors[0]);
  }
  if( std::abs(cd3d.get_H1_semi_error() - errors[1]) > eps )
  {
    ErrThrow("Program 1: H1-semi error not correct. ",
             cd3d.get_H1_semi_error() - errors[1]);
  }
}

// Here the actual computations take place
void check_cd3d(TDomain& domain, ParameterDatabase& db, std::array<double,2> errors)
{  
 
  // Choose and construct example.
  ConvectionDiffusion_AFC<3> cd3d(domain, db);
  cd3d.assemble(0);
  cd3d.solve(0);
  //maximum number of iterations for non lienar loop
  unsigned int max_it=1000;
  for(unsigned int k=1;;k++)
  {
    bool converged;
    converged = cd3d.solve(k);
    if ((converged)||(k>= max_it))
      break;
  }
  cd3d.output();
  //comparison of errors
  compareErrors(cd3d, errors); 
}

// =======================================================================
// main program
// =======================================================================
int main(int, char**)
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  Output::setVerbosity(2);
  Output::print("\ntesting with algebraic flux correction");
  
  db.merge(Example3D::default_example_database());
  db["example"] = 3; //Sharp Boundary Layer Example part of Example Database
  
  //only direct solver for CD3D class
  db.add("solver_type", std::string("direct"), "");
  db.add("refinement_n_initial_steps", (size_t) 3,"");
  db.add("multigrid_n_levels", (size_t) 0, "");
  
  //addition for AFC database
  db.merge(AlgebraicFluxCorrection::default_afc_database(), true);
  db["algebraic_flux_correction"].set("afc");
  db["diffusion_coefficient"]=1.0e-1;
  
  // default construct a domain object
  db.add("boundary_file", "Default_UnitCube", "");
  db.add("geo_file", "Default_UnitCube_Tetra", "", {"Default_UnitCube_Hexa", "Default_UnitCube_Tetra"});
  TDomain domain(db);
  
  db.info();
  
  domain.refine_and_get_hierarchy_of_collections(db);
  
  //AFC only applicable to P1 elements
  TDatabase::ParamDB->ANSATZ_ORDER = 1; //P1 elements
  //db["space_discretization_type"] = "galerkin"; //Galerkin Desicreitzation
  
  //maximum number of iterations for non linear loop
  db["afc_nonlinloop_maxit"]=1000;
 
  
  //Computation for P1 elements
  //P1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  std::array<double,2> errors = {{0.0010369623588627, 0.011894930265006}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
 
  errors = {{0.0010369624917919 , 0.011894928221049}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.0010369624396964, 0.011894928040666}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.0004164209086368, 0.010287281106288}};
  check_cd3d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.00041642068495158, 0.010287282108791}};
  check_cd3d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  errors = {{0.0004164206690843, 0.010287285812062}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  
  //P1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n ---------  P1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.0010369627652847, 0.011894929166125}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.0010369624917919, 0.011894928221049}};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.0010369624396964,0.011894928040666 }};
  check_cd3d(domain, db, errors);
  //=========================================================================
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.00041642072554058, 0.010287280993875}};
  check_cd3d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.00041642073256592, 0.010287281035637}};
  check_cd3d(domain, db, errors);
  //=========================================================================  
  Output::print("\n\n --------- P1 + BJK17 + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("BJK17");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.0004164208403835, 0.01028728162663}};
  check_cd3d(domain, db, errors);
  //=========================================================================

  //Computations for Q1 elements
  //Note: BJK17 limiter not applicable for Q1 elements.
  
  db["geo_file"].set<>("Default_UnitCube_Hexa");
  TDomain domain_quad(db);
  domain_quad.refine_and_get_hierarchy_of_collections(db);
  
  //Q1 elements and Dynamic Damping WITHOUT Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.00019851552281602, 0.0042946862986406}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.00019851553388223, 0.0042946863306704}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + no_anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("no");
  
  errors = {{0.00019851557114826, 0.0042946863131265}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  
  //Q1 elements and Dynamic Damping WITH Anderson Acceleration
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_rhs + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_rhs");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.00019851553307403, 0.0042946862979766}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + newton + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("newton");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.00019851553697763, 0.004294686306746}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  Output::print("\n\n --------- Q1 + kuzmin + fixed_point_matrix + anderson_acceleration ---------\n");
  db["afc_limiter"].set<>("kuzmin");
  db["afc_iteration_scheme"].set<>("fixed_point_matrix");
  db["afc_nonlinloop_anderson_acc"].set<>("yes");
  errors = {{0.00019851552155721, 0.0042946862932138}};
  check_cd3d(domain_quad, db, errors);
  //=========================================================================
  parmoon::parmoon_finalize();
}
