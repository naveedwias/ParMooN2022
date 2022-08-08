/**
 * @brief A test program for the solving of Darcy problems in 2D.
 *
 * This serves as a test for the solving of Darcy problems. It is intended to
 * perform Darcy calculations with different examples in different setups to 
 * test a wide variety of ParMooN core functionality.
 * So far only one such test is implemented.
 *
 * The norms of the solution are compared with reference norms.
 * If those are not approximated well enough (or something in the process goes 
 * wrong) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.
 * Or you changed some program setup (e.g. changed the default solver). Then 
 * this tests shows you how many other program parts are affected by your 
 * changes. If you are not perfectly sure how to repair this, it is a good idea
 * to describe your changes in the forum and request support.
 *
 */

#include "all_defines_external_libraries.h"
#include <Domain.h>
#include <Database.h>
#include <Darcy.h>
#include "Multigrid.h"
#include <Chrono.h>
#include "ParMooN.h"
#include <cmath>

// compare the computed errors in the Darcy object with the given ones in 
// the array
template <int d>
void compareErrors(const Darcy<d>& darcy2d, std::array<double, 5> errors)
{
  const double eps = 2e-9;
  
  // check the errors
  if( std::abs(darcy2d.getL2VelocityError() - errors[0]) > eps 
     || std::isnan(darcy2d.getL2VelocityError()) )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
             darcy2d.getL2VelocityError() - errors[0]);
  }
  if( std::abs(darcy2d.getL2DivergenceError() - errors[1]) > 2*eps )
  {
    ErrThrow("Program 1: L2 velocity divergence error not correct. ",
             darcy2d.getL2DivergenceError() - errors[1]);
  }
  if( std::abs(darcy2d.getH1SemiVelocityError() - errors[2]) > eps )
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
             darcy2d.getH1SemiVelocityError() - errors[2]);
  }
  if( std::abs(darcy2d.getL2PressureError() - errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct. ",
             darcy2d.getL2PressureError() - errors[3]);
  }
  if( std::abs(darcy2d.getH1SemiPressureError() - errors[4]) > 2*eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct. ",
             darcy2d.getH1SemiPressureError() - errors[4]);
  }
}

// Here the actual computations take place
template <int d>
void check_darcy (TDomain & domain, ParameterDatabase& db, int velocityCode, 
           std::array<double, 5> errors)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  
  Darcy<d> darcy2d(domain, db);
  darcy2d.assemble();
  darcy2d.solve();
  darcy2d.output();
  // compare computed with given errors
  compareErrors<d>(darcy2d, errors); // throws upon a difference
}

std::string solver_description(const ParameterDatabase& db)
{
  if(db["solver_type"].is("direct"))
  {
    return "direct solver";
  }
  else if (db["solver_type"].is("petsc"))
  {
    return "petsc solver";
  }
  else
  {
    if(db["preconditioner"].is("least_squares_commutator"))
    {
      return "iterative solver with LSC";
    }
    else
    {
      return "iterative solver with SIMPLE";
    }
  }
}

#ifdef __2D__
void tests_on_quads(size_t nRefinements, ParameterDatabase db, bool testall)
{
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "UnitSquare";
  db["refinement_n_initial_steps"] = nRefinements;
  // default construct a domain object
  TDomain domain(db);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  std::string sd = solver_description(db);
  std::array<double, 5> errors;
    
  Output::print("\nstarting with RT0 on quads, ", sd);
  errors = {{ 2.1136884064519, 23.120239110875, 32.350359820686, 
              0.30277518654981, 4.4428829381584 }};
  check_darcy<2>(domain, db, 1000, errors);
  
  Output::print("\nstarting with RT1 on quads, ", sd);
  errors = {{ 0.40548747487874, 4.9461605478514, 12.695872686293, 
              0.063044878850364, 1.9735667622761 }};
  check_darcy<2>(domain, db, 1001, errors);
  
  if(testall)
  {
    Output::print("\nstarting with RT2 on quads, ", sd);
    errors = {{ 0.053370236091622, 0.66178203636242, 2.7529898961154,
                0.0084084949151706, 0.43469938937457 }};
    check_darcy<2>(domain, db, 1002, errors);
  }
  
  // this is quite slow!
  if(testall && db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with RT3 on quads, ", sd);
    errors = {{ 0.0052759466028633, 0.065778892133375, 0.39754725046775,
                0.00083489202575657, 0.063038523036227 }};
    check_darcy<2>(domain, db, 1003, errors);
  }
  
  Output::print("\nstarting with BDM1 on quads, ", sd);
  errors = {{ 1.7301620785317, 23.120239110875, 25.376165640701,
              0.32075488021636, 4.4428829381584 }};
  check_darcy<2>(domain, db, 1011, errors);
  
  Output::print("\nstarting with BDM2 on quads, ", sd);
  errors = {{ 0.36989913525384, 8.7083491818683, 8.7386341355059,
              0.1118419830706, 2.6029572935706 }};
  check_darcy<2>(domain, db, 1012, errors);
  
  // this is really slow (for the preconditioner
  // semi_implicit_method_for_pressure_linked_equations)
  if(testall && (db["solver_type"].is("direct")
                 || db["preconditioner"].is("least_squares_commutator")))
  {
    Output::print("\nstarting with BDM3 on quads, ", sd);
    errors = {{ 0.073739263678954, 2.2160565937362, 2.5477742459801,
                0.028107167396319, 0.99433882928285 }};
    check_darcy<2>(domain, db, 1013, errors);
  }
  
  Output::print("\nTests on quads, solution in ansatz space\n");
  //domain.RegRefineAll(); // the lower orders lead to singular matrices
  // tests, where the solution is in the ansatz space:
  db["example"] = 5;
  errors = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
  
  if(db["solver_type"].is("direct"))
  {
    // the direct solver does not solve correctly on the coarse mesh, after one
    // more refinement it seems to be ok (except for RT3, where another 
    // refinement is necessary)
    domain.RegRefineAll();
  }
  
  db["is_RT_polynomial"] = false;
  Output::print("\nTest BDM1 on quads, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 1;
  check_darcy<2>(domain, db, 1011, errors);
  Output::print("\nTest BDM2 on quads, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 2;
  check_darcy<2>(domain, db, 1012, errors);
  if(testall)
  {
    Output::print("\nTest BDM3 on quads, solution in ansatz space, ", sd);
    db["degree_polynomial"] = 3;
    check_darcy<2>(domain, db, 1013, errors);
  }
  
  db["is_RT_polynomial"] = true;
  Output::print("\nTest RT0 on quads, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 0;
  check_darcy<2>(domain, db, 1000, errors);
  Output::print("\nTest RT1 on quads, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 1;
  check_darcy<2>(domain, db, 1001, errors);
  if(testall)
  {
    Output::print("\nTest RT2 on quads, solution in ansatz space, ", sd);
    db["degree_polynomial"] = 2;
    check_darcy<2>(domain, db, 1002, errors);
    if(db["solver_type"].is("direct"))
    {
      domain.RegRefineAll();
    }
    Output::print("\nTest RT3 on quads, solution in ansatz space, ", sd);
    db["degree_polynomial"] = 3;
    check_darcy<2>(domain, db, 1003, errors);
  }
}

void tests_on_triangles(size_t nRefinements, ParameterDatabase db, bool testall)
{
  db["boundary_file"] = "Default_UnitSquare";
  db["geo_file"] = "TwoTriangles";
  db["refinement_n_initial_steps"] = nRefinements;
  // default construct a domain object
  TDomain domain(db);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  std::string sd = solver_description(db);
  std::array<double, 5> errors;

  Output::print("\nstarting with RT0 on triangles, ", sd);
  errors = {{ 2.0190538085226, 19.176327642953, 31.034472755421,
              0.24829883521652, 4.4428829381584 }};
  check_darcy<2>(domain, db, 1000, errors);
  
  Output::print("\nstarting with RT1 on triangles, ", sd);
  errors = {{ 0.45588694172323, 5.772807889241, 14.551870959373,
              0.074463174420845, 2.3252470856326 }};
  check_darcy<2>(domain, db, 1001, errors);
  
  // this is really slow!
  if(testall && (db["solver_type"].is("direct")
                 || db["preconditioner"].is("least_squares_commutator")))
  {
    Output::print("\nstarting with RT2 on triangles, ", sd);
    errors = {{ 0.079823667222536, 1.2850584737743, 4.1870596885453,
                0.016309057144146, 0.81069387503794 }};
    check_darcy<2>(domain, db, 1002, errors);
  }
  
  // this is really slow!
  if(testall && db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with RT3 on triangles, ", sd);
    errors = {{ 0.011137525370491, 0.22596916693238, 0.83882076917212,
                0.0028837387982513, 0.2072548767011 }};
    check_darcy<2>(domain, db, 1003, errors);
  }
  
  Output::print("\nstarting with BDM1 on triangles, ", sd);
  errors = {{ 1.3463140129828, 19.176327381503, 24.553524052795,
              0.25993100617474, 4.4428829381584 }};
  check_darcy<2>(domain, db, 1011, errors);
  
  Output::print("\nstarting with BDM2 on triangles, ", sd);
  errors = {{ 0.22616033982856, 5.772807889241, 8.2103703304275,
              0.074718685375952, 2.3561242018165 }};
  check_darcy<2>(domain, db, 1012, errors);
  
  // this is really slow!
  if(testall && db["solver_type"].is("direct"))
  {
    Output::print("\nstarting with BDM3 on triangles, ", sd);
    errors = {{ 0.037364969838055, 1.2850584737743, 1.8808225137313,
                0.01629354826757, 0.81659172891947 }};
    check_darcy<2>(domain, db, 1013, errors);
  }
  
  
  Output::print("\nTests on triangles, solution in ansatz space\n");
  //domain.RegRefineAll(); // the lower orders lead to singular matrices
  // tests, where the solution is in the ansatz space:
  db["example"] = 5;
  errors = {{ 0.0, 0.0, 0.0, 0.0, 0.0 }};
  
  if(db["solver_type"].is("direct"))
  {
    // the direct solver does not solve correctly on the coarse mesh, after one
    // more refinement (3 in total) it seems to be ok for BDM1 and BDM2. For the
    // other elements even more refinements are necessary: BMD3 (6), RT0 (5),
    // RT1 (4), RT2 (5), RT3 (>6?)
    domain.RegRefineAll();
  }
  
  db["is_RT_polynomial"] = false;
  Output::print("\nTest BDM1 on triangles, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 1;
  check_darcy<2>(domain, db, 1011, errors);
  Output::print("\nTest BDM2 on triangles, solution in ansatz space, ", sd);
  db["degree_polynomial"] = 2;
  check_darcy<2>(domain, db, 1012, errors);
  if(!db["solver_type"].is("direct"))
  {
    if(testall && !db["solver_type"].is("iterative")) // only petsc
    {
      Output::print("\nTest BDM3 on triangles, solution in ansatz space, ", sd);
      db["degree_polynomial"] = 3;
      check_darcy<2>(domain, db, 1013, errors);
    }
    
    db["is_RT_polynomial"] = true;
    Output::print("\nTest RT0 on triangles, solution in ansatz space, ", sd);
    db["degree_polynomial"] = 0;
    check_darcy<2>(domain, db, 1000, errors);
    Output::print("\nTest RT1 on triangles, solution in ansatz space, ", sd);
    db["degree_polynomial"] = 1;
    check_darcy<2>(domain, db, 1001, errors);
    if(testall && !db["solver_type"].is("iterative")) // only petsc
    {
      Output::print("\nTest RT2 on triangles, solution in ansatz space, ", sd);
      db["degree_polynomial"] = 2;
      check_darcy<2>(domain, db, 1002, errors);
      Output::print("\nTest RT3 on triangles, solution in ansatz space, ", sd);
      db["degree_polynomial"] = 3;
      check_darcy<2>(domain, db, 1003, errors);
    }
  }
}
#endif // 2D

#ifdef __3D__
void tests_on_hexas(size_t nRefinements, ParameterDatabase& db, bool testall)
{
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Hexa";
  db["refinement_n_initial_steps"] = nRefinements;
  // default construct a domain object
  TDomain domain(db);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  std::array<double, 5> errors;
  std::string sd = solver_description(db);
  
  Output::print("\nstarting with RT1 on hexas, ", sd);
  if(nRefinements == 2)
  {
    errors = {{ 0.081377821454114, 0.41305643219106, 1.9531296260514, 
                0.020181597012015, 0.43928918412076 }}; // finer grid
  }
  else
  {
    errors = {{ 0.31882716748202, 1.6031495595555, 3.9689055156677,
                0.076647168024276, 0.88660895798807 }}; // coarser grid
  }
  check_darcy<3>(domain, db, 1001, errors);
  
  // for some unknown reason umfpack does not solve these systems correctly 
  // (leading to large residuals). I don't know if the error is in ParMooN or
  // not, probably in ParMooN and related to Dirichlet dofs. On a finer grid
  // there is no problem.
  if(testall && !db["solver_type"].is("direct") 
             && !db["preconditioner"].is("multigrid"))
  {
    Output::print("\nstarting with RT2 on hexas, ", sd);
    if(nRefinements == 2)
    {
       errors = {{ 0.0042461302600703, 0.027465355692105, 0.21305934724293,
                   0.00093095519313767, 0.048205238919775 }}; // finer grid
    }
    else
    {
      errors = {{ 0.034460131770028, 0.21491247109762, 0.83587773101988,
                  0.0073649888125504, 0.19044157847877 }}; // coarser grid
    }
    check_darcy<3>(domain, db, 1002, errors);
  
    Output::print("\nstarting with BDM2 on hexas, ", sd);
    if(nRefinements == 2)
    {
      errors = {{ 0.049959138131464, 0.97963312485102, 1.1592414284114,
                  0.033440305109556, 0.7293026775163 }}; // finer grid
    }
    else
    {
      errors = {{ 0.43545926716285, 3.4654312162412, 4.590648939376,
                  0.11734467078918, 1.3177499679789 }}; // coarser grid
    }
    check_darcy<3>(domain, db, 1012, errors);
  }
  
  // this is really slow and/or not working
  if(testall && !db["solver_type"].is("iterative"))
  {
    Output::print("\nstarting with BDM3 on hexas, ", sd);
    //errors = {{ 0.0089705896498808, 0.17583623566595, 0.22177144929476,
    //            0.0060922146271361, 0.18113542308142 }}; // finer grid
    errors = {{ 0.15931272192534, 1.2466839070804, 2.1075592487452,
                0.044811802526853, 0.64542130527488 }}; // coarser grid
    check_darcy<3>(domain, db, 1013, errors);
  }
  
  // for RT0 and BDM1 we need one more refinement step
  domain.RegRefineAll();
  Output::print("\nstarting with RT0 on hexas, ", sd);
  if(nRefinements == 2)
  {
    errors = {{ 0.30780429929604, 2.0370781691994, 8.6278645931962,
                0.06894139176197,  1.9238247452428 }}; // finer grid
  }
  else
  {
    errors = {{ 0.61129533426472, 3.966863809562, 8.8488505583356, 
                0.13495292835257, 1.9238247452428 }}; // coarser grid
  }
  check_darcy<3>(domain, db, 1000, errors);
  
  Output::print("\nstarting with BDM1 on hexas, ", sd);
  if(nRefinements == 2)
  {
    errors = {{ 0.087304340192265, 2.0370781691994, 2.5106627579722,
                0.069621362924773, 1.9238247452428 }}; // finer grid
  }
  else
  {
    errors = {{ 0.3387619514894, 3.9668638095605, 4.9771676096652, 
                0.1386528493848, 1.9238247452428 }}; // coarser grid
  }
  check_darcy<3>(domain, db, 1011, errors);
}

void tests_on_tetras(size_t nRefinements, ParameterDatabase& db,
                     bool /*testall*/)
{
  db["boundary_file"] = "Default_UnitCube";
  db["geo_file"] = "Default_UnitCube_Tetra";
  db["refinement_n_initial_steps"] = nRefinements;
  // default construct a domain object
  TDomain domain(db);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  std::array<double, 5> errors;
  std::string sd = solver_description(db);

  /// @todo tests for darcy in tetrahedra for higher order H(div)-elements
  
  Output::print("\nstarting with RT0 on tetras, ", sd);
  if(nRefinements == 2)
  {
    errors = {{ 0.59317475308108, 2.9288678573989, 8.7129926338484,
                0.10463892759756, 1.9240537502598 }}; // finer grid
  }
  else
  {
    errors = {{ 1.1967906038354, 5.3236531054449, 9.0834898178534,
                0.23103198049002, 1.9154448537963 }}; // coarser grid
  }
  check_darcy<3>(domain, db, 1000, errors);
  
  //Output::print("\nstarting with RT1 on tetras, ", sd);
  //errors = {{ 0.31882716748202, 1.6031495595555, 3.9689055156677, 
  //            0.076647168024276, 0.88660895798807 }};
  //check_darcy<3>(domain, db, 1001, errors);
  
  //Output::print("\nstarting with RT2 on tetras, ", sd);
  //errors = {{ 0., 0., 0., 0., 0. }};
  //check_darcy<3>(domain, db, 1002, errors);
  
  //Output::print("\nstarting with RT3 on tetras, ", sd);
  //errors = {{ 0., 0., 0., 0., 0. }};
  //check_darcy<3>(domain, db, 1003, errors);
  
  Output::print("\nstarting with BDM1 on tetras, ", sd);
  if(nRefinements == 2)
  {
    errors = {{ 0.23560384814218, 2.9288678573989, 4.7433456980206, 
                0.10489178373724, 1.9240537502598 }}; // finer grid
  }
  else
  {
    errors = {{ 0.78119356117933, 5.3236531054449, 8.0239757519786,
                0.21286638433987, 1.9154448537963 }}; // coarser grid
  }
  check_darcy<3>(domain, db, 1011, errors);
  
  //Output::print("\nstarting with BDM2 on tetras, ", sd);
  //errors = {{ 0., 0., 0., 0., 0. }};
  //check_darcy<3>(domain, db, 1012, errors);
  
  //Output::print("\nstarting with BDM3 on tetras, ", sd);
  //errors = {{ 0., 0., 0., 0., 0. }};
  //check_darcy<3>(domain, db, 1013, errors);
}
#endif // 3D

void tests_on_different_grids(size_t nRefinements, ParameterDatabase& db,
                              bool testall)
{
#ifdef __2D__
  tests_on_quads(nRefinements, db, testall);
  tests_on_triangles(nRefinements, db, testall);
#else
  if(!db["preconditioner"].is("multigrid"))
    nRefinements--; // save some computational time
  tests_on_hexas(nRefinements, db, testall);
  tests_on_tetras(nRefinements, db, testall);
#endif
}

// =======================================================================
// main program
// =======================================================================
int main(int argc, char*argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  bool testall = false;
  if(argc > 1 && argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }
  
#ifdef __2D__
  using Example = Example2D;
#else
  using Example = Example3D;
#endif
  
  // velocity space code for Raviart-Thomas (RT) and 
  // Brezzi-Douglas-Marini(BDM) elements:
  // 1000    RT_0
  // 1001    RT_1
  // 1002    RT_2
  // 1003    RT_3
  // 1011    BDM_1
  // 1012    BDM_2
  // 1013    BDM_3
  TDatabase::ParamDB->VELOCITY_SPACE = 1000;
  // automatically choose pressure space to get inf-sup stable pair
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
    // high order quadrature for computing errors
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  TDatabase::ParamDB->SIGMA_PERM = 1.; // permeability

  size_t nRefinements = 2;
  
  Output::setVerbosity(2);
  
  db.merge(Solver<>::default_solver_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(TDomain::default_domain_parameters());
  db.merge(Multigrid::default_multigrid_database());
  db.merge(Example::default_example_database());
  db["problem_type"] = 0; // problem type is not needed
  db["example"] = 0; // known sin-cos solution
  db["residual_tolerance"] = 1.0e-11;
  db["boundary_file"] = "Default_UnitSquare";
  
  db["output_compute_errors"] = true;
  
  Chrono timer;
  
  Output::print("\n\n ----------- direct solver -----------\n");
  db["solver_type"] = "direct";

  tests_on_different_grids(nRefinements, db, testall);
  timer.restart_and_print("all tests, the direct solver");
  
  
  db["max_n_iterations"] = 10000;
#ifdef PARMOON_WITH_PETSC
  Output::print("\n\n ----------- PETSc solver -----------\n");
  db["solver_type"] = "petsc";

  std::string petsc_args = "-ksp_monitor"
      " -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur"
      " -fieldsplit_0_ksp_atol 1.0e-12 -fieldsplit_0_ksp_rtol 0."
      " -fieldsplit_1_ksp_atol 1.0e-12 -fieldsplit_1_ksp_rtol 0.";

  db["petsc_arguments"].impose(Parameter("petsc_arguments", petsc_args, ""));

  tests_on_different_grids(nRefinements, db, testall);
  timer.restart_and_print("all tests, the petsc solver");
#endif // PARMOON_WITH_PETSC
  
  db["solver_type"] = "iterative";
  
  Output::print("\n\n --------- fgmres+lsc solver ---------\n");
  db["preconditioner"] = "least_squares_commutator";
  tests_on_different_grids(nRefinements, db, testall);
  timer.restart_and_print("all tests, fgmres with lsc preconditioning");
  
  /// @todo test fgmres+simple for darcy problems in 3D does not work
#ifdef __2D__
  Output::print("\n\n -------- fgmres+simple solver -------\n");
  db["preconditioner"] = "semi_implicit_method_for_pressure_linked_equations";
  tests_on_different_grids(nRefinements, db, testall);
  timer.restart_and_print("all tests, fgmres with simple preconditioning");
#endif
  
  //Output::print("\n\n --------- fgmres+multigrid solver ---------\n");
  //db["preconditioner"] = "multigrid";
  //db["multigrid_n_levels"] = nRefinements;
  //db["multigrid_type"] = "standard";
  //db["multigrid_cycle_type"] = "V";
  //db["multigrid_smoother"] = "cell_vanka";
  //tests_on_different_grids(nRefinements, db, testall);
  //timer.restart_and_print("all tests, fgmres with multigrid preconditioning");
  
  timer.print_total_time("all darcy tests");
  parmoon::parmoon_finalize();
}
