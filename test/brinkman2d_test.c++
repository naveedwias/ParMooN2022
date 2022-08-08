/**
 * @brief A test program for solving 2D Brinkman problems.
 *
 * This serves as a test for solving Brinkman problems in 2D. It is intended to
 * perform Brinkman calculations with different examples in different setups to test
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
 * @date 25/07/2016
 * @author Laura Blank
 *
 */

#include <Domain.h>
#include <Database.h>
#include <Brinkman2D.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <Example_Brinkman2D.h>
#include "LocalAssembling.h"
#include <BoundaryAssembling2D.h>
#include <MainUtilities.h> //for error measuring
#include <Chrono.h>
#include "ParMooN.h"
#include <algorithm>

#include <ParMooN_repository_info.h>

const std::string path = parmoon::source_directory;
const std::string path_to_repo = path + "/data/mesh/";

//----------------------------------------------------------------------------
// compare the computed errors in the Brinkman2D object with the given ones in
// the array
void compareErrors(const Brinkman2D& brinkman2d, std::array<double, 5> reference_errors)
{
  const double eps = 2e-9;

  Output::print(setprecision(14),reference_errors[0]);
  Output::print(setprecision(14),reference_errors[1]);
  Output::print(setprecision(14),reference_errors[2]);
  Output::print(setprecision(14),reference_errors[3]);
  Output::print(setprecision(14),reference_errors[4]);

  // check the errors
  if( std::abs(brinkman2d.getL2VelocityError() - reference_errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
        brinkman2d.getL2VelocityError(), " and ",  reference_errors[0]);
  }
  if( std::abs(brinkman2d.getL2DivergenceError() - reference_errors[1]) > eps)
  {
    ErrThrow("Program 1: L2 divergence velocity error not correct. ",
        brinkman2d.getL2DivergenceError() - reference_errors[1]);
  }
  if( std::abs(brinkman2d.getH1SemiVelocityError() - reference_errors[2]) > eps)
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
        brinkman2d.getH1SemiVelocityError() - reference_errors[2]);
  }
 if( std::abs(brinkman2d.getL2PressureError() - reference_errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct.",
        brinkman2d.getL2PressureError() - reference_errors[3]);
  }
  if( std::abs(brinkman2d.getH1SemiPressureError() - reference_errors[4]) > eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct.",
        brinkman2d.getH1SemiPressureError() - reference_errors[4]);
  }
}
//----------------------------------------------------------------------------
// compare the computed boundary errors in the Brinkman2D object with the given ones in
// the array
void compareAllErrors(const Brinkman2D& brinkman2d, std::array<double, 8> reference_errors)
{
  const double eps = 2e-9;

  Output::print(setprecision(14),reference_errors[0]);
  Output::print(setprecision(14),reference_errors[1]);
  Output::print(setprecision(14),reference_errors[2]);
  Output::print(setprecision(14),reference_errors[3]);
  Output::print(setprecision(14),reference_errors[4]);
  Output::print(setprecision(14),reference_errors[5]); //boundary error of u
  Output::print(setprecision(14),reference_errors[6]); // boundary error of u.n


  // check the errors
  if( std::abs(brinkman2d.getL2VelocityError() - reference_errors[0]) > eps )
  {
    ErrThrow("Program 1: L2 velocity error not correct. ",
        brinkman2d.getL2VelocityError(), " and ",  reference_errors[0]);
  }
  if( std::abs(brinkman2d.getL2DivergenceError() - reference_errors[1]) > eps)
  {
    ErrThrow("Program 1: L2 divergence velocity error not correct. ",
        brinkman2d.getL2DivergenceError() - reference_errors[1]);
  }
  if( std::abs(brinkman2d.getH1SemiVelocityError() - reference_errors[2]) > eps)
  {
    ErrThrow("Program 1: H1-semi velocity error not correct. ",
        brinkman2d.getH1SemiVelocityError() - reference_errors[2]);
  }
 
  if( std::abs(brinkman2d.getL2BoundaryError() - reference_errors[5]) > eps)
  {
    ErrThrow("Program 1: L2 velocity error at the boundary is not correct. ",
        brinkman2d.getL2BoundaryError() - reference_errors[5]);
  }
 if( std::abs(brinkman2d.getL2NormNormalComponentError() - reference_errors[6]) > eps)
  {
    ErrThrow("Program 1: L2 error of the normal velocity at the boundayr is not correct. ",
        brinkman2d.getL2NormNormalComponentError() - reference_errors[6]);
  }

 
 if( std::abs(brinkman2d.getL2PressureError() - reference_errors[3]) > eps )
  {
    ErrThrow("Program 1: L2 pressure error not correct.",
        brinkman2d.getL2PressureError() - reference_errors[3]);
  }
  if( std::abs(brinkman2d.getH1SemiPressureError() - reference_errors[4]) > eps )
  {
    ErrThrow("Program 1: H1-semi pressure error not correct.",
        brinkman2d.getH1SemiPressureError() - reference_errors[4]);
  }
}


//----------------------------------------------------------------------
// Here the actual computations take place
void check_brinkman2d(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
    std::array<double, 5> reference_errors, unsigned int nRefinements)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  TDatabase::ParamDB->PRESSURE_SPACE = pressureCode;

  Output::print("FEDATABASE");
  Output::set_outfile("Test.out");
  TDatabase::WriteParamDB((char*)"FEDATA");

  Brinkman2D brinkman2d(domain, db);
  brinkman2d.assemble();

/////////////// Routines for periodic boundary conditions /////////////////
if(db["example"].is(11))
{
	brinkman2d.findPeriodicDOFs();
	brinkman2d.checkPeriodicDOFs();
	brinkman2d.makePeriodicBoundary();
}

  brinkman2d.solve();
  brinkman2d.output(nRefinements);

  // compare computed with given errors
  compareErrors(brinkman2d, reference_errors); 
}

//----------------------------------------------------------------------
// Here the actual computations take place
// This version includes boundary errors
void check_brinkman2d_New(TDomain & domain, ParameterDatabase& db, int velocityCode,int pressureCode,
    std::array<double, 8> reference_errors, unsigned int nRefinements)
{
  TDatabase::ParamDB->VELOCITY_SPACE = velocityCode;
  TDatabase::ParamDB->PRESSURE_SPACE = pressureCode;

  Output::print("FEDATABASE");
  Output::set_outfile("Test.out");
  TDatabase::WriteParamDB((char*)"FEDATA");

  Brinkman2D brinkman2d(domain, db);
  brinkman2d.assemble();

/////////////// Routines for periodic boundary conditions /////////////////
if(db["example"].is(11))
{
	brinkman2d.findPeriodicDOFs();
	brinkman2d.checkPeriodicDOFs();
	brinkman2d.makePeriodicBoundary();
}

  brinkman2d.solve();
  brinkman2d.output(nRefinements);

  // compare computed with given errors
  compareAllErrors(brinkman2d, reference_errors); 
}



//################# EXAMPLE 1 - Exponential Flow (Hannukainen & Co) ####################################

void tests_on_triangles_P2P1_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ // default construct a domain object
  TDomain domain(db);
  // Initialization of the default parameters
  TDatabase::SetDefaultParameters(); 
  // refine grid up to the coarsest level
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;

  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
   
  db["viscosity"] = 0.004;
  db["PkPk_stab"] = false;
  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "symmetric Galerkin formulation";
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_nitsche_boundary = 0; //2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles (symmetric Galerkin formulation), Example 1, P2/P1, with Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.0058568178410394, 0.01431796219032, 0.14460835896882, 5.2469877709531e-05, 0.00061083280165148}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);

  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  Output::print("\nstarting with Brinkman2D on TwoTriangles (nonsymmetric Galerkin formulation), Example 1, P2/P1, with Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.0058568178410394, 0.01431796219032, 0.14460835896882, 5.2469877709531e-05, 0.00061083280165148}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements);
}

//======================================================================================================
void tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche_Example1(unsigned int nRefinements, ParameterDatabase& db)
{
  TDomain domain(db);
  TDatabase::SetDefaultParameters();

  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
 
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;
  db["PkPk_stab"] = false;
  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2; //2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 1, P2/P1, with penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{0.049394982619311, 0.049746071781476, 0.20444890666219, 9.3291684594484e-05, 0.0010263458082477}};
  check_brinkman2d(domain, db, 2,1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{
  TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
  
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 0; //2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{  1.2424692993215, 0.031419104874062, 16.668762604057, 0.0018335349292131, 0.034019179940349 }}; 
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Dirichlet and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{ 1.2347001456568, 1.3166645102027, 16.750021539147,0.018674085154587, 0.22010507412489}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ 
  TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;

  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";
  db["GradDiv_stab"] = false;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.192177714504, 5.8853778774783, 25.146773888043, 0.16699522666312, 3.1493307465155 }};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{9.599443219491, 4.1492978030163,  23.384929664058, 0.022130923948982, 0.27177748647831}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab100_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;
 
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.01;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.01;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 100;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by h_T";


  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1,3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;                     
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{ 10.869461373038, 0.0015555188075453, 24.270540284686, 0.1885807925005, 2.3578866662985 }};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.581773844742, 0.29752211831635, 24.022855667394, 0.10196998348134, 0.59373617246339}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;

  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 1.;
  db["effective_viscosity"] = 0.004;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{10.666137970649, 1.2062306735909, 24.166125351829, 0.21575466767375, 3.3261411856147}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example1, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and Neumann bcs and with visc_eff = visc = 0.004, perm = 1");
  reference_errors = {{9.8125894591994, 8.617299671234, 25.051206032843, 0.045806976649214, 0.32576244958451}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}



//############################## EXAMPLE 8 - SinCos_BadiaCodina ###########################################
//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;

  db["example"] = 8; // Poiseuille_Hannukainen
  db["permeability"] = 1;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 1");
  reference_errors = {{5.967728815198, 27.218645176007, 45.026548460105, 2.5622671365087, 33.995673079637}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (nonsymmetric GLS), Grad-Div stab (0.01), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 1");
  reference_errors = {{2.3166504721005, 24.851364206268, 27.233327688498, 0.00084168519671418, 0.010828366364366}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 
}

//======================================================================================================
void tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_smallK_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 5>  reference_errors;

  db["example"] = 8; // Poiseuille_Hannukainen
  db["permeability"] = 0.00001;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 0.004;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  TDatabase::ParamDB->l_T = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 84.168519671385, 1082.8366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001 ");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 93.875214873957, 2044.9301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  //l_T=1
  TDatabase::ParamDB->l_T = 1;
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";


  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 84.168519671385, 1082.8366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 93.875214873957, 2044.9301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.00001");
  reference_errors = {{5.6528468593088, 33.609960101853, 43.609769730799, 362.86296558258, 2403.2148842124}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 



  db["permeability"] = 0.01;
  //l_T=-1
  TDatabase::ParamDB->l_T = -1;

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01 ");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 0.093875214873957, 2.0449301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = -1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  //l_T=1
  TDatabase::ParamDB->l_T = 1;
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";


  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{2.3166504721003, 24.851364206268, 27.233327688498, 0.084168519671385, 1.0828366364365}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by L_0, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["equal_order_stab_scaling"] = "by h_T";
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{1.3167265706952, 21.848230410711, 28.026777395638, 0.093875214873957, 2.0449301910272}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 

  db["Galerkin_type"] = "symmetric Galerkin formulation";
  db["EqualOrder_PressureStab_type"] = "symmetric GLS";

  Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (symmetric GLS), Grad-Div stab (0.1), l_T = 1, scaling by h_T, penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 0.004, perm = 0.01");
  reference_errors = {{6.9403153847088, 27.940690277924, 68.678905054303, 3.0862494660051, 35.636962892309}};
  check_brinkman2d(domain, db, 1, 1, reference_errors, nRefinements); 


}

// includes boundary errors
void tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;

  db["example"] = 8; // SinCos BadiaCodina
  db["permeability"] = 0.001;
  db["effective_viscosity"] = 0.;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 1;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 0;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 4;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 1, 2, 3};
  TDatabase::ParamDB->nitsche_penalty = {0, 0, 0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 8, P1/P1-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");

  reference_errors = {{0.42679932017637, 11.721665354721, 15.463776092232, 36.359814532241, 1717.712382925, 1.3868446186387, 1.0439538851369}};

  check_brinkman2d_New(domain, db, 1, 1, reference_errors, nRefinements); 

}

// includes boundary errors
void tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;
 
  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 10;
  db["effective_viscosity"] = 0.00001;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 0;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 1, P1/P1-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");

  reference_errors = {{0.83048830850642, 6.418943841256, 79.759909017578, 0.00097536094080783, 0.028344716213816, 14.017340685057, 1.5519014816554}};

  check_brinkman2d_New(domain, db, 1, 1, reference_errors, nRefinements); 

}

// includes boundary errors
void tests_on_triangles_P2P2_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(unsigned int nRefinements, ParameterDatabase& db)
{ TDomain domain(db);
  TDatabase::SetDefaultParameters();
  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }
  std::array<double, 8>  reference_errors;

  db["example"] = 1; // Poiseuille_Hannukainen
  db["permeability"] = 10;
  db["effective_viscosity"] = 0.00001;
  db["viscosity"] = 1;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";
  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.1;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.1;
  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;
  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 0;

  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 2;
  TDatabase::ParamDB->neumann_boundary_id = {1, 3};
  TDatabase::ParamDB->neumann_boundary_value = {-0.5, 0.5};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {0, 2};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on TwoTriangles, Example 1, P2/P2-Stab (non-symmetric GLS) (0.1), Grad-Div stab (0.1), corner stab (1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with visc_eff = 0, visc = 1, perm = 0.001");

  reference_errors = {{0.81412759047597, 1.0492228233603, 89.89334987863, 9.1990083634641e-05, 0.0046073708563482, 12.146138140402, 0.092490981212866}};

  check_brinkman2d_New(domain, db, 2, 2, reference_errors, nRefinements); 

}


//################# EXAMPLE 11 - Riverbed ####################################
// includes boundary errors
void tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example11(unsigned int nRefinements, ParameterDatabase& db)
{ 

  db["boundary_file"].set(path_to_repo + "doubleRiverbed.PRM", false); // ( "../ParMooN/data/mesh/doubleRiverbed.PRM", false);
  db["geo_file"].set(path_to_repo + "doubleRiverbed3.mesh", false); // ("../ParMooN/data/mesh/doubleRiverbed3.mesh", false);

 db["example"] = 11; // Riverbed

  TDomain domain(db);
  TDatabase::SetDefaultParameters();

  for(unsigned int i = 0; i < nRefinements; i++)
  {
    domain.RegRefineAll();
  }

  std::array<double, 8>  reference_errors;

 
  db["permeability"] = -2.;
  db["effective_viscosity"] = 1.;
  db["viscosity"] = 1.;

  db["coefficient_function_type"] = 0;

  //Note that the parameters below have to be set in db AND TDatabase↲ 
  db["Galerkin_type"] = "nonsymmetric Galerkin formulation";

  db["PkPk_stab"] = true;
  db["equal_order_stab_weight_PkPk"] = 0.5;
  TDatabase::ParamDB->equal_order_stab_weight_PkPk = 0.5;

  db["GradDiv_stab"] = true;
  TDatabase::ParamDB->grad_div_stab_weight = 0.1;

  db["EqualOrder_PressureStab_type"] = "nonsymmetric GLS";
  db["equal_order_stab_scaling"] = "by L_0";
  TDatabase::ParamDB->L_0 = 0.1;

  db["corner_stab_weight"] = 0.1;

  //TDatabase::ParamDB->INPUT_QUAD_RULE = 99;
  
  TDatabase::ParamDB->n_neumann_boundary = 4;
  TDatabase::ParamDB->neumann_boundary_id = {0, 2, 3, 5};
  TDatabase::ParamDB->neumann_boundary_value = {0.001, 0,  0,  0.001};

  TDatabase::ParamDB->n_nitsche_boundary = 2;
  TDatabase::ParamDB->nitsche_boundary_id = {1, 4};
  TDatabase::ParamDB->nitsche_penalty = {0, 0};
  TDatabase::ParamDB->s1 = -1;
  TDatabase::ParamDB->s2 = -1;

  //l_T=-1
  //TDatabase::ParamDB->l_T = -1;

 Output::print("\nstarting with Brinkman2D on Riverbed, Example 11, P1/P1-Stab (non-symmetric GLS) (0.5), Grad-Div stab (0.1), corner stab (0.1), scaling by L_0 (0.1), penalty-free non-symmetric Nitsche approach and with physical parameters as set in the example file.");

  reference_errors = {{0.0045336169529337, 0.01367163604282, 0.041823496986892,  0.00057172495047007,  0.0010079231509394,  1.8300267399471e-05,  1.1929246910012e-06}};

  check_brinkman2d_New(domain, db, 1, 1, reference_errors, nRefinements); 

}



// ========================================================================
// =======================================================================
// main program
// =======================================================================
// ========================================================================
int main(int, char**)
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  // high order quadrature for computing errors
  TDatabase::ParamDB->INPUT_QUAD_RULE = 99;

  Output::setVerbosity(2);

  //---------------------------------------------------------------------
  // merge TDatabase with Problem specific ParameterDatabase db
  db.merge(ParameterDatabase::default_output_database(),true);
  db.merge(Example2D::default_example_database(),true);
 // db.merge(TDomain::default_domain_parameters(), true);  

  db.merge(Brinkman2D::get_default_Brinkman2D_parameters(),true);

  db["example"] = 1;
  db["output_compute_errors"] = true;
  db["output_write_vtk"] = false;

  db.add("boundary_file", "Default_UnitSquare", "");
  db.add("geo_file", "TwoTriangles", "", {"UnitSquare","Default_UnitSquare", "TwoTriangles"});

  //db.add("P1P1_stab", (bool) false, "" );
  //db.add("PkPk_stab", (bool) false, "", {true, false} );
  //db.add("equal_order_stab_weight_PkPk", (double) 0., "", (double) -1000, (double) 1000 );
  //db.add("refinement_n_initial_steps", (size_t) 2.0 , "", (size_t) 0, (size_t) 10000);



  tests_on_triangles_P2P1_Example1(2, db);

  tests_on_triangles_P2P1_PenaltyFreeNonSymmetricNitsche_Example1(2, db);

  tests_on_triangles_P1P1_GLSStab_Example1(2, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_Example1(2, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab100_Example1(2, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(2, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(2, db);

  tests_on_triangles_P1P1_GLSStab_PenaltyFreeNonSymmetricNitsche_GradDivStab_smallK_Example8(2, db);



// Tests including boundary errors
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example8(3, db);
tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(3, db);
tests_on_triangles_P2P2_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example1(3, db);

tests_on_triangles_P1P1_GLSStab_cornerstab_PenaltyFreeNonSymmetricNitsche_GradDivStab_Example11(0, db);

  parmoon::parmoon_finalize();
}
