/**
 * @brief A test program for a time dependent CD2D problem,
 * algorithm "FEM-FCT"
 *
 * @author Clemens Bartsch, Naveed Ahmed, Ondrej Partl
 */
#include <algorithm>

#include <cmath>

#include "all_defines_external_libraries.h"
#include "TimeConvectionDiffusion.h"
#include <Database.h>
#include <TimeDiscRout.h>
#include <MainUtilities.h>
#include <ConvDiff.h>
#include <AuxParam2D.h>
#include <TimeDiscretizations.h>
#include "ParMooN.h"

#include <ParMooN_repository_info.h>
const std::string path = parmoon::source_directory;
const std::string rectangle_PRM_file = path + "/data/mesh/" + "rectangle_3x1.PRM";
const std::string rectangle_mesh_file = path + "/data/mesh/" + "rectangle_3x1.mesh";

enum TestNames
{
  LinearWithDirichletCrankNicolson,
  LinearWithDirichletNeumannCrankNicolson,
  LinearWithDirichletEuler,
  LinearWithDirichletNeumannEuler,
  NonLinearWithDirichletCrankNicolson,
  NonLinearWithDirichletNeumannCrankNicolson,
  NonLinearWithDirichletEuler,
  NonLinearWithDirichletNeumannEuler,
  RotatingBodies, 
  RotatingBodiesWithPETSC,
  MonolithicWithDirichlet,
  MonolithicWithDirichletNeumann,
  MolNonlinearWithDirichletCrankNicolson,
  MolNonLinearWithDirichletNeumannCrankNicolson,
  ZalExplicitWithDirichlet,
  ZalExplicitWithDirichletNeumann
};

void testCN(TimeConvectionDiffusion<2> &tcd, int m);
void timeIntegration(TimeConvectionDiffusion<2>& tcd,
                     ParameterDatabase errorDatabase);
void setCommonParameters(ParameterDatabase& db);
void runTest(TestNames testName,
             ParameterDatabase& db,
             ParameterDatabase& errorDatabase,
             std::vector<std::string>& errorNames);
void setTestSpecificParameters(TestNames testName,
                               ParameterDatabase& db,
                               ParameterDatabase& errorDatabase);
void prepareErrorDatabase(ParameterDatabase& errorDatabase,
                          std::vector<std::string>& errorNames);
void checkErrors(TimeConvectionDiffusion<2> &tcd,
                 ParameterDatabase& errorDatabase,
                 std::vector<std::string>& errorNames);
bool differentErrors(double measuredError, double expectedError);


int main(int, char**)
{
  parmoon::parmoon_initialize();
  ParameterDatabase db = TimeConvectionDiffusion<2>::default_tcd_database(true);

  std::vector<std::string>
  errorNames = { "L2(0,T;L2(c))", "L2(0,T;H1(c))", "L_inf(0,T,L_inf(c))" };
  
  ParameterDatabase errorDatabase("database of errors");
  prepareErrorDatabase(errorDatabase, errorNames);
  
  setCommonParameters(db);
  
  runTest(LinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletEuler, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletNeumannEuler, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletEuler, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannEuler, db, errorDatabase, errorNames);
  
  runTest(MonolithicWithDirichlet, db, errorDatabase, errorNames);
  runTest(MonolithicWithDirichletNeumann, db, errorDatabase, errorNames);
  runTest(MolNonlinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(MolNonLinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  
  runTest(ZalExplicitWithDirichlet, db, errorDatabase, errorNames);
  runTest(ZalExplicitWithDirichletNeumann, db, errorDatabase, errorNames);
  
  runTest(RotatingBodies, db, errorDatabase, errorNames);
  
  #ifdef PARMOON_WITH_PETSC
  runTest(RotatingBodiesWithPETSC, db, errorDatabase, errorNames);
  #endif
  parmoon::parmoon_finalize();
}


void prepareErrorDatabase(ParameterDatabase& errorDatabase,
                          std::vector<std::string>& errorNames)
{
  errorDatabase.add("checkOnlyFinalErrors", true, "", {true, false});
  
  for( unsigned int i = 0; i < errorNames.size(); i++)
    errorDatabase.add(errorNames.at(i), 0.0, "", 0.0, 1.0e3);
}

void setCommonParameters(ParameterDatabase& db)
{
  db["algebraic_flux_correction"].set("fem-fct-cn");
  db["space_discretization_type"] = "galerkin";
  TDatabase::ParamDB->ANSATZ_ORDER = 1;
}

void runTest(TestNames testName,
             ParameterDatabase& db,
             ParameterDatabase& errorDatabase,
             std::vector<std::string>& errorNames)
{
  setTestSpecificParameters(testName, db, errorDatabase);
  
  TDatabase::TimeDB->TIMESTEPLENGTH = db["time_step_length"];
  TDatabase::TimeDB->CURRENTTIME = db["time_start"];
  
  TDomain domain(db);
  SetTimeDiscParameters(0);
  domain.refine_and_get_hierarchy_of_collections(db);
  
  TimeConvectionDiffusion<2> tcd(domain, db);
  
  timeIntegration(tcd, errorDatabase);
  
  if(errorDatabase["checkOnlyFinalErrors"])
    checkErrors(tcd, errorDatabase, errorNames);
}

void setTestSpecificParameters(TestNames testName,
                               ParameterDatabase& db,
                               ParameterDatabase& errorDatabase)
{
  switch(testName)
  {
    case RotatingBodies:
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(0);
    db["boundary_file"].set("Default_UnitSquare", false);
    db["geo_file"].set("UnitSquare", false);
    db["refinement_n_initial_steps"].set(5);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.;
    db["time_end"] = 0.02;
    db["time_step_length"] = 0.001;
    errorDatabase["checkOnlyFinalErrors"] = false;
    break;
    
    case RotatingBodiesWithPETSC:
    db["example"] = 3;
    db["solver_type"] = "petsc";
    db["afc_prelimiter"].set(0);
    db["boundary_file"].set("Default_UnitSquare", false);
    db["geo_file"].set("UnitSquare", false);
    db["refinement_n_initial_steps"].set(5);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.;
    db["time_end"] = 0.02;
    db["time_step_length"] = 0.001;
    errorDatabase["checkOnlyFinalErrors"] = false;
    break;
    
    case LinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(0);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.0076402599936465;
    errorDatabase["L2(0,T;H1(c))"] = 0.027778375069014;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.020803672151631;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1/4.0;
    errorDatabase["L2(0,T;L2(c))"] = 0.044128398660717;
    errorDatabase["L2(0,T;H1(c))"] = 0.40873110156187;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.24555341708065;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(0);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.0046988533818225;
    errorDatabase["L2(0,T;H1(c))"] = 0.017090987764126;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.010031317180298;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletNeumannEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1/4.0;
    errorDatabase["L2(0,T;L2(c))"] = 0.035608714678536;
    errorDatabase["L2(0,T;H1(c))"] = 0.40968163793062;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.25572324713511;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 1.2388213044174e-08;
    errorDatabase["L2(0,T;H1(c))"] = 1.1574788676608e-07;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 3.2199690033289e-16;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.050719397170556;
    errorDatabase["L2(0,T;H1(c))"] = 0.58834431753712;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.22594131547883;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.01;
    errorDatabase["L2(0,T;L2(c))"] = 6.6579918021874e-06;
    errorDatabase["L2(0,T;H1(c))"] = 4.7044953971594e-05;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.396084381744e-05;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.056454364548808;
    errorDatabase["L2(0,T;H1(c))"] = 0.62178733853338;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.28532896602492;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MonolithicWithDirichlet:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("monolithic");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.05;
    errorDatabase["L2(0,T;L2(c))"] = 0.037542358263816;
    errorDatabase["L2(0,T;H1(c))"] = 0.2400110876656;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.10967297224676;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MonolithicWithDirichletNeumann:
    db["diffusion_coefficient"] = 1e-3;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("monolithic");
    db["admissible_minimum"] = 0.0;
    db["admissible_maximum"] = 3.0;
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.025;
    errorDatabase["L2(0,T;L2(c))"] = 0.044760538496678;
    errorDatabase["L2(0,T;H1(c))"] = 0.32752167395021;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.15683455981275;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MolNonlinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("monolithic");
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.05;
    errorDatabase["L2(0,T;L2(c))"] = 0.037508462972126;
    errorDatabase["L2(0,T;H1(c))"] = 0.23993275914316;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.10967301210719;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MolNonLinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-3;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("monolithic");
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.05;
    errorDatabase["L2(0,T;L2(c))"] = 0.044773933802125;
    errorDatabase["L2(0,T;H1(c))"] = 0.32771109884545;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.15683406503625;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case ZalExplicitWithDirichlet:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 8;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("zalesak");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.05;
    errorDatabase["L2(0,T;L2(c))"] = 0.00035085459758354;
    errorDatabase["L2(0,T;H1(c))"] = 0.0063954701516334;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.0063004913475045;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case ZalExplicitWithDirichletNeumann:
    db["diffusion_coefficient"] = 1e-3;
    db["example"] = 9;
    db["solver_type"] = "direct";
    db["boundary_file"].set(rectangle_PRM_file, false);
    db["geo_file"].set(rectangle_mesh_file, false);
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("zalesak");
    db["admissible_minimum"] = 0.0;
    db["admissible_maximum"] = 3.0;
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.025;
    errorDatabase["L2(0,T;L2(c))"] = 0.015297147298357;
    errorDatabase["L2(0,T;H1(c))"] = 0.25604935246255;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.12055078702737;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    default:
      ErrThrow("Unspecified test");
  }
}

void timeIntegration(TimeConvectionDiffusion<2>& tcd,
                     ParameterDatabase errorDatabase)
{
  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.set_time_disc_parameters();
  
  tcd.assemble_initial_time();
  int step = 0;
  
  bool checkOnlyFinalErrors = errorDatabase["checkOnlyFinalErrors"];
  
  if(checkOnlyFinalErrors)
    tcd.output();
  else
    testCN(tcd, step);

  double end_time = tcd.get_db()["time_end"];
  while(TDatabase::TimeDB->CURRENTTIME < end_time -1e-10)
  {
    tss.nonlin_iteration = 0;
    step ++;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);

    TDatabase::TimeDB->CURRENTTIME += TDatabase::TimeDB->TIMESTEPLENGTH;
    tss.current_time_ = TDatabase::TimeDB->CURRENTTIME;

    Output::print<1>("\nCURRENT TIME: ",
         TDatabase::TimeDB->CURRENTTIME);
    
    tcd.assemble();
    tcd.solve();
    
    if(checkOnlyFinalErrors)
      tcd.output();
    else
      testCN(tcd, step);
  }

}

void checkErrors(TimeConvectionDiffusion<2> &tcd,
                      ParameterDatabase& errorDatabase,
                      std::vector<std::string>& errorNames)
{
  std::array<double, 3> measuredErrors(tcd.get_errors());
  
  for(unsigned int i = 0; i < measuredErrors.size(); i++)
  {
    if( differentErrors(measuredErrors.at(i), errorDatabase[errorNames.at(i)]) )
      ErrThrow( "Test failed.\nError ",
                errorNames.at(i), " is too different.\n",
                "Its value is ",std::setprecision(14), measuredErrors.at(i),
                ", but it should be ", errorDatabase[errorNames.at(i)] );
  }
}

bool differentErrors(double measuredError, double expectedError)
{
  const double tolerance = 1.0e-7;
  
  if( std::abs(measuredError - expectedError) > tolerance )
    return true;
  else
    return false;
    
}

void testCN(TimeConvectionDiffusion<2> &tcd, int m)
{
  const int n_errors = 5;
  double errors[n_errors];
  errors[0]=errors[1]=errors[2]=errors[3]=errors[4]=0.;
  TAuxParam2D aux;
  MultiIndex2D AllDerivatives[3] = {MultiIndex2D::D00, MultiIndex2D::D10,
                                    MultiIndex2D::D01};
  const TFEFunction2D& function = tcd.get_function();
  const TFESpace2D* space = function.GetFESpace2D().get();

  function.GetErrors(tcd.get_example().get_exact(0), 3, AllDerivatives,
                     n_errors, conv_diff_l2_h1_linf_error<2>,
                     tcd.get_example().get_coeffs(), &aux, 1, &space, errors);
  double eps1 = 1E-6;
  double eps2 = 1E-5;
  if(m==0)
  {
    if( std::abs(errors[0] - 0.0705721) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct. ", errors[0]);
    if( std::abs(errors[1] - 6.91068) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");

  }
  else if(m==1)
  {
    if( std::abs(errors[0] - 0.0707957) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.85906) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==2)
  {
    if( std::abs(errors[0] - 0.0711781) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.81221) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==3)
  {
    if( std::abs(errors[0] - 0.0713873) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.76706) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==18)
  {
    if( std::abs(errors[0] - 0.0791431) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.17222) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==19)
  {
    if( std::abs(errors[0] - 0.0798896) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.13685) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
  else if(m==20)
  {
    if( std::abs(errors[0] - 0.0799563) > eps1 )
      ErrThrow("test Crank-Nicolson: L2 norm not correct.");
    if( std::abs(errors[1] - 6.10199) > eps2 )
      ErrThrow("test Crank-Nicolson: H1 norm not correct.");
  }
}
