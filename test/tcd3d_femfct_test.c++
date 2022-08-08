/**
 * @brief A test program for a time dependent CD3D problem
 * algorithm "FEM-FCT", sequential and parallel versions
 *
 * @author Clemens Bartsch, Naveed Ahmed, Ondrej Partl
 */
#include "TimeConvectionDiffusion.h"
#include <Database.h>
#include <TimeDiscRout.h>
#include <TimeDiscretizations.h>
#include <MainUtilities.h>
#include <ConvDiff.h>
#include "AuxParam3D.h"
#include "ParMooN.h"
#include <cmath>

#include <ParMooN_repository_info.h>
const std::string path = parmoon::source_directory;
const std::string block_mesh_file = path + "/data/mesh/" + "block_1x2x3.mesh";

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
#endif

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
  ConcentrationSpecies,
  NonLinearWithDirichletNeumannEulerMesh,
  NonLinearWithDirichletNeumannCrankNicolsonMesh,
  MonolithicWithDirichlet,
  MonolithicWithDirichletNeumann,
  MolNonlinearWithDirichletCrankNicolson,
  MolNonLinearWithDirichletNeumannCrankNicolson,
  ZalExplicitWithDirichlet,
  ZalExplicitWithDirichletNeumann
};

void prepareErrorDatabase(ParameterDatabase& errorDatabase,
                          std::vector<std::string>& errorNames);
void setCommonParameters(ParameterDatabase& db);
void setTestSpecificParameters(TestNames testName,
                               ParameterDatabase& db,
                               ParameterDatabase& errorDatabase);
void timeIntegration(TimeConvectionDiffusion<3>& tcd,
                     ParameterDatabase& db,
                     ParameterDatabase errorDatabase);
bool differentErrors(double measuredError, double expectedError);
void checkErrors(TimeConvectionDiffusion<3> &tcd,
                 ParameterDatabase& errorDatabase,
                 std::vector<std::string>& errorNames);
void runTest(TestNames testName,
             ParameterDatabase& db,
             ParameterDatabase& errorDatabase,
             std::vector<std::string>& errorNames);
void check_solution_norms(TimeConvectionDiffusion<3> &tcd, int m);


int main()
{
  parmoon::parmoon_initialize();
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
#endif
  
  ParameterDatabase db = TimeConvectionDiffusion<3>::default_tcd_database(true);
  
  setCommonParameters(db);
  
  std::vector<std::string>
  errorNames = { "L2(0,T;L2(c))", "L2(0,T;H1(c))", "L_inf(0,T,L_inf(c))" };
  
  ParameterDatabase errorDatabase("database of errors");
  prepareErrorDatabase(errorDatabase, errorNames);
  
  runTest(ConcentrationSpecies, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletEuler, db, errorDatabase, errorNames);
  runTest(LinearWithDirichletNeumannEuler, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletEuler, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannCrankNicolsonMesh, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannEuler, db, errorDatabase, errorNames);
  runTest(NonLinearWithDirichletNeumannEulerMesh, db, errorDatabase, errorNames);
  runTest(MonolithicWithDirichlet, db, errorDatabase, errorNames);
  runTest(MonolithicWithDirichletNeumann, db, errorDatabase, errorNames);
  runTest(MolNonlinearWithDirichletCrankNicolson, db, errorDatabase, errorNames);
  runTest(MolNonLinearWithDirichletNeumannCrankNicolson, db, errorDatabase, errorNames);
  runTest(ZalExplicitWithDirichlet, db, errorDatabase, errorNames);
  runTest(ZalExplicitWithDirichletNeumann, db, errorDatabase, errorNames);
  
  parmoon::parmoon_finalize();
}

void setCommonParameters(ParameterDatabase& db)
{
  db["algebraic_flux_correction"].set("fem-fct-cn");
  db["space_discretization_type"] = "galerkin";
  TDatabase::ParamDB->ANSATZ_ORDER = 1;
}

void prepareErrorDatabase(ParameterDatabase& errorDatabase,
                          std::vector<std::string>& errorNames)
{
  errorDatabase.add("checkOnlyFinalErrors", true, "", {true, false});
  
    for( unsigned int i = 0; i < errorNames.size(); i++)
      errorDatabase.add(errorNames.at(i), 0.0, "", 0.0, 1.0e3);
}

void runTest(TestNames testName,
             ParameterDatabase& db,
             ParameterDatabase& errorDatabase,
             std::vector<std::string>& errorNames)
{
  #ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0)
  #endif 
    Output::print("\nStarting test number ", (int)testName, "\n");
  
  setTestSpecificParameters(testName, db, errorDatabase);
  
  TDatabase::TimeDB->TIMESTEPLENGTH = db["time_step_length"];
  TDatabase::TimeDB->CURRENTTIME = db["time_start"];
  
  TDomain domain(db);
  SetTimeDiscParameters(0);
  domain.refine_and_get_hierarchy_of_collections(db);
  TimeConvectionDiffusion<3> tcd(domain, db);
  
  timeIntegration(tcd, db, errorDatabase);
  
  if(errorDatabase["checkOnlyFinalErrors"])
  {
    #ifdef _MPI
    if(rank == 0)
    #endif 
      checkErrors(tcd, errorDatabase, errorNames);
  }
}

void setTestSpecificParameters(TestNames testName,
                               ParameterDatabase& db,
                               ParameterDatabase& errorDatabase)
{
  switch(testName)
  {
    case ConcentrationSpecies:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 0;
    db["solver_type"] = "iterative";
    db["iterative_solver_type"] = "fgmres";
    db["preconditioner"] = "jacobi";
    db["residual_tolerance"] = 1.0e-13;
    db["residual_reduction"] =  0.0;
    db["afc_prelimiter"].set(0);
    db["boundary_file"].set("Default_UnitCube");
    db["geo_file"].set("Default_UnitCube_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.01;
    errorDatabase["checkOnlyFinalErrors"] = false;
    break;
    
    case LinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.035914027505992;
    errorDatabase["L2(0,T;H1(c))"] = 0.17860811359288;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.085134085563176;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Tetra");
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.55398984224944;
    errorDatabase["L2(0,T;H1(c))"] = 2.379345453276;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.941369333252;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 8.8055574248847e-16;
    errorDatabase["L2(0,T;H1(c))"] = 6.5734837542128e-15;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 3.1086244689504e-15;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case LinearWithDirichletNeumannEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.10804978615762;
    errorDatabase["L2(0,T;H1(c))"] = 0.81313512105758;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.2528267239211;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(1);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 8.1927982446379e-16;
    errorDatabase["L2(0,T;H1(c))"] = 6.1801625955912e-15;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 2.6645352591004e-15;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.296944423616;
    errorDatabase["L2(0,T;H1(c))"] = 1.4710455398687;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.7573470201824;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannCrankNicolsonMesh:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set("");
    db["geo_file"].set(block_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.40278680027633;
    errorDatabase["L2(0,T;H1(c))"] = 2.7373489868535;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.8120390944888;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 8.7309879728458e-16;
    errorDatabase["L2(0,T;H1(c))"] = 6.51670407249e-15;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.2166812511657e-30;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannEuler:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Tetra");
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.63648491707003;
    errorDatabase["L2(0,T;H1(c))"] = 2.5904088669627;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.737748120519;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case NonLinearWithDirichletNeumannEulerMesh:
    db["diffusion_coefficient"] = 1e-6;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["afc_prelimiter"].set(2);
    db["boundary_file"].set("");
    db["geo_file"].set(block_mesh_file, false);
    db["refinement_n_initial_steps"].set(1);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("zalesak");
    db["solDotUsedInZalesak"] = false;
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "backward_euler";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.40929436572019;
    errorDatabase["L2(0,T;H1(c))"] = 2.6360317633841;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.9358730843547;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MonolithicWithDirichlet:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("monolithic");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.025;
    errorDatabase["L2(0,T;L2(c))"] = 0.1120377565396;
    errorDatabase["L2(0,T;H1(c))"] = 0.61275694556238;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.18740975662506;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MonolithicWithDirichletNeumann:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Tetra");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("monolithic");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.04;
    errorDatabase["L2(0,T;L2(c))"] = 0.20893863957259;
    errorDatabase["L2(0,T;H1(c))"] = 1.5405249101034;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.82946169315684;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MolNonlinearWithDirichletCrankNicolson:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("monolithic");
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.05;
    errorDatabase["L2(0,T;L2(c))"] = 0.11225426870716;
    errorDatabase["L2(0,T;H1(c))"] = 0.6142481896612;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.18765027365623;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case MolNonLinearWithDirichletNeumannCrankNicolson:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("non-linear");
    db["afc_limiter"].set("monolithic");
    db["afc_nonlinloop_maxit"].set(1000);
    db["time_discretization"] = "crank_nicolson";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.1;
    errorDatabase["L2(0,T;L2(c))"] = 0.15353636183791;
    errorDatabase["L2(0,T;H1(c))"] = 0.63906771308993;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.51956544729604;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case ZalExplicitWithDirichlet:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 2;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Hexa");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("zalesak");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.025;
    errorDatabase["L2(0,T;L2(c))"] = 0.017220093090019;
    errorDatabase["L2(0,T;H1(c))"] = 0.083065888822555;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 0.037909515114179;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    case ZalExplicitWithDirichletNeumann:
    db["diffusion_coefficient"] = 1e-1;
    db["example"] = 3;
    db["solver_type"] = "direct";
    db["boundary_file"].set("Default_Block_1x2x3");
    db["geo_file"].set("Default_Block_1x2x3_Tetra");
    db["refinement_n_initial_steps"].set(2);
    db["afc_fct_scheme"].set("explicit");
    db["afc_limiter"].set("zalesak");
    db["time_discretization"] = "Runge_Kutta_Heun";
    db["time_start"] = 0.0;
    db["time_end"] = 1.0;
    db["time_step_length"] = 0.04;
    errorDatabase["L2(0,T;L2(c))"] = 0.20284977632633;
    errorDatabase["L2(0,T;H1(c))"] = 1.4279145289208;
    errorDatabase["L_inf(0,T,L_inf(c))"] = 1.103462904841;
    errorDatabase["checkOnlyFinalErrors"] = true;
    break;
    
    default:
      ErrThrow("Unspecified test");
  }
}


void timeIntegration(TimeConvectionDiffusion<3>& tcd,
                     ParameterDatabase& db,
                     ParameterDatabase errorDatabase)
{
  TimeDiscretization& tss = tcd.get_time_stepping_scheme();
  tss.current_step_ = 0;
  tss.current_time_ = db["time_start"];
  
  int rank = 0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  tcd.assemble_initial_time();
  
  bool checkOnlyFinalErrors = errorDatabase["checkOnlyFinalErrors"];
  
  if(checkOnlyFinalErrors)
    tcd.output();

  double end_time = tcd.get_db()["time_end"];
  while(TDatabase::TimeDB->CURRENTTIME < end_time -1e-10)
  {
//     tss.nonlin_iteration = 0; // partl, not needed, set in for loop
    tss.current_step_++;
    tss.set_time_disc_parameters();
    SetTimeDiscParameters(1);
    
    tss.current_time_ += tss.get_step_length();
    TDatabase::TimeDB->CURRENTTIME += tss.get_step_length();

    if(rank == 0)
      Output::print<1>("CURRENT TIME: ", tss.current_time_);
    
    tcd.assemble();
    tcd.solve();
    
    if(checkOnlyFinalErrors)
      tcd.output();
    else
      check_solution_norms(tcd, tss.current_step_);
  }

}

void checkErrors(TimeConvectionDiffusion<3> &tcd,
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
// these are the previous errors - you should get them if
// you do not change the values at node 59 by prescribing
// the Dirichlet boundary conditions - comment by Ondrej Partl
// std::vector<std::vector<double>> target_norms=
// {
//   {0.007686841218, 0.09028123769, 4.647793003e-09}, //t=0.1
//   {0.01597853147, 0.179622979, 2.226841822e-07}, //t=0.2
//   {0.02482957109, 0.2680028932, 1.593788446e-06}, //t=0.3
//   {0.03394574615, 0.3538330922, 5.448744344e-06}, //t=0.4
//   {0.04302352518, 0.4349243671, 1.281619336e-05}, //t=0.5
//   {0.05175860945, 0.5089748268, 2.462124805e-05}, //t=0.6
//   {0.05986967233, 0.5738383662, 4.054345229e-05}, //t=0.7
//   {0.06710294614, 0.6276691113, 5.919686464e-05},  //t=0.8
//   {0.07323577128, 0.6689863334, 7.873596671e-05},  //t=0.9
//   {0.07808004918, 0.6967105354, 9.708502041e-05}   //t=1
// };
std::vector<std::vector<double>> target_norms=
{
  {0.007686841351, 0.09028123756, 0.1098068789}, //t=0.1
  {0.01597853364, 0.1796229804, 0.2190507937}, //t=0.2
  {0.02482957961, 0.2680028848, 0.324693393}, //t=0.3
  {0.0339457792, 0.3538330112, 0.4234723442}, //t=0.4
  {0.04302370833, 0.4349238981, 0.5125560048}, //t=0.5
  {0.0517596078, 0.5089738488, 0.5894766554}, //t=0.6
  {0.0598735036, 0.5738399255, 0.6521594774}, //t=0.7
  {0.06711401485, 0.6276851275, 0.6989447884},  //t=0.8
  {0.07326161907, 0.6690452317, 0.7286082832},  //t=0.9
  {0.07813175239, 0.6968646354, 0.7403766544}   //t=1
};

void check_solution_norms(TimeConvectionDiffusion<3> &tcd, int m)
{

  if(m%10 == 0)
  {
    MultiIndex3D allDerivatives[4] = { MultiIndex3D::D000, MultiIndex3D::D100,
                                       MultiIndex3D::D010, MultiIndex3D::D001 };
    TAuxParam3D aux(1, 0, 0, 0, nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr);
    const int n_errors = 5;
    std::array<double, n_errors> locError = {};
    const TFESpace3D* space = tcd.get_function().GetFESpace3D().get();

    tcd.get_function().GetErrors(tcd.get_example().get_exact(0),
                                 4, allDerivatives,
                                 n_errors, conv_diff_l2_h1_linf_error<3>,
                                 tcd.get_example().get_coeffs(), &aux, 1, &space,
                                 locError.data());

  #ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> errorsReduced(n_errors);
    MPI_Reduce(locError.data(), errorsReduced.data(), n_errors-1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&locError[n_errors-1], &errorsReduced[n_errors-1], 1, MPI_DOUBLE,
               MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank == 0)
  {//this is only performed on root - just root willl have the correct values
    for(int j = 0; j < n_errors-1; j++)
      locError[j] = std::sqrt(errorsReduced[j]);
    locError[n_errors-1] = errorsReduced[n_errors-1];
  #endif // _MPI
    
    Output::print("Norms in step ", m);
    Output::dash(std::setprecision(10), locError[0], ", ", locError[1], ", ", locError[n_errors-1]);

    int control_step = m/10 - 1;
    double tol = 1e-7;

    cout<<TDatabase::TimeDB->CURRENTTIME<<endl;
    if( std::abs(locError[0] - target_norms[control_step][0]) > tol )
      ErrThrow("L2 norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ", "calculated norm: ",
               locError[0], ",  reference norm: ", target_norms[control_step][0]);

    if( std::abs(locError[1]  - target_norms[control_step][1]) > tol )
      ErrThrow("H1 norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ", "calculated norm: ",
               locError[1], ",  reference norm: ", target_norms[control_step][1]);

    if( std::abs(locError[n_errors-1]  - target_norms[control_step][2]) > tol )
      ErrThrow("L^inf norm at timestep ", TDatabase::TimeDB->CURRENTTIME,  " is not correct, ", "calculated norm: ",
               locError[n_errors-1], ",  reference norm: ", target_norms[control_step][2]);
    Output::print("Solution norms checked succesfully.");
 #ifdef _MPI
  }//end if rank == 0
#endif
  }
}
