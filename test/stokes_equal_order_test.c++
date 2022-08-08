#include <cmath>
#include <Domain.h>
#include <Database.h>
#include <ParameterDatabase.h>
#include "NavierStokes.h"
#include <Multigrid.h>
#include <Chrono.h>
#include <algorithm>
#include "LocalAssembling.h"
#include "ParMooN.h"

double accuracy = 1e-6;
bool testall = false;

#ifdef __2D__
typedef NavierStokes<2> Stokes;
#else // 3D
typedef NavierStokes<3> Stokes;
#endif

void compare(const Stokes& stokes, std::array<double, int(5)> errors)
{
  auto computed_errors = stokes.get_errors();
  
  // check the L2-error of the velcoity
  if( std::abs(computed_errors[0]-errors[0]) > accuracy )
  {
    ErrThrow("L2 norm of velocity: computed ", computed_errors[0],
             "; expected ", errors[0]);
  }
  // check the L2-error of the divergence of the velocity
  if( std::abs(computed_errors[1] - errors[1]) > accuracy )
  {
    ErrThrow("L2 norm of divergence: computed ", computed_errors[1],
             "; expected ", errors[1]);
  }
  // check the H1-error of the velcoity
  if( std::abs(computed_errors[2] - errors[2]) > accuracy )
  {
    ErrThrow("H1 norm of velocity: computed ", computed_errors[2],
             "; expected ", errors[2]);
  }
  // check the L2-error of the pressure
  if( std::abs(computed_errors[3] - errors[3]) > accuracy)
  {
    ErrThrow("L2 norm of pressure: computed ", computed_errors[3],
             "; expected ", errors[3]);
  }
  // check the H1-error of the pressure
  if(std::abs(computed_errors[4] - errors[4]) > accuracy )
  {
    ErrThrow("H1 norm of pressure: computed ", computed_errors[4],
             "; expected ", errors[4]);
  }
}

void compute(const TDomain& domain, ParameterDatabase& db,
             std::array<double, int(5)> errors)
{
  Stokes stokes(domain, db);
  stokes.assemble_linear_terms();
  stokes.solve();
  stokes.output();
  // compare now the errors
  compare(stokes, errors);
}

void check(const TDomain& domain, ParameterDatabase db,
           int velocity_order, int laplace_type,
           std::array<double, int(5)> errors)
{
  Output::print("\n\nCalling check with velocity_order=", velocity_order,
                ", laplace_type=", laplace_type);
  db.merge(Solver<>::default_solver_database());
  // db["problem_type"] = 5;
  db["solver_type"] = "direct";
  db["iterative_solver_type"] = "fgmres";
  db["residual_tolerance"] = 1.e-12;
  db["preconditioner"] = "sor";
  db["max_n_iterations"] = 2000;
  
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->LAPLACETYPE = laplace_type;
  db["nse_nonlinear_form"] = "convective";
  
  Chrono timer;
  compute(domain, db, errors);
  timer.restart_and_print("stokes direct solver,                    velocity "
                          + std::to_string(velocity_order));
  
#ifdef __3D__
  // too slow, needs more testing (e.g. which smoother is useful)
  if(velocity_order > 2)
    return;
#endif
  
  // we have to reset the space codes because they are changed in nse2d/nse3d
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  
  db["solver_type"] = "iterative";
  compute(domain, db, errors);
  timer.restart_and_print("stokes fgmres(lsc preconditioner),       velocity "
                          + std::to_string(velocity_order));
  
#ifdef __3D__
  // too slow with multigrid, needs more testing (e.g. which smoother is useful)
  if(velocity_order > 1)
    return;
#endif
  
  // we have to reset the space codes because they are changed in nse2d/nse3d
  TDatabase::ParamDB->PRESSURE_SPACE = velocity_order;
  TDatabase::ParamDB->VELOCITY_SPACE = velocity_order;
  
  db["preconditioner"] = "multigrid";
  //choose smoother on fine grid according to element
  std::vector<int> disc_p = {12,13,14,15,22,23,24};
  if(std::find(disc_p.begin(), disc_p.end(), velocity_order) != disc_p.end())
    db["multigrid_smoother"] = "cell_vanka_store";
  else
    db["multigrid_smoother"] = "patch_vanka_store";

  db["multigrid_type"] = "standard";
  db["multigrid_smoother_coarse"] = "direct_solve";
  db["multigrid_n_pre_smooth"] = 0;
  db["multigrid_n_post_smooth"] = 1;
  db["multigrid_correction_damp_factor"] = 1.0;
  db["multigrid_vanka_damp_factor"] = 1.0;
  compute(domain, db, errors);
  timer.restart_and_print("stokes fgmres(multigrid preconditioner), velocity "
                          + std::to_string(velocity_order));
}

#ifdef __2D__
void pspg_on_triangles(ParameterDatabase db)
{
  db["geo_file"] = "TwoTriangles";
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the P1/P1 elements on triangles");
  errors = {{ 0.14425380968049, 0.89062163325473, 1.5272900521314,
              1.3353534209894, 8.0717803059975 }};
  check(domain, db, 1, 0, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the P2/P2 elements on triangles");
  errors = {{ 0.0096576356355381, 0.092510617191165, 0.16390635960963,
              0.10816420781699, 1.9898981980016 }};
  check(domain, db, 2, 0, errors);

  if(testall)
  {
    //=========================================================================
    Output::print<1>("\nTesting the P3/P3 elements on triangles");
    errors = {{ 0.0013655488005913, 0.012260506129698, 0.017792417110209,
                0.022236536026839, 0.3674535115527 }};
    check(domain, db, 3, 0, errors);

    //=========================================================================
    Output::print<1>("\nTesting the P4/P4 elements on triangles");
    errors = {{ 5.7237636133624e-05, 0.0015288367247389, 0.0016535091698868,
                0.0017097099552422, 0.072004774733194 }};
    check(domain, db, 4, 0, errors);
  }
}

void pspg_on_quads(ParameterDatabase db)
{
  db["geo_file"] = "UnitSquare";
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // refine grid
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the Q1/Q1 elements on quads");
  errors = {{ 0.14770518128781, 0.5612829754172, 1.1794990495543,
              1.3441843485325, 8.0847736329225 }};
  check(domain, db, 1, 0, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the Q2/Q2 elements on quads");
  errors = {{ 0.0067110682559855, 0.065269023656333, 0.11630741869309,
              0.092921606881003, 1.5925769236787 }};
  check(domain, db, 2, 0, errors);

  if(testall)
  {
    //=========================================================================
    Output::print<1>("\nTesting the Q3/Q3 elements on quads");
    errors = {{ 0.0012077200634523, 0.011283264757702, 0.015848841146417,
                0.023611185292096, 0.33408390051687 }};
    check(domain, db, 3, 0, errors);

    //=========================================================================
    Output::print<1>("\nTesting the Q4/Q4 elements on quads");
    errors = {{ 2.9945160214221e-05, 0.00041913940235019, 0.00061350090960876,
                0.00083988300949934, 0.023863686042054 }};
    check(domain, db, 4, 0, errors);
  }
}

void check_other_stabilizations(ParameterDatabase db)
{
  db["geo_file"] = "TwoTriangles";
  db["example"] = 5;
  db["solver_type"] = "direct";
  TDatabase::ParamDB->VELOCITY_SPACE = 1;
  TDatabase::ParamDB->PRESSURE_SPACE = 1;
  check_parameters_consistency_NSE(db);
  // default construct a domain object
  TDomain domain(db);
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors{{0.89286842997152, 2.9279300657466,
                                     11.015597277738, 3.3072741836106,
                                     24.547961810913}};
  
  db["space_discretization_type"] = "local_projection";
  compute(domain, db, errors);
}
#endif // 2D

#ifdef __3D__
void pspg_on_tetrahedra(ParameterDatabase db)
{
  db["geo_file"] = "Default_UnitCube_Tetra";
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // unfortunately the domain refinement depends on whether multigrid is used 
  // or not. So the corresponding values are set here
  db["solver_type"] = "iterative";
  db["preconditioner"] = "multigrid";
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the P1/P1 elements on tetrahedra");
  errors = {{ 0.29285842764476, 1.7431669708558, 2.9734925486895,
              1.8843195210647, 15.547490346215 }};
  check(domain, db, 1, 0, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the P2/P2 elements on tetrahedra");
  errors = {{ 0.08968982096761, 0.66817478276149, 0.94421281958215,
              1.0325260352197, 9.6560759201588 }};
  check(domain, db, 2, 0, errors);

  //=========================================================================
//   Output::print<1>("\nTesting the P3/P3 elements on tetrahedra");
//   errors = {{  }};
//   check(domain, db, 3, 0, errors);
// 
  // no P4 elements on tetrahedra implemented
}

void pspg_on_hexahedra(ParameterDatabase db)
{
  db["geo_file"] = "Default_UnitCube_Hexa";
  // possibly parameters in the database
  check_parameters_consistency_NSE(db);
  
  // default construct a domain object
  TDomain domain(db);
  // unfortunately the domain refinement depends on whether multigrid is used 
  // or not. So the corresponding values are set here
  db["solver_type"] = "iterative";
  db["preconditioner"] = "multigrid";
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors;
  
  //=========================================================================
  Output::print<1>("\nTesting the Q1/Q1 elements on hexahedra");
  errors = {{ 0.22856023489804, 0.89092443247775, 1.6452219971921,
              1.973590165051, 15.739699775035 }};
  check(domain, db, 1, 0, errors);
  
  //=========================================================================
  Output::print<1>("\nTesting the Q2/Q2 elements on hexahedra");
  errors = {{ 0.0069455335607669, 0.082635103449365, 0.13147215146863,
              0.14765503705692, 2.3797784866446 }};
  check(domain, db, 2, 0, errors);

  if(testall)
  {
    //=========================================================================
    Output::print<1>("\nTesting the Q3/Q3 elements on hexahedra");
    errors = {{ 0.0022527831353322, 0.027168669298461, 0.030403766396944,
                0.035154057676842, 0.82782568660416 }};
    check(domain, db, 3, 0, errors);

    //=========================================================================
    Output::print<1>("\nTesting the Q4/Q4 elements on hexahedra");
    errors = {{ 4.8691604234513e-05, 0.0007091120130668, 0.0011725745728817,
                0.0020906154874585, 0.044064601854262 }};
    check(domain, db, 4, 0, errors);
  }
}

void check_other_stabilizations(ParameterDatabase db)
{
  db["geo_file"] = "Default_UnitCube_Hexa";
  db["example"] = 0; // linear velocity, constant pressure example
  db["solver_type"] = "direct";
  TDatabase::ParamDB->VELOCITY_SPACE = 1;
  TDatabase::ParamDB->PRESSURE_SPACE = 1;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  db["nse_nonlinear_form"] = "convective";
  check_parameters_consistency_NSE(db);
  // default construct a domain object
  TDomain domain(db);
  // refinement
  domain.refine_and_get_hierarchy_of_collections(db);
  std::array<double, int(5)> errors{{0., 0., 0.,0., 0.}};
  
  // note: for this example f=0 and using Q1/Q1 the laplacians vanish so that
  // each of the methods used below should give the same solution. Furthermore 
  // that solution is in the ansatz space and the computed errors are therefore
  // (close to) zero.
  db["space_discretization_type"] = "symm_gls";
  accuracy = 1.e-12;
  Chrono timer;
  compute(domain, db, errors);
  timer.restart_and_print("stokes direct solver, symmetric GLS");
  db["space_discretization_type"] = "nonsymm_gls";
  compute(domain, db, errors);
  timer.restart_and_print("stokes direct solver, non-symmetric GLS");
  db["space_discretization_type"] = "brezzi_pitkaeranta";
  compute(domain, db, errors);
  timer.restart_and_print("stokes direct solver, Brezzi-Pitkaeranta");
}

#endif // 3D

// =======================================================================
// main program
// =======================================================================
int main(int, char* argv[])
{
  ParameterDatabase db = parmoon::parmoon_initialize();
  if(argv[1])
    testall = (std::string(argv[1]).compare("testall") == 0);
  db.merge(ParameterDatabase::default_nonlinit_database());
  db.merge(ParameterDatabase::default_output_database());
  db.merge(Stokes::default_nse_database());
#ifdef __2D__
  db.merge(Example2D::default_example_database(), true);
#else // 3D
  db.merge(Example3D::default_example_database(), true);
#endif // 2D
  db.merge(TDomain::default_domain_parameters());
  db.merge(Multigrid::default_multigrid_database());

  db["problem_type"].set<size_t>(3); // Stokes
  db["example"] = 2;
  db["reynolds_number"] = 1;
  db["refinement_n_initial_steps"] = 2;
  db["multigrid_n_levels"] = db["refinement_n_initial_steps"].get<size_t>();
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE=3;
  db["space_discretization_type"] = "pspg";
  db["laplace_type_deformation"] = false;
  TDatabase::ParamDB->NSTYPE = 14;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  db["boundary_file"] = "Default_UnitSquare";
  
  
#ifdef __2D__
  pspg_on_triangles(db);
  pspg_on_quads(db);
#else // 3D
  db["boundary_file"] = "Default_UnitCube";
  pspg_on_tetrahedra(db);
  pspg_on_hexahedra(db);
#endif
  check_other_stabilizations(db);
  parmoon::parmoon_finalize();
}
