/** @brief A test program for solving NSE problems with DG methods in 2D.
 *
 * This serves as a test for solving NSE problems with discontinuous Galerkin
 * methods in 2D. It uses H(div)-conforming methods, namely Raviart-Thomas and
 * Brezzi-Douglas-Marini elements. It is intended to perform NSE calculations
 * with different examples in different setups to test the DG implementation.
 *
 * First, consistency is checked in the sense, that if the exact solution is
 * contained in the ansatz space, then the errors should be zero. Afterwards,
 * the norms of some solutions to different problems are compared with reference
 * norms. If those norms are not approximated well enough (or something in the
 * process goes wrong) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake. Or
 * you changed some program setup (e.g. changed the default solver). Then this
 * tests shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea to
 * describe your changes in the forum and request support.
 *
 * See also: cd2d_test.c++
 *
 *
 * @date 2021/10
 * @author Derk Frerichs-Mihov, Cristina Melnic
 */

#include "Database.h"
#include "NavierStokes_DG.h"
#include "ParMooN.h"
#include <ParMooN_repository_info.h>
#ifdef __3D__
ErrThrow("Neither the test nor the assemble routine is implemented yet in 3D.");
#else
constexpr int d = 2;
using NavierStokesND = NavierStokes_DG<2>;
#endif

#include <cmath>

/* ########################################################################## */
/* Declare functions and constant expressions, see below for definition */
/* ########################################################################## */
void set_solver(ParameterDatabase& db, int type);
void test_all(const bool& testall, ParameterDatabase parmoon_db);
void test_consistency(const bool& testall, ParameterDatabase& parmoon_db);
void check_nse(const TDomain& domain, ParameterDatabase& parmoon_db, const int&
    element_code, const std::array<double, 5>& errors);
void compareErrors(const NavierStokesND& nse, const std::array<double, 5>&
    errors);

/* ########################################################################## */
/* Declare some often used variables */
/* ########################################################################## */
const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";


/* ########################################################################## */
/* Main program */
/* ########################################################################## */
int main(int, char* argv[])
{
  parmoon::parmoon_initialize();
  int my_rank = 0;
  if(my_rank==0)
  {
    Output::print("\n        START ",d,"D TEST");
    Output::print("------------------------------\n\n");
  }

  bool testall = false;
  if (argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }

  Output::setVerbosity(2);

  // Construct the ParMooN database and set some parameters
  ParameterDatabase parmoon_db =
    NavierStokesND::default_NavierStokes_DG_database(true);

  TDatabase::ParamDB->USE_ISOPARAMETRIC = 0; // not implemented yet
  // Stokes
  TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 3;
  // Navier-Stokes stationary
  //TDatabase::ParamDB->FLOW_PROBLEM_TYPE = 5; // implementation not finished yet
  parmoon_db["space_discretization_type"].set("dg", true);

  parmoon_db["output_compute_errors"] = true;

  int max_solver_types = (testall) ? 2 : 1;
  max_solver_types = 1;
  for (auto solver_type = 1; solver_type <= max_solver_types; ++solver_type)
  {
    Output::print("Testing solver ", solver_type);
    //set_solver(parmoon_db, solver_type);
    set_solver(parmoon_db, solver_type);

    // Perform all test and print information if everything worked.
    test_all(testall, parmoon_db);
  }

  if(my_rank==0)
  {
    Output::print("\n\n------------------------------");
    Output::print("        ", d, "D TEST PASSED\n");
  }
  parmoon::parmoon_finalize();
}

void set_solver(ParameterDatabase& db, int type = 1)
{
  switch(type)
  {
    case 1:
      db["solver_type"] = "direct";
      break;
    case 2:
      db["solver_type"] = "iterative";
      db["solver_verbosity"] = 2;
      db["max_n_iterations"] = 10000;
      db["iterative_solver_type"] = "fgmres";
      db["residual_tolerance"] = 1.e-12;
      db["preconditioner"] = "least_squares_commutator";
      break;
    default:
      ErrThrow("call 'set_solver' with type 1 or 2. ", type);
      break;
  }
}

/* ########################################################################## */
/* Declaration of functions */
/* ########################################################################## */
// Perform the tests
void test_all(const bool& testall, ParameterDatabase parmoon_db)
{
  // Consistency check
  test_consistency(testall, parmoon_db);
}

void test_consistency(const bool& testall, ParameterDatabase& parmoon_db)
{
  // This test checks the consistency of the method, i.e. if the exact solution
  // is contained in the DG space, then the method should compute the exact
  // solution, and hence the error should be zero.
  int my_rank = 0;

  parmoon_db["example"] = 20; // Polynomial solution of a certain degree

  // If the solution is calculated exactly all errors should be zero
  // L2(u), L2(div(u)), H1-semi(u), L2(p), H1-semi(p)
  std::array<double, 5> errors = { 0, 0, 0, 0, 0 };

  // Define diffusion coefficient
  TDatabase::ParamDB->RE_NR = 10; // viscosity = 1/RE_NR

  std::array<std::string, 4> geometries;
  std::array<std::string, 4> boundary_files;
  if (d == 2)
  {
    // define geometries and initial refinements
    geometries = {"UnitSquare_quads.mesh", "UnitSquareCrissCross.GEO",
      "UnitSquareIrregular.GEO", "UnitSquareWithSphericalInscribedRegion.mesh"};
    boundary_files.fill("UnitSquare.PRM");
  }
  else
  {
  }

  // Test on all geometries
  int n_geoms_to_test = (!testall) ? 2 : (d == 2) ? geometries.size() : 0;
  for (int geom_i = 0; geom_i < n_geoms_to_test; ++geom_i)
  {
    // Define geometry
    std::string geo_file;
    std::string boundary_file;
    if (d == 2)
    {
      geo_file = path_to_meshes + geometries[geom_i];
      boundary_file = path_to_meshes + boundary_files[geom_i];
    }
    else
    {
    }
    parmoon_db["geo_file"].set(geo_file, false);
    parmoon_db["boundary_file"].set(boundary_file, false);

    check_parameters_consistency_NSE(parmoon_db);

    auto max_ref = (testall) ? 5 : 2;
    for (auto ref_i = 1; ref_i < max_ref; ++ref_i)
    {
      // refine at least once, to have interior edges
      parmoon_db["refinement_n_initial_steps"].set((size_t) ref_i);

      // Construct a domain and refine
      TDomain domain(parmoon_db);
      domain.refine_and_get_hierarchy_of_collections(parmoon_db);

      // Try symmetric, incomplete and non-symmetric IPG
      //for (int symmetry = 1; symmetry > -2; --symmetry)
      for (int symmetry = 1; symmetry > 0; --symmetry)
      {
      
      	
        parmoon_db["symmetry_DG"].set(symmetry, true);

        // Test elements with order up to some max_degree
        int max_degree = 1;
        if (testall)
        {
          if (d == 2)
          {
            // Actually higher orders are implemented but because of reasons that
            // Ulrich and Derk don't understand higher orders do not work
            // properly.
            // They believe that there are issues with the solver since the exact
            // solution solves the system but the solver is not able to compute
            // this solution.
            max_degree = (geom_i < 2) ? 2 : max_degree;
          }
          else
          {
          }
        }

        for (int degree_poly = 1; degree_poly <= max_degree; ++degree_poly)
        {
          parmoon_db["face_sigma_DG"].set(.5 * degree_poly*degree_poly *
              std::pow(10, symmetry), true);

          // Output
          if (my_rank == 0)
          {
            Output::print("\n\n---------- NEW TEST ----------\n\n"
                "Test with following properties:",
                "\nDimension:                ", d,
                "\nSolver:                   ", parmoon_db["solver_type"],
                "\nProblem:                  Polynomial of degree ", degree_poly,
                "\nDomain:                   ", parmoon_db["geo_file"],
                "\nViscosity coefficient     ", 1 / TDatabase::ParamDB->RE_NR,
                "\nsymmetry_DG:              ", parmoon_db["symmetry_DG"],
                "\nface_sigma_DG:            ", parmoon_db["face_sigma_DG"]
                );
          }

          // Compute solution and compare errors
          if (my_rank == 0)
          {
            Output::print("\nChecking RT elements of order ", degree_poly - 1);
          }
          // Check RTk
          parmoon_db["degree_polynomial"] = degree_poly - 1;
          parmoon_db["is_RT_polynomial"] = true;
          check_nse(domain, parmoon_db, 1000 + degree_poly - 1, errors);
          if (degree_poly > 3)
          {
            continue;
          }
          if (my_rank == 0)
          {
            Output::print("\nChecking BDM elements of order ", degree_poly);
          }
          // Check BDM
          parmoon_db["degree_polynomial"] = degree_poly;
          parmoon_db["is_RT_polynomial"] = false;
          check_nse(domain, parmoon_db, 1010 + degree_poly, errors); // BDM space

        } //endfor polynomial degree
      } //endfor symmetry of method
    } //endfor geom_i (different geometries)
  }
} // end test_consistency


void check_nse(const TDomain& domain, ParameterDatabase& parmoon_db, const int&
    element_code, const std::array<double, 5>& errors)
{ // Here the actual computations take place

  TDatabase::ParamDB->VELOCITY_SPACE = element_code;
  TDatabase::ParamDB->PRESSURE_SPACE = -4711;
  //TDatabase::ParamDB->NSTYPE = 3;
//  TDatabase::ParamDB->NSTYPE = 5;
  TDatabase::ParamDB->LAPLACETYPE = 0;
  parmoon_db["nse_nonlinear_form"] = "convective";

  NavierStokesND nse(domain, parmoon_db);
  if (nse.get_size() > 5e4)
  {
    // Only compute problems with at most 50000 dof
    return;
  }

  nse.assemble();
  //nse.stop_it(0);
  
  if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 3)
  {
  	nse.solve();
  }
  else if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5)
  {	
  	  Output::print("\n\n Warning: DG with H(div) for stationary Navier-Stokes not finished yet.\n\n");
  	  
  	  
  	  for(unsigned int k=1;k<10; k++)
	  {
	    //Output::print<1>("nonlinear step " , setw(3), k-1, "\t",
	      //               nse.get_residuals());
	    nse.solve();
	    // checking the first nonlinear iteration    
	    nse.assemble_nonlinear_term();;
	  //  if(nse.stop_it(k))
	    //  break;
	  }
  }
  
  //nse.solve();
  //nse.assemble_nonlinear_term();
  nse.output();
  
  // compare computed with given errors
  int my_rank = 0;
  if (my_rank == 0)
  {
    compareErrors(nse, errors); // throws upon a difference
  }

}

void compareErrors(const NavierStokesND& nse, const std::array<double, 5>&
    errors)
{ // compare the computed errors in the NSE object with the given ones in the
  // array
  const double eps = 3e-07;

  // check the errors
  {
    if( std::abs(nse.getL2VelocityError() - errors[0]) > eps )
    {
      ErrThrow("L2 velocity error not correct. Computed: ",
          nse.getL2VelocityError(), ", reference: ", errors[0]);
    }
    if( std::abs(nse.getL2DivergenceError() - errors[1]) > eps )
    {
      ErrThrow("L2 velocity divergence error not correct. Computed: ",
          nse.getL2DivergenceError(), ", reference: ", errors[1]);
    }
    if( std::abs(nse.getH1SemiVelocityError() - errors[2]) > eps )
    {
      ErrThrow("H1-semi velocity error not correct. Computed: ",
          nse.getH1SemiVelocityError(), ", reference: ", errors[2]);
    }
    if( std::abs(nse.getL2PressureError() - errors[3]) > eps )
    {
      ErrThrow("L2 pressure error not correct. Computed: ",
          nse.getL2PressureError(), ", reference: ", errors[3]);
    }
    if( std::abs(nse.getH1SemiPressureError() - errors[4]) > eps )
    {
      ErrThrow("H1-semi pressure error not correct. Computed: ",
          nse.getH1SemiPressureError(), ", reference: ", errors[3]);
    }
  }
}
