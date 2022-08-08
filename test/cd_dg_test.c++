/**
 * @brief A test program for the solving of CD problems with DG methods in 2D
 * and 3D.
 *
 * This serves as a test for the solving of CD problems with discontinuous
 * Galerkin methods in 2D and 3D. It is intended to perform CD calculations with
 * different examples in different setups to test the DG implementation.
 *
 * First, consistency is checked in the sense, that if the exact solution is
 * contained in the DG space, then the errors should be zero. Afterwards, the
 * norms of some solutions to different problems are compared with reference
 * norms.  If those norms are not approximated well enough (or something in the
 * process goes wrong) the test fails.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.  Or
 * you changed some program setup (e.g. changed the default solver). Then this
 * tests  shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea to
 * describe your changes in the forum and request support.
 *
 * See also: cd2d_test.c++
 *
 *
 * @date 2020/11 and 2021/01
 * @author Derk Frerichs
 */

#include "Domain.h"
#include "Database.h"
#include "ConvectionDiffusion.h"
#include "Assemble_DG.h"
#include "MooNMD_Io.h"
#include <ParMooN_repository_info.h>
#include "ParMooN.h"
#include <cmath>
#include <iomanip>
#include <string>
#ifdef __3D__
constexpr int d = 3;
using ConvectionDiffusionND = ConvectionDiffusion<3>;
#else
constexpr int d = 2;
using ConvectionDiffusionND = ConvectionDiffusion<2>;
#endif
#include <Chrono.h>

#ifdef _MPI
#include <mpi.h>
#include <MeshPartition.h>
#endif

/* ########################################################################## */
/* Declare functions and constant expressions, see below for definition */
/* ########################################################################## */
void test_all(const bool& testall, ParameterDatabase parmoon_db);
void compareErrors(const ConvectionDiffusionND& cd,
    const std::array<double, 5>& errors);
void check_cd(const TDomain& domain, ParameterDatabase& parmoon_db,
    const int& element_code, const std::array<double, 5>& errors);
void test_consistency( const bool& testall, ParameterDatabase& parmoon_db);
void test_reference_setting(const bool& testall, ParameterDatabase& parmoon_db,
    const int& problem_number);
void get_setting_parameter(const int& problem_number, ParameterDatabase&
    parmoon_db, int symmetries [3],  double face_sigmas [3]);
void get_reference_values(const int& problem_number, const int& symmetry, const
    unsigned int& refinement_level, const unsigned int& degree,
    std::array<double, 5>& errors);

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
#ifdef _MPI
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
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

  ParameterDatabase parmoon_db =
    ConvectionDiffusionND::default_cd_database(true);
  parmoon_db["solver_type"] = "direct";
  TDatabase::ParamDB->USE_ISOPARAMETRIC = 0;
  parmoon_db["space_discretization_type"].set("dg", true);
  parmoon_db["eta_upwind_DG"].set(1, true);

  parmoon_db["residual_tolerance"] = 1.0e-13;
  parmoon_db["output_compute_errors"] = true;

  // Perform all test and print information if everything worked.
  test_all(testall, parmoon_db);

  if(my_rank==0)
  {
    Output::print("\n\n------------------------------");
    Output::print("        ",d,"D TEST PASSED\n");
  }
  parmoon::parmoon_finalize();
}

/* ########################################################################## */
/* Declaration of functions */
/* ########################################################################## */
void test_all(const bool& testall, ParameterDatabase parmoon_db)
{
  // Perform the tests

  // Consistency check
  test_consistency(testall, parmoon_db);

  // Check reference values for some configurations
  std::vector<int> problem_numbers; // list of tested geometries
  if (d == 2)
  {
    if (testall)
    {
      problem_numbers = {0, 6, 7, 11};
    }
    else
    {
      problem_numbers = {6, 11};
    }
  }
  else
  {
    if (testall)
    {
      problem_numbers = {0, 4, -4, 6};
    }
    else
    {
      problem_numbers = {6, -4};
    }
  }
  for (auto problem_i : problem_numbers)
  { // Test several reference settings
    test_reference_setting(testall, parmoon_db, problem_i);
  }
}

void test_consistency( const bool& testall, ParameterDatabase& parmoon_db)
{ // This test checks the consistency of the method, i.e. if the exact solution
  // is contained in the DG space, then the method should compute the exact
  // solution, and hence the error should be zero.
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  // polynomial example should be calculated exactly
  if (d == 2)
  {
    parmoon_db["example"] = 10; // Polynomial solution of a certain degree
  }
  else
  {
    parmoon_db["example"] = -4; // Polynomial solution of a certain degree
  }

  // If the solution is calculated exactly all errors should be zero
  std::array<double, 5> errors = { 0, 0, 0, 0, 0 }; // L2, H1semi, SD, DG, L_inf

  // Define diffusion coefficient
  parmoon_db["diffusion_coefficient"] = 0.1;

  std::array<std::string, 5> geometries;
  std::array<std::string, 5> boundary_files;
  if (d == 2)
  {
    // define geometries and initial refinements
    geometries = {"UnitSquare_quad.mesh",
      "UnitSquareCrissCross.GEO", "UnitSquareIrregular.GEO",
      "UnitSquare_quads.mesh", "UnitSquareWithSphericalInscribedRegion.mesh"};
    boundary_files = {"UnitSquare.PRM", "UnitSquare.PRM",
      "UnitSquare.PRM", "UnitSquare.PRM", "UnitSquare.PRM"};
  }
  else
  {
    // define geometries and initial refinements
    geometries = {"Default_UnitCube_Tetra", "Default_UnitCube_Hexa", path_to_meshes + "UnitCube_TwoHexahedra.mesh"};
    boundary_files = {"Default_UnitCube", "Default_UnitCube", ""};
  }

  // refine at least once, to have interior edges
  size_t n_init_refinements = 1;
  parmoon_db["refinement_n_initial_steps"].set((size_t) n_init_refinements);

  // Test on all geometries
  int n_geoms_to_test = (!testall) ? 2 : (d == 2) ? geometries.size() : 3;
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
      geo_file = geometries[geom_i];
      boundary_file = boundary_files[geom_i];
    }
    parmoon_db["geo_file"].set(geo_file, false);
    parmoon_db["boundary_file"].set(boundary_file, false);

    // Construct a domain and refine
    TDomain domain(parmoon_db);
    domain.refine_and_get_hierarchy_of_collections(parmoon_db);

    // Try symmetric, incomplete and non-symmetric IPG
    for (int symmetry = 1; symmetry > -2; --symmetry)
    {
      parmoon_db["symmetry_DG"].set(symmetry, true);

      // Test polynomials up to degrees that are implemented
      int max_degree = 1;
      if (testall)
      {
        if (d == 2)
        {
          max_degree = 4;
        }
        else if (d == 3 && geom_i > 0)
        {
          // In 3D only Q1 and Q2 are implemented up to now
          max_degree = 2;
        }
        else
        {
          // In 3D only P1, P2 and P3 are implemented up to now
          max_degree = 3;
        }
      }

      for (int degree_poly = 1; degree_poly <= max_degree; ++degree_poly)
      {
        parmoon_db["face_sigma_DG"].set(0.5 * degree_poly*degree_poly *
            std::pow(10, symmetry-1), true);
        parmoon_db["degree_polynomial"] = degree_poly;

        // Output
        if (my_rank == 0)
        {
        Output::print("\n\n---------- NEW TEST ----------\n\n"
            " Test with following properties:",
            "\nDimension:                ", d,
            "\nSolver:                   ", parmoon_db["solver_type"],
            "\nProblem:                  Polynomial of degree ", degree_poly,
            "\nDomain:                   ", parmoon_db["geo_file"],
            "\nDiffusion coefficient     ", parmoon_db["diffusion_coefficient"],
            "\nsymmetry_DG:              ", parmoon_db["symmetry_DG"],
            "\nface_sigma_DG:            ", parmoon_db["face_sigma_DG"]
            );
        Output::print("\nstarting with discontinuous elements of degree ",
            degree_poly);
        }

        // Compute solution and compare errors
        check_cd(domain, parmoon_db, -10-degree_poly, errors);

      } //endfor polynomial degree
    } //endfor symmetry of method
  } //endfor geom_i (different geometries)
} // end test_consistency


void test_reference_setting(const bool& testall, ParameterDatabase& parmoon_db,
    const int& problem_number)
{ // This test checks for some reference setting whether the correct errors and
  // convergence rates are computed.
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  parmoon_db["example"] = problem_number;

  // Get details of current setting, i.e. domain, diffusion coefficient,
  // symmetry and face_sigma
  int symmetries [3];
  double face_sigmas [3];
  get_setting_parameter(problem_number, parmoon_db, symmetries, face_sigmas);

  // Output
  std::string problem_name;
  double diffusion_coefficient = parmoon_db["diffusion_coefficient"];

  if (d == 2)
  {
    switch (problem_number)
    {
      case 0:
        problem_name = "SineLaplace"; break;
      case 5:
        problem_name = "HMM86"; break;
      case 6:
        problem_name = "hump"; break;
      case 7:
        problem_name = "boundary_layer_known"; break;
      case 11:
        problem_name = "convection-reaction (pure convection-reaction)";
        break;
      default:
        ErrThrow("cd_dg_test.C : problem_number ", problem_number,
            " not known.");
    }
  }
  else
  {
    switch (problem_number)
    {
      case 0:
        problem_name = "SineLaplace"; break;
      case -4:
        parmoon_db["degree_polynomial"] = 6;
        problem_name = "Polynomial of degree 6";
        break;
      case 4:
        problem_name = "Boundary layer known 3D"; break;
      case 6:
        problem_name = "Pure Convection"; break;
      default:
        ErrThrow("cd_dg_test.C : problem_number ", problem_number,
            " not known.");
    }
  }

  // Compute solutions on a series of refined grids
  // unsigned int init_n_refs = (testall) ? 0 : 2;  // initial number refinements
  unsigned int init_n_refs = (testall) ? 1 : 2;  // initial number refinements
  unsigned int max_n_refs = (testall) ? (d == 2) ? 5 : 4 : 3;  // maximal number refinements
  // Test polynomials up to degrees that are implemented
  int max_degree = (d == 2) ? 4 : 2;
  if (testall && d == 3 && (problem_number == 0 || problem_number == 6))
  {
    // In 3D only Q1 and Q2 are implemented up to now
    max_degree = 3;
  }
  max_degree = (testall) ? max_degree : 1;

  for (unsigned int refine_lvl = init_n_refs; refine_lvl < max_n_refs;
      ++refine_lvl)
  { // refine grid

    for (int symmetry_i = 0; symmetry_i < 3; ++symmetry_i)
    { // test SIPG, IIPG and NIPG

      // Fix setting parameters that have not been fixed in
      parmoon_db["symmetry_DG"].set(symmetries[symmetry_i], true);
      parmoon_db["face_sigma_DG"].set(face_sigmas[symmetry_i], true);

      parmoon_db["refinement_n_initial_steps"].set((size_t) refine_lvl);
      // Output
      if (my_rank == 0)
      {
      Output::print("\n\n---------- NEW TEST ----------\n\nTest with following",
          " properties:",
          "\nDimension:                 ", d,
          "\nSolver:                    ", parmoon_db["solver_type"],
          "\nProblem:                   ", problem_name,
          "\nDomain:                    ", parmoon_db["geo_file"],
          "\nDiffusion coefficient      ", diffusion_coefficient,
          "\nDG_SYMMETRY:               ", parmoon_db["symmetry_DG"],
          "\nFACE_SIGMA:                ", parmoon_db["face_sigma_DG"],
          "\nRefinement level:          ", parmoon_db["refinement_n_initial_steps"]
          );
      }


      // Construct a domain, refine and get diameter of elements
      TDomain domain(parmoon_db);
      domain.refine_and_get_hierarchy_of_collections(parmoon_db);

      // Compute solutions for all FE-Spaces of order 1,2,...,max_degree
      for (int degree_space = 1; degree_space <= max_degree;
          ++degree_space)
      {

        // Output
        if (my_rank == 0)
        {
          Output::print("\nDiscontinuous elements of degree ", degree_space);
        }

        // Solve the problem and compare the errors
        std::array<double, 5> errors = { 0, 0, 0, 0, 0 };
        get_reference_values(problem_number, symmetries[symmetry_i],
            refine_lvl, degree_space, errors);
        check_cd(domain, parmoon_db, -10-degree_space, errors);

      } // endfor degree_space
      if ( (d == 2 && problem_number == 11) || (d == 3 && problem_number == 6))
      { // pure convection-reaction does not depend on symmetry. Hence, one
        // setting is enough
        break;
      }
    } //endfor SIPG, IIPG, NIPG
  } // endfor refinements
} // end test_reference_setting

void check_cd(const TDomain& domain, ParameterDatabase& parmoon_db,
    const int& element_code, const std::array<double, 5>& errors)
{ // Here the actual computations take place
  TDatabase::ParamDB->ANSATZ_ORDER = element_code;

  ConvectionDiffusionND cd(domain, parmoon_db);
  cd.assemble();
  cd.solve();
  cd.output();
  // compare computed with given errors
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  if (my_rank == 0)
  {
    compareErrors(cd, errors); // throws upon a difference
  }
}

void compareErrors(const ConvectionDiffusionND& cd,
    const std::array<double, 5>& errors)
{ // compare the computed errors in the CD object with the given ones in the
  // array
  const double eps = 1e-9;

  // check the errors
  {
    if( std::abs(cd.get_L2_error() - errors[0]) > eps )
    {
      ErrThrow("L2 error not correct. Computed: ",
          cd.get_L2_error(),", reference: ", errors[0]);
    }
    if( std::abs(cd.get_H1_semi_error() - errors[1]) > eps )
    {
      ErrThrow("H1-semi error not correct. Computed: ",
          cd.get_H1_semi_error(),", reference: ", errors[1]);
    }
    if( std::abs(cd.get_SD_error() - errors[2]) > eps )
    {
      ErrThrow("SD error not correct. Computed: ",
          cd.get_SD_error(),", reference: ", errors[2]);
    }
    if( std::abs(cd.get_DG_error() - errors[3]) > eps )
    {
      ErrThrow("DG error not correct. Computed: ",
          cd.get_DG_error(),", reference: ", errors[3]);
    }
    if( std::abs(cd.get_L_inf_error() - errors[4]) > eps )
    {
      ErrThrow("L_inf error not correct. Computed: ",
          cd.get_L_inf_error(),", reference: ", errors[3]);
    }
  }
}

void get_setting_parameter(const int& problem_number, ParameterDatabase& parmoon_db,
    int symmetry_values [3], double face_sigmas [3])
{ // Get the parameters of the current setting, i.e. domain, diffusion
  // coefficient, symmetry and face_sigma
  symmetry_values[0] = 1;
  symmetry_values[1] = 0;
  symmetry_values[2] = -1;
  std::string mesh = path_to_meshes;
  std::string boundary = path_to_meshes;
  double diffusion_coefficient;
  if (d == 2)
  {
    switch (problem_number)
    {
      case 0: // SineLaplace
        mesh += "UnitSquareIrregular.GEO";
        boundary += "UnitSquare.PRM";
        diffusion_coefficient= 1;
        face_sigmas[0] = 75;
        face_sigmas[1] = 5;
        face_sigmas[2] = .5;
        break;
      case 5:
        mesh += "UnitSquare_tria.mesh";
        boundary += "UnitSquare.PRM";
        diffusion_coefficient= 1e-8;
        face_sigmas[0] = 25;
        face_sigmas[1] = 5;
        face_sigmas[2] = .5;
        break;
      case 6:
        mesh += "UnitSquare_quads.mesh";
        boundary += "UnitSquare.PRM";
        diffusion_coefficient= 1;
        face_sigmas[0] = 25;
        face_sigmas[1] = 5;
        face_sigmas[2] = .5;
        break;
      case 7:
        mesh += "UnitSquareIrregular.GEO";
        boundary += "UnitSquare.PRM";
        diffusion_coefficient= 1;
        face_sigmas[0] = 45;
        face_sigmas[1] = 10;
        face_sigmas[2] = 1;
        break;
      case 11:
        mesh += "UnitSquare_quads.mesh";
        boundary += "UnitSquare.PRM";
        diffusion_coefficient= 0;
        face_sigmas[0] = 0;
        face_sigmas[1] = 0;
        face_sigmas[2] = 0;
        break;
      default:
        ErrThrow("cd_dg_test.C : get_setting_parameter() : problem_number ",
            problem_number, " not known.");
    }
  }
  else
  {
    switch (problem_number)
    {
      case 0: // SineLaplace
        mesh = "Default_UnitCube_Tetra";
        boundary = "Default_UnitCube";
        diffusion_coefficient= 1e-8;
        face_sigmas[0] = 75 * diffusion_coefficient;
        face_sigmas[1] = 5 * diffusion_coefficient;
        face_sigmas[2] = .5 * diffusion_coefficient;
        break;
      case -4: // Smooth polynomial of a degree
        mesh += "UnitCube_TwoHexahedra.mesh";
        boundary = "";
        diffusion_coefficient= 1e-1;
        face_sigmas[0] = 75 * diffusion_coefficient;
        face_sigmas[1] = 5 * diffusion_coefficient;
        face_sigmas[2] = .5 * diffusion_coefficient;
        break;
      case 4: // Boundary_layer_known_3D
        mesh += "UnitCube_TwoHexahedra.mesh";
        boundary = "";
        diffusion_coefficient = 1;
        face_sigmas[0] = 75 * diffusion_coefficient;
        face_sigmas[1] = 5 * diffusion_coefficient;
        face_sigmas[2] = .5 * diffusion_coefficient;
        break;
      case 6: // Pure_Convection
        mesh = "Default_UnitCube_Tetra";
        boundary = "Default_UnitCube";
        diffusion_coefficient= 0;
        face_sigmas[0] = 0;
        face_sigmas[1] = 0;
        face_sigmas[2] = 0;
        break;
      default:
        ErrThrow("cd_dg_test.C : get_setting_parameter() : problem_number ",
            problem_number, " not known.");
    }

  }
  parmoon_db["diffusion_coefficient"].set(diffusion_coefficient, false);
  parmoon_db["geo_file"].set(mesh, false);
  parmoon_db["boundary_file"].set(boundary, false);
}

/* ########################################################################## */
/* Reference values */
/* ########################################################################## */
void get_reference_values(const int& problem_number, const int& symmetry, const
    unsigned int& refinement_level, const unsigned int& degree,
    std::array<double, 5>& errors)
{ // Get reference values for all problems
  errors = {0,0,0,0,0};

  if (d == 2)
  {
    switch (problem_number)
    {
      case 0: // SineLaplace
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.12503701023662, 1.0580999415939, 1.0580999415939, 1.0740127694631, 0.28894076609974};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.052213994724722, 0.63189315447449, 0.63189315447449, 0.67728397079195, 0.1271492641473};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.002964024379374, 0.05362346219136, 0.05362346219136, 0.055358646763649, 0.0077875973695272};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.0014447736503394, 0.032268721580794, 0.032268721580794, 0.036808286866512, 0.0046148278801493};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.068089180617228, 0.74933709314943, 0.74933709314943, 0.76824427437795, 0.17568911145844};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0055409574655101, 0.14136407771959, 0.14136407771959, 0.15144603662601, 0.016784599115697};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00053637815924901, 0.016539954368888, 0.016539954368888, 0.01769326216064, 0.001600761395211};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.6709279509995e-05, 0.0016461001788205, 0.0016461001788205, 0.0018451336862926, 0.00012338204200671};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.019108221498528, 0.39589437178357, 0.39589437178357, 0.40584425097209, 0.060730495199741};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00073485132455766, 0.036790515550466, 0.036790515550466, 0.03951400190366, 0.0025705022650089};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {3.5307479877206e-05, 0.0021501657099222, 0.0021501657099222, 0.0022780322536957, 0.00012984962354834};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.2217109647723e-06, 0.00010701839286344, 0.00010701839286344, 0.00011942565209152, 4.7667975176802e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0050785818901811, 0.20199224921228, 0.20199224921228, 0.20668130223611, 0.019192405585065};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {9.3465519896578e-05, 0.0093445903313886, 0.0093445903313886, 0.01007048161087, 0.00033542327098299};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.2573886186184e-06, 0.00027258744549454, 0.00027258744549454, 0.00028674482283868, 8.78796720679e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.8986423865677e-08, 6.7933183967624e-06, 6.7933183967624e-06, 7.5959933079527e-06, 1.5551660664803e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0013031276550355, 0.10171900884347, 0.10171900884347, 0.10394180865623, 0.0057668918168994};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {1.1772981989102e-05, 0.0023526804373947, 0.0023526804373947, 0.0025398197945058, 4.2366403425413e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.4244717072173e-07, 3.4262168619715e-05, 3.4262168619715e-05, 3.5920437711806e-05, 5.6119361514784e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.2291920183738e-09, 4.2739401238296e-07, 4.2739401238296e-07, 4.7839232566858e-07, 4.8968061550081e-09};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.083977109377329, 0.9901881021392, 0.9901881021392, 1.1238980942085, 0.21109791813736};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.051228486408744, 0.56731684973577, 0.56731684973577, 0.81096077261257, 0.30280508622388};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0044323319040733, 0.052423026330272, 0.052423026330272, 0.065776706905395, 0.012566760873833};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.0019382614847568, 0.03351491083466, 0.03351491083466, 0.04514695497563, 0.030912126877778};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.035118000000724, 0.67888212771801, 0.67888212771801, 0.82938681088859, 0.095046668417497};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.011449738498237, 0.12322379174934, 0.12322379174934, 0.18067000714501, 0.030367409717608};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00087809667412559, 0.01519397687967, 0.01519397687967, 0.020201843494511, 0.0038768221919007};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.00011004736107595, 0.0016682725247329, 0.0016682725247329, 0.002297652413429, 0.0007480005262684};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0089315308289908, 0.34944171546552, 0.34944171546552, 0.41981260347941, 0.029745644682801};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0028709290232949, 0.031747220779055, 0.031747220779055, 0.047223450773604, 0.0063556998004923};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {5.7787447212004e-05, 0.0019414881844673, 0.0019414881844673, 0.0024935011415385, 0.00034300762987949};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {7.0145390033469e-06, 0.0001073022274249, 0.0001073022274249, 0.00014701266863273, 2.7512256368767e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0022740017766698, 0.17612146675339, 0.17612146675339, 0.20897132480524, 0.0095035438302642};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00071213974110974, 0.0080427539578396, 0.0080427539578396, 0.012088018943794, 0.0014781577409959};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {3.7189896562948e-06, 0.00024314194635372, 0.00024314194635372, 0.00030605242291051, 2.3396763002899e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {4.3902901571534e-07, 6.804616331063e-06, 6.804616331063e-06, 9.30595306694e-06, 1.1643127334482e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00057495080051447, 0.088310501720478, 0.088310501720478, 0.10413280017845, 0.0027824738835769};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00017748743719284, 0.0020214647179691, 0.0020214647179691, 0.0030529258678427, 0.0003500889198147};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.3551770270933e-07, 3.0352181774451e-05, 3.0352181774451e-05, 3.7858892086294e-05, 1.494927380759e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {2.7432376810659e-08, 4.2795113506559e-07, 4.2795113506559e-07, 5.8456806405933e-07, 6.3337413314102e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.1396336108712, 1.0761055615951, 1.0761055615951, 1.2416601280596, 0.34386264100162};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.069797115550398, 0.76674862941549, 0.76674862941549, 0.82718564483121, 0.4039426950479};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0065480685027187, 0.061259559360473, 0.061259559360473, 0.063406062691509, 0.015739940868661};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.0017626198126518, 0.037079573097843, 0.037079573097843, 0.038075455553361, 0.019977139624739};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.075742690321906, 0.78521215555536, 0.78521215555536, 0.92835452004339, 0.29053708486747};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.018021244195293, 0.17705806535985, 0.17705806535985, 0.19541573137294, 0.054312421816263};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0010674459687666, 0.018064732872127, 0.018064732872127, 0.018892449650416, 0.0033532462036121};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {8.3562998579764e-05, 0.0019646654504585, 0.0019646654504585, 0.002035888155947, 0.00053791019459314};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.01845824744019, 0.38790494396047, 0.38790494396047, 0.46866601191764, 0.083513318383076};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0044310452469873, 0.046119657039376, 0.046119657039376, 0.051770027304905, 0.010833427754286};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {7.1574806625237e-05, 0.0022283937552316, 0.0022283937552316, 0.0023538710061029, 0.00025950959558072};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {4.9392241613013e-06, 0.00012565640985951, 0.00012565640985951, 0.00013040954232974, 2.2315392046002e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0044802345783798, 0.19162569297644, 0.19162569297644, 0.23416687406776, 0.023309522782526};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0011155524844182, 0.011796117571444, 0.011796117571444, 0.013373811251923, 0.0024842537806508};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {4.7040571050592e-06, 0.00027373876483486, 0.00027373876483486, 0.00029123921981027, 1.7376069430197e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.0786906008602e-07, 7.8941230493156e-06, 7.8941230493156e-06, 8.1942482048381e-06, 9.3686060709253e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0011012478565411, 0.095216812028053, 0.095216812028053, 0.11696135651472, 0.0067407437042476};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00028076624185732, 0.002979532282937, 0.002979532282937, 0.0033945445507943, 0.00060221325195264};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {3.0207932194353e-07, 3.3873893536615e-05, 3.3873893536615e-05, 3.6193825859078e-05, 1.1534065291752e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.9362525319955e-08, 4.937907929525e-07, 4.937907929525e-07, 5.1261381991028e-07, 4.8041089883988e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case 6: // hump (moderate diffusion)
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.05598960539821, 0.502654311611746, 0.56651757959506, 0.515596526363158, 0.104871534533839};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00951855853016425, 0.14904962456944, 0.166851755545912, 0.160580177873316, 0.0215261024275195};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00100993636554085, 0.021552483664138, 0.0240960988347004, 0.0229327506068576, 0.00282384152105308};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {9.12729589711679e-05, 0.00332714650600316, 0.00371541312356028, 0.00413613315031279, 0.000207346905593986};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0145113217448486, 0.266477833019623, 0.275889956728542, 0.268629589695974, 0.0493074589108964};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00122951074891228, 0.0371627192581546, 0.0384082587154306, 0.0406494143601038, 0.00272802318048915};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {5.92110335051727e-05, 0.00236774288785203, 0.00244603228165149, 0.00243649472258649, 0.000230965994946708};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {6.10655665017985e-06, 0.000397000557886483, 0.000409973793822247, 0.000475730161539489, 1.38198496752451e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0036537202151042, 0.134888741971183, 0.136116135731423, 0.135244616433973, 0.0134377965332126};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000150361067062002, 0.00909596356681683, 0.00917459050904204, 0.0100044874746773, 0.000347691225183827};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {4.2863338573665e-06, 0.000341416900562107, 0.00034432891592876, 0.000354465193107174, 1.49376704756676e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.70551776962111e-07, 2.35809802598802e-05, 2.37796039927551e-05, 2.89709272335634e-05, 5.09533547202179e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000914810386715719, 0.0676056387887004, 0.067760644624282, 0.0676683805906934, 0.00337861156901409};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {1.87215323171807e-05, 0.00226665695954036, 0.00227159432766991, 0.00249673131787328, 4.20137712242308e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.78432726542267e-07, 4.34990765059434e-05, 4.35926452821543e-05, 4.45696198831107e-05, 9.47776421400537e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.1503328562008e-09, 1.45351877280207e-06, 1.45660384167193e-06, 1.79894845349868e-06, 1.49274581862091e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000228856074042296, 0.0338219850580569, 0.0338414113347459, 0.0338341920416948, 0.000844158294103003};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {2.33750441419278e-06, 0.000566256525027858, 0.00056656550004935, 0.000624036286765981, 5.13056996548356e-06};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.76886882009713e-08, 5.45210200475992e-06, 5.45504141629537e-06, 5.52560385067297e-06, 5.93817713889777e-08};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.59533368221956e-10, 9.06379327015383e-08, 9.06861235254559e-08, 1.1243434355227e-07, 4.81735332180655e-10};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0446304621138935, 0.500491679747327, 0.562208345783243, 0.531237052551703, 0.128925053802036};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0136050913548196, 0.146724683530567, 0.164896749698358, 0.186545581966799, 0.0329740987356719};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00172074901439671, 0.0215243764494506, 0.0241443926573147, 0.0256607386571741, 0.00371888102103318};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.000413686374099547, 0.00310121774395449, 0.00351003922778561, 0.00484956394895954, 0.000657941569290887};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0117687505866684, 0.265981797559421, 0.275126377700174, 0.2723360831902, 0.0401385334045724};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00310497580092992, 0.0370561216565311, 0.038513586905901, 0.0503452436184457, 0.00656181694847324};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {8.85937045760261e-05, 0.00236614060827754, 0.00244619075802975, 0.00258644923289987, 0.000380585106788511};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.86379176543096e-05, 0.000377849827402775, 0.000393949463584997, 0.000568427954019852, 7.04914393844236e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00295668575781197, 0.1348220410157, 0.136015445081859, 0.136053278011607, 0.0111697894537452};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00075724412790703, 0.00911950178217348, 0.00925833538157029, 0.0127074139668863, 0.0015908487058518};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {7.65887738984335e-06, 0.000340190198875278, 0.000343210173291652, 0.000384675487393518, 2.45692880988524e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {2.26400164973436e-06, 2.23107102799249e-05, 2.27249157020525e-05, 3.5803252593455e-05, 3.86947122321679e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000738154214927607, 0.0675968596777393, 0.0677475520384835, 0.0678401559318476, 0.00283970703918501};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000189786597258574, 0.00227624442422684, 0.00229680733178503, 0.00319793396790213, 0.000391671451033582};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {5.01780379836779e-07, 4.33415515085383e-05, 4.34388111923828e-05, 4.6924978925116e-05, 1.53721987966593e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.39913036727089e-07, 1.37528500077781e-06, 1.39233150639724e-06, 2.25920176350467e-06, 2.3590136127849e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000184291194533661, 0.0338207924820582, 0.033839674534022, 0.0338722616583047, 0.000710470213490577};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {4.76680995500489e-05, 0.00056892909646227, 0.00057320925989911, 0.000801842028662549, 9.71678599269721e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {3.17630092761798e-08, 5.43965826055893e-06, 5.44271914405467e-06, 5.68372392405349e-06, 9.85082776253088e-08};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {8.75649392992632e-09, 8.57722028221502e-08, 8.6706631530593e-08, 1.42059233922897e-07, 1.46446706428272e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0505015028661856, 0.516051939281752, 0.580827682443755, 0.572541348149995, 0.280569629229764};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0427321618177727, 0.182558281296701, 0.212788539778318, 0.210772725629267, 0.0758175469215297};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00311280739526581, 0.0243478790822661, 0.0275317245320847, 0.0253699401183943, 0.0063293627296736};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.000516535414320916, 0.00401418343867881, 0.00453846714509459, 0.00417498515892399, 0.00084026329789777};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0140537660377143, 0.271153490267926, 0.28069226002266, 0.282216335679423, 0.0863426336486305};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0152337248402966, 0.0523181411740378, 0.0582157836655675, 0.0652645597429588, 0.030834518219419};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.000154547766676123, 0.00254741731609847, 0.00263908089008523, 0.0025979002785563, 0.000553418262481509};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.52210998113096e-05, 0.000495479156159049, 0.000517508956469608, 0.000513356507622652, 0.000109218474110662};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00415127935925971, 0.136216382555481, 0.137483503817822, 0.139186903145695, 0.0204961837271932};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00448421876787009, 0.0136422452253054, 0.0151528386019601, 0.0176584027725763, 0.00916189167431675};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.51332525890028e-05, 0.000384519778072069, 0.000388324354928048, 0.000392657362655645, 3.49828049479167e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.40800713153052e-06, 3.13058173593628e-05, 3.19343652616928e-05, 3.24762118695535e-05, 6.00927627603109e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00111261490070375, 0.067948508809397, 0.0681100992089138, 0.0686497769098231, 0.00523385114689101};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00120542461090572, 0.00350677074282729, 0.00390615671418439, 0.00457032514593376, 0.00244900620760108};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {8.70989199997387e-07, 4.67757603727069e-05, 4.68906037726932e-05, 4.72689221341667e-05, 2.81218750387141e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {2.14335112611679e-07, 1.98721675917377e-06, 2.01435856625047e-06, 2.06080172801977e-06, 3.63274693992999e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000285618890857759, 0.033906031139915, 0.0339263598478095, 0.0340732081602604, 0.00132548427666779};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000311098128012209, 0.000890552064038958, 0.000993743885725773, 0.00115951438990603, 0.000629750983094479};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {5.00785952412292e-08, 5.66872423624932e-06, 5.67216263835661e-06, 5.69869322705275e-06, 1.89731529605386e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.34689747983545e-08, 1.25313852875165e-07, 1.26818951945955e-07, 1.29850170707445e-07, 2.25617116655918e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case 7: // boundary_layer_known (small diffusion)
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00399252407684254, 0.0385029856027784, 0.0439562455831579, 0.0749712521606381, 0.011107358536679};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000818224169422652, 0.0145235228578446, 0.0163204706983112, 0.0371939271393764, 0.00631316724086535};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.000292487752291059, 0.00744291438479758, 0.00843312737324695, 0.0107215404584998, 0.00111983114427596};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.20829535916964e-05, 0.0013262068045327, 0.00149455480834785, 0.00190703256433307, 0.000130060819179071};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00124756871259017, 0.0212371721851839, 0.0221934830009818, 0.0506483252123234, 0.00696435177251216};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000179902353809638, 0.00622123190802229, 0.00641884433456066, 0.010843459944037, 0.00138592522466196};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.34729074873717e-05, 0.00113354312562239, 0.00117690953416791, 0.00150736337431041, 0.000122392898199083};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.38731371941147e-06, 0.000108827858230469, 0.000113260396075851, 0.000150485673541057, 7.84625586328427e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000354677448396554, 0.0123176991325445, 0.0124596355212581, 0.0218221440405182, 0.00249596936304048};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {2.54624739369482e-05, 0.00178049249853369, 0.00179546921199812, 0.00259282432979103, 0.000231926019756935};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.60406769065816e-06, 0.000150368825107543, 0.00015188294603669, 0.00018698335487125, 1.04604220355613e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {4.9424322991487e-08, 7.7834522908449e-06, 7.86980571851273e-06, 1.05027976925295e-05, 3.76865518087605e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {9.48886484268027e-05, 0.00653217774296812, 0.00655086935308085, 0.00917120707211657, 0.000770490668408614};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {3.25640578197726e-06, 0.000464280637412514, 0.000465307357009532, 0.000607306121636699, 3.3776656675168e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.03417226526271e-07, 1.91065213706288e-05, 1.91556519239776e-05, 2.27702069736909e-05, 8.47896273156821e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.62351219585697e-09, 5.25596713578357e-07, 5.27087459495275e-07, 7.03331623179527e-07, 1.61593218181433e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {2.42616024846997e-05, 0.00332894782157319, 0.00333130813586195, 0.00406314717615907, 0.000229103009229912};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {4.09158103056109e-07, 0.000117794522375749, 0.000117861649761691, 0.000145082869982963, 4.56771113175716e-06};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {6.53797607250395e-09, 2.39849984467158e-06, 2.40005618355125e-06, 2.79275249298205e-06, 6.14718363484059e-08};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.17808232132323e-11, 3.42340407823315e-08, 3.42584482785896e-08, 4.57108541603342e-08, 5.97245719893835e-10};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00275778813576876, 0.0362634653578502, 0.0412600123000492, 0.0522550758686127, 0.0109792057565159};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000789068030501833, 0.0138005387021564, 0.0156186396845492, 0.0229300744429551, 0.00727873265996717};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.000309068218048751, 0.00698817187646164, 0.0080684814316604, 0.0103056399356863, 0.00267946397909814};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {4.80325133422138e-05, 0.00117634287301511, 0.00137465222419525, 0.00175987461053237, 0.000370164727739481};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00094505811266541, 0.0200785968685527, 0.0209595543620884, 0.031520717706692, 0.00736693125618708};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000187522586865911, 0.00564491093819837, 0.00585148199335319, 0.00834053198563398, 0.00177415323293208};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.90339128989017e-05, 0.00106548179189472, 0.00111411000030284, 0.00150829810669704, 0.000219345757845427};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {2.18767126443211e-06, 8.89101506329419e-05, 9.33572825654199e-05, 0.000129797469005172, 1.80585421993684e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000254417912197727, 0.0110775553894933, 0.0112001509516025, 0.0154370221392356, 0.00282228081466142};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {3.15255367146925e-05, 0.00158352459234299, 0.00159930209612142, 0.00225535539783242, 0.000309162837904983};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.97533388425727e-06, 0.000140401895199302, 0.000142011910993935, 0.000193538690922338, 2.47548287770055e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {9.85434579123517e-08, 5.88431744164692e-06, 5.96233413390905e-06, 8.48205141186388e-06, 9.87685440545071e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {6.49802507016075e-05, 0.00576781109122451, 0.0057838264103826, 0.00736380981842644, 0.000869589912338885};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {5.95708224341287e-06, 0.000407085233557674, 0.000408171567924917, 0.000572673569785418, 5.26095292871709e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.26995259780254e-07, 1.7733492994949e-05, 1.7782668347844e-05, 2.40343187298718e-05, 2.10594130574783e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.23433649973618e-09, 3.74208617503407e-07, 3.75493488153598e-07, 5.35970272703173e-07, 4.26473673906073e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {1.64463392447631e-05, 0.00292990309543439, 0.0029319389984662, 0.0035613335167395, 0.000241031703720642};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {1.33207667971846e-06, 0.000102456270160142, 0.000102532287396982, 0.000143564146258699, 8.21638531848254e-06};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {8.03556842543166e-09, 2.21778842102583e-06, 2.21928298379743e-06, 2.97709429457268e-06, 1.53643754397578e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.10278171542812e-10, 2.35262858582361e-08, 2.35481083509883e-08, 3.36056631525443e-08, 1.56414281353434e-09};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0047004958465528, 0.037819307618775, 0.0439065898983516, 0.0529549684105854, 0.0217958899954419};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00099652358895443, 0.0150065485138558, 0.0167927548049406, 0.01823613725217, 0.00915731608004333};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.000552441400906674, 0.00937419326424781, 0.0108918160844542, 0.0110106787665676, 0.00414014121335049};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {9.44077726001272e-05, 0.00159959477212554, 0.0018617753678949, 0.0018103664666223, 0.000697716050882248};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00136058050713842, 0.0204625471542998, 0.0213737297839001, 0.0272038311269896, 0.0101945784023904};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.000299482435408251, 0.00666162227856117, 0.00694051316259804, 0.00838133612902073, 0.00329862323848113};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {4.68988375713146e-05, 0.00135720924729253, 0.00142445070361666, 0.00158911580528107, 0.000553154841305503};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {3.39157504039242e-06, 0.000118126491042891, 0.000124471972655852, 0.000131782773032741, 4.33902088494336e-05};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.000418644297613721, 0.0109363549099457, 0.0110539723121714, 0.0143698443703247, 0.00401472963976737};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {5.84981913357418e-05, 0.00192063706928566, 0.00194424188287016, 0.00240708599434316, 0.000746550465039501};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {3.1858168718074e-06, 0.000174077027305137, 0.000176418053249593, 0.000203368614830373, 5.18290629162437e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.11366710061207e-07, 7.71981419294849e-06, 7.83315753991785e-06, 8.52423867114636e-06, 1.94509771464258e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00010848687194095, 0.00549977297205817, 0.00551450874133599, 0.00718668152220264, 0.00131809103908296};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {1.21928682638184e-05, 0.000503041661332088, 0.000504784780595233, 0.000627627246667824, 0.000135089262233343};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.02397682942518e-07, 2.16470057688955e-05, 2.1722257689949e-05, 2.52603263179395e-05, 3.97164057699425e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {4.08097883799623e-09, 4.87573089755155e-07, 4.89433754485615e-07, 5.34654743149213e-07, 7.28014478541027e-08};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {2.72882959661677e-05, 0.0027390722902625, 0.00274091819367906, 0.00356618056613611, 0.000376189704141007};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {2.83312188517108e-06, 0.000128233345829166, 0.000128370084370022, 0.000159576616148521, 2.30259216634826e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.26822700643821e-08, 2.68411157783231e-06, 2.68647546633154e-06, 3.13249879552234e-06, 2.74987834938084e-07};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.87629987824608e-10, 3.05506804431144e-08, 3.05806095330498e-08, 3.33674879057267e-08, 2.74284649593587e-09};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case 11: // convection_reaction (pure convection-reaction)
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.087423482379532, 1.0979652764761, 0.58896892318376, 0.34016084681932, 0.51640100183487};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.01332225334409, 0.2904429871326, 0.15793912702429, 0.061572524540975, 0.083897436880195};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0020562129904171, 0.059673539813059, 0.029168779262219, 0.010103310635975, 0.01973002206286};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {0.00061012776597731, 0.025562829083026, 0.011787324868779, 0.0032923824505361, 0.0058625734194051};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.024499405895532, 0.5811418983634, 0.22755186317643, 0.1223954241572, 0.20289810781337};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0017614975947925, 0.073014750155763, 0.027954202066227, 0.010680829596382, 0.016409507738718};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00024592929968886, 0.015078411729212, 0.0051275836924041, 0.0016007990710562, 0.0038388120618769};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.8542229779542e-05, 0.001638049478061, 0.00058846823592249, 0.0001378153382152, 0.00020032651286606};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0064906391870304, 0.29637883898507, 0.085227710255762, 0.043559505724286, 0.061837747588871};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00025511813738949, 0.020909274551153, 0.0055824234009777, 0.0020819278582188, 0.0033871284375167};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.6505921791508e-05, 0.0020653261393961, 0.00052807465455455, 0.00014780227782602, 0.00031896852934254};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.7520143666239e-07, 0.0001017598713327, 2.6436894723621e-05, 5.8972209422985e-06, 8.3617009851711e-06};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0016778865575375, 0.15018427975698, 0.031232519130537, 0.015536133067538, 0.016540411242226};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {3.3418429854439e-05, 0.0054756835566635, 0.0010475191276892, 0.00037768624517814, 0.00063457803633549};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {1.0559287940318e-06, 0.000266253452525, 4.937683771072e-05, 1.317836214738e-05, 2.1072716364312e-05};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {1.8364322812337e-08, 6.5124579038875e-06, 1.2089857730432e-06, 2.6233886363512e-07, 2.6907574235402e-07};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00042707936362827, 0.075664097481038, 0.011267335631879, 0.0055231776218506, 0.0043671766666171};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {4.2501203269782e-06, 0.0013931038162513, 0.00019015861765895, 6.7304175558393e-05, 9.4815654063218e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {6.6699950657639e-08, 3.3759749666875e-05, 4.4795236539695e-06, 1.1687178958549e-06, 1.364858176256e-06};
                    break;
                  case 4: // P/Q4 polynomials
                    errors = {5.7970017572841e-10, 4.1189080367232e-07, 5.4438144279915e-08, 1.1629638804544e-08, 8.7482652233462e-09};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      default:
        ErrThrow("cd_dg_test.C : get_reference_values() : problem_number ",
            problem_number, " not known.");
    } //endswitch problem_number
  }
  else if (d == 3)
  {
    switch (problem_number)
    {
      case 0: // SineLaplace
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.35206175004723, 1.856589649811, 0.0001856589649811, std::nan(""), 0.86325509718742};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.14073395371412, 1.0478532758664, 0.00010478532758664, std::nan(""), 0.24151088147567};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.084837676548199, 0.84473294829606, 8.4473294829606e-05, std::nan(""), 0.30768579249563};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.21215208434548, 1.4432911172939, 0.00014432911172939, std::nan(""), 0.4636981343076};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.046506179851747, 0.58307203050061, 5.8307203050061e-05, std::nan(""), 0.15814669020822};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0172761613414, 0.2832537619133, 2.832537619133e-05, std::nan(""), 0.077792167454603};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.088232602102826, 0.89499050026551, 8.9499050026551e-05, std::nan(""), 0.30244409161787};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0085147699020816, 0.22640400940985, 2.2640400940985e-05, std::nan(""), 0.02911859773228};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0010230582084041, 0.037823790652798, 3.7823790652798e-06, std::nan(""), 0.0074665387661494};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.030904096143647, 0.51129419102147, 5.1129419102148e-05, std::nan(""), 0.11770126839608};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00093068878752139, 0.061090154309443, 6.1090154309443e-06, std::nan(""), 0.0038524760326016};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {6.329249762322e-05, 0.0050674164979956, 5.0674164979956e-07, std::nan(""), 0.00056702409329457};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.22583812525174, 1.4798117750743, 0.00014798117750743, std::nan(""), 0.54692008827167};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.083676349380838, 0.94525455599173, 9.4525455599173e-05, std::nan(""), 0.20227541070425};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.081009124441816, 0.74362212959555, 7.4362212959555e-05, std::nan(""), 0.55843074507599};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.08806126862545, 1.0635960482825, 0.00010635960482825, std::nan(""), 0.2578087308804};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.030825487244349, 0.42247077894972, 4.2247077894972e-05, std::nan(""), 0.075775626260727};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.014616714651172, 0.1998015174011, 1.998015174011e-05, std::nan(""), 0.30031185440172};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.025792541639708, 0.5898272474631, 5.898272474631e-05, std::nan(""), 0.10109302992495};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0092627458663032, 0.1385118783565, 1.385118783565e-05, std::nan(""), 0.024715670189033};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0011013968162084, 0.026609901824959, 2.6609901824959e-06, std::nan(""), 0.02942474991151};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0074316919469813, 0.30674679846986, 3.0674679846986e-05, std::nan(""), 0.033989179905379};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0022819046532897, 0.035875423235305, 3.5875423235305e-06, std::nan(""), 0.006392834684992};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {7.489112239036e-05, 0.0035217148794558, 3.5217148794558e-07, std::nan(""), 0.0020810841183464};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.23794188695837, 1.4165808651474, 0.00014165808651474, std::nan(""), 0.64477837980164};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.15520945057699, 1.1993287772449, 0.00011993287772449, std::nan(""), 0.3666525901821};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.092155920694458, 0.88177943797002, 8.8177943797002e-05, std::nan(""), 0.95774530222004};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.10159435780273, 1.088014306024, 0.0001088014306024, std::nan(""), 0.45410163272544};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.033551938231291, 0.49615661818172, 4.9615661818172e-05, std::nan(""), 0.11317054999584};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.015042076868901, 0.23714543059344, 2.3714543059344e-05, std::nan(""), 0.41299164069895};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.030978851184942, 0.60511044192728, 6.0511044192728e-05, std::nan(""), 0.20655602847024};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0067925597131422, 0.15330242007534, 1.5330242007534e-05, std::nan(""), 0.022936552131587};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0011466333567921, 0.030994614423613, 3.0994614423613e-06, std::nan(""), 0.038136189584328};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0079843692424019, 0.30878661446384, 3.0878661446384e-05, std::nan(""), 0.07120029653151};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0014292484947638, 0.039251660168181, 3.9251660168181e-06, std::nan(""), 0.0048207171507069};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {7.7963544586857e-05, 0.0040589429302455, 4.0589429302455e-07, std::nan(""), 0.0027554425821701};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case -4: // Polynomial of degree 6
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {1.3011691743936, 10.886917505714, 10.357323550552, std::nan(""), 4.3889855692852};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.59240808643566, 5.1718810768483, 4.3978537395822, std::nan(""), 1.3774700706174};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.64649393925691, 7.2546114160447, 4.9399722791038, std::nan(""), 1.8466636612532};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.10377439208612, 1.7995326125343, 1.3576812305608, std::nan(""), 0.25795906101779};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.17408183758355, 3.8859432577068, 2.0613047836222, std::nan(""), 0.58974849245104};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.01534946529265, 0.48990328486649, 0.29297729163909, std::nan(""), 0.04653494742622};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.040765609817759, 1.9104969300077, 0.80209901114738, std::nan(""), 0.18604401204674};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0020342493509971, 0.11966824574972, 0.054841575820163, std::nan(""), 0.0074766506257984};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {1.018641497064, 10.445326008318, 9.9335491897595, std::nan(""), 3.4195908315533};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.40398530023032, 5.0546999425252, 4.0453664614788, std::nan(""), 1.0858117667663};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.43493403833075, 6.7403254529374, 4.4823192393927, std::nan(""), 1.4952396387053};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.07308046702667, 1.7096818103635, 1.2627939728814, std::nan(""), 0.18801308690813};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.12242583600358, 3.6628224890545, 1.9021090644146, std::nan(""), 0.48451643708648};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.010491336447183, 0.47710777498581, 0.28031651818226, std::nan(""), 0.036963329168083};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.030485253393299, 1.8637063076442, 0.77613837596702, std::nan(""), 0.14939348192348};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0016288658827753, 0.12367989807642, 0.05637821782058, std::nan(""), 0.0072652248019858};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.98838503598489, 10.418794118575, 9.9518503309285, std::nan(""), 3.2377488971465};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.40358781395673, 5.156360770924, 4.0629581370236, std::nan(""), 1.2877610242125};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.4189615722019, 6.695494527964, 4.477392959172, std::nan(""), 1.7453047017679};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.079657004931151, 1.7658116825703, 1.2893327984127, std::nan(""), 0.24251584962315};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.1529677917408, 3.6382388318724, 1.8992114840245, std::nan(""), 0.72950543835017};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.017720314029408, 0.52414306915693, 0.30327822243358, std::nan(""), 0.071590484700422};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.063792136809077, 1.8764662977353, 0.78933432464204, std::nan(""), 0.3609352256056};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0049396206432514, 0.14940971513666, 0.067684554509783, std::nan(""), 0.021630345577049};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case 4: // Boundary Layer known
        switch (symmetry)
        {
          case 1: // SIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0028624386375823, 0.016252984398517, 0.01680595816383, std::nan(""), 0.0073821916476072};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00026127222355276, 0.0024946935435422, 0.0026279923065599, std::nan(""), 0.00090174891487296};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00082423407183164, 0.0080062758539485, 0.0081021033287357, std::nan(""), 0.0022540734051107};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {3.8968102393672e-05, 0.00061376622324048, 0.00062699015945462, std::nan(""), 0.0001441890970653};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00020103152402943, 0.0038137580654387, 0.0038266167065231, std::nan(""), 0.00068499158334266};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {5.3328150551052e-06, 0.00015137987658823, 0.00015234538712972, std::nan(""), 2.3728831342914e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {4.9836492656735e-05, 0.0018808868702312, 0.0018825201472975, std::nan(""), 0.00018726353543166};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {6.8211075832759e-07, 3.7673897635021e-05, 3.7737225559961e-05, std::nan(""), 3.7636750872215e-06};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case 0: // IIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.001780160742726, 0.01519452361984, 0.015782819680921, std::nan(""), 0.0049261786284417};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00015762124803009, 0.0022540040335378, 0.0024265795289309, std::nan(""), 0.00057401587415014};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00047982837224909, 0.0076694320302836, 0.0077661082954207, std::nan(""), 0.0025402872469123};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {2.5299889699908e-05, 0.00063497142796854, 0.00065042439739935, std::nan(""), 9.6485805857296e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00012139322850528, 0.0037722884324888, 0.0037851322948375, std::nan(""), 0.00058365132435916};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {3.9955444856857e-06, 0.00016471631285793, 0.00016580397991368, std::nan(""), 1.7212598351787e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {3.0403596570815e-05, 0.0018757806163871, 0.0018774129033979, std::nan(""), 0.00014088208230054};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {7.6814433120387e-07, 4.1665829342724e-05, 4.1736179666161e-05, std::nan(""), 3.5031432672294e-06};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0065838621245357, 0.019578651034851, 0.020124582941165, std::nan(""), 0.014263275406382};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00021404274286434, 0.0025750384456379, 0.002748776483968, std::nan(""), 0.0010769386206616};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0018097355169578, 0.0095348551453583, 0.0096268476788042, std::nan(""), 0.0078095801779484};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {5.3126341443012e-05, 0.00081364267282381, 0.00083060941304846, std::nan(""), 0.00034376711364279};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00053842683795964, 0.0043420602463418, 0.0043547709785429, std::nan(""), 0.0021345795331762};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {1.5326039284534e-05, 0.0002367871602668, 0.00023810903306703, std::nan(""), 6.957044863419e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00014700776149376, 0.0020292338953325, 0.0020308887249832, std::nan(""), 0.00070368998859656};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {4.1661468107237e-06, 6.3023369119504e-05, 6.3113052983547e-05, std::nan(""), 1.8228997480024e-05};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      case 6: // Pure Convection
        switch (symmetry)
        {
          case 1: // SIPG
          case 0: // IIPG
          case -1: // NIPG
            switch (refinement_level)
            {
              case 0: // 0th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.090667369192631, 0.86511418951215, 1.1688130649559, std::nan(""), 0.21237662484046};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.084330021329121, 1.0088848597704, 1.1411363630383, std::nan(""), 0.21010664206363};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.039275544691743, 0.64060116485732, 0.71106675510599, std::nan(""), 0.33335574518472};
                    break;
                  default:
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 1: // 1st refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.0600613962730523,0.789554620169025,0.738299925876695,std::nan(""),0.276727825316605};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.012924358681841, 0.35666903251546, 0.29879067947767, std::nan(""), 0.066140689908981};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.0040183927716804, 0.12673810713789, 0.09484716172786, std::nan(""), 0.11892026182763};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 2: // 2nd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.018941688919937, 0.48761052029666, 0.30350895035078, std::nan(""), 0.20424639840832};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.0022728739838021, 0.12849774974202, 0.07511773020682, std::nan(""), 0.022061557026999};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {0.00032284225743228, 0.020568391745844, 0.011602502470302, std::nan(""), 0.019234709641132};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 3: // 3rd refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {0.00519706435817057,0.262434764076368,0.11474822441519,std::nan(""),0.112843313062794};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {0.00029208448845407, 0.033241602885352, 0.013738183587122, std::nan(""), 0.0033757471295268};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {2.5676869031305e-05, 0.0032817652131559, 0.0012243985935782, std::nan(""), 0.0036669473311217};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              case 4: // 4th refinement
                switch (degree)
                {
                  case 1: // P/Q1 polynomials
                    errors = {};
                    break;
                  case 2: // P/Q2 polynomials
                    errors = {};
                    break;
                  case 3: // P/Q3 polynomials
                    errors = {};
                    break;
                    ErrThrow("cd_dg_test.C : get_reference_values() : ",
                        "Degree not implemented yet.");
                } // endswitch degree
                break;
              default:
                ErrThrow("cd_dg_test.C : get_reference_values() : Reference ",
                    "values are only known up to refinement level 5. You chose ",
                    refinement_level);
            } // endswitch refinement_level
            break;
          default:
            ErrThrow("cd_dg_test.C : get_reference_values() : DG_SYMMETRY ",
                symmetry, " is not in {1, 0, -1}.");
        }
        break;
      default:
        ErrThrow("cd_dg_test.C : get_reference_values() : problem_number ",
            problem_number, " not known.");
    } //endswitch problem_number
  }
  else
  {
    ErrThrow("Unknown dimension for reference values");
  }
} // end get_reference_values
