/**
 * @brief A test program for the testing the slope limiter class.
 *
 * This serves as a test for the slope limiter class. It is intended to test all
 * features of the class.
 *
 * All the basic features of the class as constructors, getters and so on are
 * tested. In the end the different limiters are tested to behave correctly.
 *
 * Should this test fail, there are two possibilities: either you made a mistake
 * which broke the programs functionality. Then you must find the mistake.  Or
 * you changed some program setup (e.g. changed the default solver). Then this
 * tests  shows you how many other program parts are affected by your changes.
 * If you are not perfectly sure how to repair this, it is a good idea to
 * describe your changes in the forum and request support.
 *
 *
 * @date 2021/08
 * @author Derk Frerichs
 */

#include <cmath>
#include <memory>

#include "MainUtilities.h"
#include "MooNMD_Io.h"
#include "Domain.h"
#include "Database.h"
#include "SlopeLimiter.h"
#include "ParMooN.h"
#include "ParMooN_repository_info.h"

#ifdef __3D__
constexpr int d = 3;
#else
constexpr int d = 2;
#endif

void test_default_SL_database();
void test_constructor();
void test_info();
void test_getters();
void test_different_limiters(const bool& testall);

void check_pw_constant_function(ParameterDatabase parmoon_db);
void check_linear_function(ParameterDatabase parmoon_db);
void linear_function(double x, double y, double* values);
void check_quadratic_function(ParameterDatabase parmoon_db);
void quadratic_function(double x, double y, double* values);
void check_2x_limiter(ParameterDatabase parmoon_db);
void modulus_function(double x, double y, double* values);

bool are_vectors_equal(std::vector<double> values_orig,
    std::vector<double> values_new);

const std::string path_to_repo = parmoon::source_directory;
const std::string path_to_meshes = path_to_repo + "/data/mesh/";
const int max_n_refinements = 4;


int main(int, char* argv[])
{
  parmoon::parmoon_initialize();
  Output::print("\n        START ",d,"D TEST");
  Output::print("------------------------------\n\n");

  bool testall = false;
  if (argv[1])
  {
    testall = (std::string(argv[1]).compare("testall") == 0);
  }

  // Increase verbosity to get more information about the tests
  Output::setVerbosity(1);

  TDatabase::ParamDB->USE_ISOPARAMETRIC = 0;

  // test different features of slope limiter class
  Output::print("Start checking default slope limiter database");
  test_default_SL_database();
  Output::print("-- Passed\n");

  Output::print("Start checking constructor");
  test_constructor();
  Output::print("-- Passed\n");

  Output::print("Start checking info");
  test_info();
  Output::print("-- Passed\n");

  Output::print("Start checking getters");
  test_getters();
  Output::print("-- Passed\n");

  Output::print("Start checking behaviour of limiters");
  test_different_limiters(testall);
  Output::print("-- Passed");

  Output::print("\n\n------------------------------");
  Output::print("        ",d,"D TEST PASSED\n");

  parmoon::parmoon_finalize();
}


// Check if all parameters needed for slope limiters are in the default slope
// limiter database
void test_default_SL_database()
{
  auto db = TSlopeLimiter::default_slope_limiter_database();
  std::array<std::string, 9> parameters = {"apply_limiter", "limiter_name",
    "M_lim", "gamma_limiter", "alpha_ref", "characteristic_length",
    "characteristic_solution_scale", "C0_CJM", "C1_CJN"};
  for (auto para : parameters)
  {
    if (!db.contains(para))
    {
      ErrThrow("\"", para, "\" does not exist in database!");
    }
  }
}

// Test if the constructor works properly
void test_constructor()
{

  auto db = TSlopeLimiter::default_slope_limiter_database();
  TSlopeLimiter limiter(db);

  if (!db["limiter_name"].is( limiter.get_limiter_name() ))
  {
    ErrThrow("Constructed object is limiter ",limiter.get_limiter_name(),
        " but should be ", db["limiter_name"], "." );
  }
  if (!db["M_lim"].is( limiter.get_m_lim() ))
  {
    ErrThrow("Constructed object has M_lim ", limiter.get_m_lim(),
        " but should be ", db["M_lim"], "." );
  }
  if (!db["gamma_limiter"].is( limiter.get_gamma() ))
  {
    ErrThrow("Constructed object has gamma_limiter ",limiter.get_gamma(),
        " but should be ", db["gamma_limiter"], "." );
  }
  if (!db["alpha_ref"].is( limiter.get_alpha_ref() ))
  {
    ErrThrow("Constructed object has alpha_ref ",limiter.get_alpha_ref(),
        " but should be ", db["alpha_ref"], "." );
  }
  if (!db["C0_CJM"].is( limiter.get_C0_CJM() ))
  {
    ErrThrow("Constructed object has C0_CJM ",limiter.get_C0_CJM(),
        " but should be ", db["C0_CJM"], "." );
  }
  if (!db["characteristic_length"].is( limiter.get_char_length() ))
  {
    ErrThrow("Constructed object has characteristic_length ",
        limiter.get_char_length(), " but should be ",
        db["characteristic_length"], "." );
  }
  if (!db["characteristic_solution_scale"].is( limiter.get_char_sol_scale() ))
  {
    ErrThrow("Constructed object alpha_ref ",limiter.get_char_sol_scale(),
        " but should be ", db["characteristic_solution_scale"], "." );
  }
  if (!db["C1_CJN"].is( limiter.get_C1_CJN() ))
  {
    ErrThrow("Constructed object has C1_CJN ",limiter.get_C1_CJN(),
        " but should be ", db["C1_CJN"], "." );
  }
}

void test_getters()
{
  auto db = TSlopeLimiter::default_slope_limiter_database();
  double m_lim = 42;
  db["M_lim"] = m_lim;
  auto gamma = 3.14159265358979323846264338327950288419716939937510582097494459;
  db["gamma_limiter"] = gamma;
  auto alpha_ref = 2.7182818284590452353602874713526624977572470936999595749669;
  db["alpha_ref"] = alpha_ref;
  auto C0_CJM = 1.6180339887498948482045868343656381177203091798057628621354486;
  db["C0_CJM"] = C0_CJM;
  auto characteristic_length = 43252003274489856000.;
  db["characteristic_length"] = characteristic_length;
  auto characteristic_solution_scale = 15;
  db["characteristic_solution_scale"] = characteristic_solution_scale;
  auto C1_CJN = 37;
  db["C1_CJN"] = C1_CJN;
  TSlopeLimiter limiter(db);

  if ( m_lim != limiter.get_m_lim() )
  {
    ErrThrow("Constructed object has M_lim ", limiter.get_m_lim(),
        " but should be ", m_lim, "." );
  }
  if (gamma != limiter.get_gamma() )
  {
    ErrThrow("Constructed object has gamma_limiter ", limiter.get_gamma(),
        " but should be ", gamma, "." );
  }
  if (alpha_ref != limiter.get_alpha_ref() )
  {
    ErrThrow("Constructed object has alpha_ref ", limiter.get_alpha_ref(),
        " but should be ", alpha_ref, "." );
  }
  if (C0_CJM != limiter.get_C0_CJM() )
  {
    ErrThrow("Constructed object has C0_CJM ", limiter.get_C0_CJM(),
        " but should be ", C0_CJM, "." );
  }
  if (characteristic_length != limiter.get_char_length() )
  {
    ErrThrow("Constructed object has characteristic_length ",
        limiter.get_char_length(), " but should be ", characteristic_length,
        ".");
  }
  if (characteristic_solution_scale != limiter.get_char_sol_scale() )
  {
    ErrThrow("Constructed object has characteristic_solution_scale ",
        limiter.get_char_sol_scale(), " but should be ",
        characteristic_solution_scale, "." );
  }
  if (C1_CJN != limiter.get_C1_CJN() )
  {
    ErrThrow("Constructed object has C1_CJN ",
        limiter.get_C1_CJN(), " but should be ", C1_CJN, "." );
  }

  TFEFunction2D fe_func;
  auto features = limiter.get_features(&fe_func);
  if (!features.empty())
  {
    ErrThrow("Feature vector should be empty but is not.");
  }
  auto is_cell_to_limit = limiter.get_is_cell_to_limit(&fe_func);
  if (!is_cell_to_limit.empty())
  {
    ErrThrow("is_cell_to_limit vector should be empty but is not.");
  }

}

void test_info()
{
  // The info method prints information about the limiter. Depending on if the
  // limiter was already used or not it prints different information. Therefore
  // we first prepare a FE function such on which the limiter can be applied.
  auto db = TSlopeLimiter::default_slope_limiter_database();
  db["apply_limiter"] = true;
  db["limiter_name"] = "LinTriaReco";

  db.merge(TDomain::default_domain_parameters());
  db["geo_file"] = "TwoTriangles";
  db["refinement_n_initial_steps"] = 2;
  Output::set_script_mode(true);
  TDomain domain(db);
  Output::set_script_mode(false);
  domain.refine_and_get_hierarchy_of_collections(db);
  std::shared_ptr<TFESpace2D> fespace (new
      TFESpace2D(domain.GetCollection(It_Finest, 0), "",
        BoundConditionNoBoundCondition, -11));
  std::vector<double> values(fespace->get_n_dof(), 0);
  TFEFunction2D fefunc(fespace, "", values.data());
  fefunc.Interpolate(quadratic_function);

  // Create limiter and see what info method does.
  TSlopeLimiter limiter(db);
  limiter.info();
  auto dummy = limiter.get_is_cell_to_limit(&fefunc);
  limiter.info();
}

void test_different_limiters(const bool& testall)
{
  auto db = TSlopeLimiter::default_slope_limiter_database();
  db.merge(TDomain::default_domain_parameters());
  db["apply_limiter"].set(true, true);
  db["M_lim"].set(1, true);
  db["gamma_limiter"].set(1, true);
  db["alpha_ref"].set(4, true);
  db["C0_CJM"].set(1, true);
  db["characteristic_length"].set(1, true);
  db["characteristic_solution_scale"].set(1, true);
  db["C1_CJN"].set(-1, true);

  // This function calls all the single tests.
  std::vector<std::string> geometries = {"UnitSquare_quads.mesh",
    "UnitSquare_tria.mesh", "UnitSquareIrregular.GEO", "UnitSquareWithSpherical"
      "InscribedRegion.mesh"};
  std::vector<std::string> boundaries(geometries.size(), "UnitSquare.PRM");
  if (!testall)
  {
    geometries.resize(2);
    boundaries.resize(geometries.size());
  }
  int max_approx_degree = (testall) ? 4 : 2;
  db.add("max_degree", max_approx_degree, " ");

  int start_geom_i = 0;
  db.add("geometry_nr", start_geom_i, " ", start_geom_i, (int)
      geometries.size());

  for (unsigned int geom_i = start_geom_i; geom_i < geometries.size(); ++geom_i)
  {
    db["geometry_nr"].set((int) geom_i, true);
    db["geo_file"].set(path_to_meshes+geometries[geom_i], false);
    db["boundary_file"].set(path_to_meshes+boundaries[geom_i], false);
    Output::print<1>("\nStart checking geometry ", geom_i, "\n");

    std::vector<std::string> limiters = {"Galerkin", "LinTriaReco",
      "ConstTriaReco", "LinQuadReco", "ConstQuadReco", "LinQuadDeriv",
      "ConstQuadDeriv", "ConstJump", "ConstJumpMod", "ConstJumpL1Norm",
      "ConstJumpL2Norm", "ConstJumpLinftyNorm"};

    for (auto method : limiters)
    {
      if ((method == "LinQuadDeriv" || method == "ConstQuadDeriv") &&
          geom_i != 0)
      {
        // These limiter do only make sense on quadliteral meshes
        continue;
      }
      Output::print<2>("Start checking limiter ", method);
      db["limiter_name"].set(method, true);

      if (method != "ConstJumpL1Norm" && method != "ConstJumpL2Norm" && method
          != "ConstJumpLinftyNorm")
      {
      Output::print<2>("\n", method, ": Start checking piecewise constant ",
          "functions");
      check_pw_constant_function(db);
      Output::print<2>("Checked a piecewise constant function");
      }

      if (! (geom_i > 1 && (method == "ConstTriaReco" ||
              method == "ConstQuadReco")))
      {
        Output::print<2>("\n", method, ": Start checking linear function");
        check_linear_function(db);
        Output::print<2>("Checked a linear function");
      }

      if (method == "ConstJump" || method ==  "ConstJumpMod" || method ==
          "ConstJumpL1Norm" || method == "ConstJumpL2Norm" || method ==
          "ConstJumpLinftyNorm")
      {
        Output::print<2>("\n", method, ": Start checking quadratic function");
        check_quadratic_function(db);
        Output::print<2>("Checked a quadratic function");
      }
      else
      {
        Output::print<2>("\n", method, ": Start checking 2x limiter");
        check_2x_limiter(db);
        Output::print<2>("Checked 2x limiter function");
      }
      Output::print<2>('\n');
    }
  }
}

void check_pw_constant_function(ParameterDatabase parmoon_db)
{
  // This function takes a piecewise constant function that should not be
  // limited and checks if it gets limited
  Output::set_script_mode(true);
  TDomain domain(parmoon_db);
  Output::set_script_mode(false);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  domain.RegRefineAll();
  for (int ref_i = 1; ref_i < max_n_refinements; ++ref_i)
  {
    for (int degree = 1; degree <=parmoon_db["max_degree"].get<int>(); ++degree)
    {
      std::shared_ptr<TFESpace2D> fespace (new
          TFESpace2D(domain.GetCollection(It_Finest, 0), "",
            BoundConditionNoBoundCondition, -10-degree));
      std::vector<double> values(fespace->get_n_dof(), 0);
      TFEFunction2D fefunc(fespace, "", values.data());
      auto n_loc_dof = fespace->get_n_local_dof(0);
      auto n_cells = domain.GetCollection(It_Finest, 0)->GetN_Cells();
      for (int cell_i = 0; cell_i < n_cells; ++cell_i)
      { // piecewise constant function
        double divisor = (n_cells > 1) ? n_cells - 1 : 1;
        values[n_loc_dof * cell_i] =
          std::pow(-1,cell_i) * (double) cell_i/divisor * (ref_i+1);
      }

      auto val_orig = values;
      TSlopeLimiter limiter(parmoon_db);
      limiter.limit_function(&fefunc);
      if (!are_vectors_equal(values, val_orig))
      {
        ErrThrow("The limiter has limited a piecewise constant function!");
      }
    }
    domain.RegRefineAll();
  }
}

void check_linear_function(ParameterDatabase parmoon_db)
{
  // This function takes a linear function and checks if it gets limited (should
  // not)
  Output::set_script_mode(true);
  TDomain domain(parmoon_db);
  Output::set_script_mode(false);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  domain.RegRefineAll();
  for (int ref_i = 1; ref_i < max_n_refinements; ++ref_i)
  {
    for (int degree = 1; degree <=parmoon_db["max_degree"].get<int>(); ++degree)
    {
      std::shared_ptr<TFESpace2D> fespace (new
          TFESpace2D(domain.GetCollection(It_Finest, 0), "",
            BoundConditionNoBoundCondition, -10-degree));
      std::vector<double> values(fespace->get_n_dof(), 0);
      TFEFunction2D fefunc(fespace, "", values.data());
      fefunc.Interpolate(linear_function);

      auto fe_function_val_orig = values;
      TSlopeLimiter limiter(parmoon_db);
      limiter.limit_function(&fefunc);
      if (!are_vectors_equal(fe_function_val_orig, values))
      {
        ErrThrow("The limiter has limited a linear function!");
      }
    }
    domain.RegRefineAll();
  }
}

void linear_function(double x, double y, double *values)
{
  values[0] = 4 * x + 2 * y - 3.14;
}

void check_quadratic_function(ParameterDatabase parmoon_db)
{
  // This function takes a quadratic function and checks if it gets limited
  // (should not)
  Output::set_script_mode(true);
  TDomain domain(parmoon_db);
  Output::set_script_mode(false);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  domain.RegRefineAll();
  for (int ref_i = 1; ref_i < max_n_refinements; ++ref_i)
  {
    for (int degree = 2; degree <=parmoon_db["max_degree"].get<int>(); ++degree)
    {
      std::shared_ptr<TFESpace2D> fespace (new
          TFESpace2D(domain.GetCollection(It_Finest, 0), "",
            BoundConditionNoBoundCondition, -10-degree));
      std::vector<double> values(fespace->get_n_dof(), 0);
      TFEFunction2D fefunc(fespace, "", values.data());
      fefunc.Interpolate(quadratic_function);

      auto fe_function_val_orig = values;
      TSlopeLimiter limiter(parmoon_db);
      limiter.limit_function(&fefunc);
      if (!are_vectors_equal(fe_function_val_orig, values))
      {
        ErrThrow("The limiter has limited a quadratic function!");
      }
    }
    domain.RegRefineAll();
  }
}

void quadratic_function(double x, double y, double *values)
{
  values[0] = 10 * (x * x + y * y + x + y) - 10;
}

void check_2x_limiter(ParameterDatabase parmoon_db)
{
  // This function checks what happens if the limiter is applied twice. The
  // second application should not change the solution anymore
  Output::set_script_mode(true);
  TDomain domain(parmoon_db);
  Output::set_script_mode(false);
  domain.refine_and_get_hierarchy_of_collections(parmoon_db);
  for (int ref_i = 0; ref_i < max_n_refinements; ++ref_i)
  {
    for (int degree = 1; degree <=parmoon_db["max_degree"].get<int>(); ++degree)
    {
      std::shared_ptr<TFESpace2D> fespace (new
          TFESpace2D(domain.GetCollection(It_Finest, 0), "",
            BoundConditionNoBoundCondition, -10-degree));
      std::vector<double> values(fespace->get_n_dof(), 0);
      TFEFunction2D fefunc(fespace, "", values.data());
      fefunc.Interpolate(modulus_function);

      TSlopeLimiter limiter(parmoon_db);
      limiter.limit_function(&fefunc);
      auto fe_function_val_orig = values;
      limiter.limit_function(&fefunc);
      if (!are_vectors_equal(fe_function_val_orig, values))
      {
        ErrThrow("The limiter has limited an already limited function!");
      }
    }
    domain.RegRefineAll();
  }
}

void modulus_function(double x, double y, double *values)
{
  values[0] = std::abs(x + y - 1);
  if (x > 0.5 && y > 0.5)
  {
    values[0] = - values[0];
  }
}

bool are_vectors_equal(std::vector<double> vec_1, std::vector<double> vec_2)
{
  // This function checks whether vectors are equal up to a tolerance of 1e-10
  bool are_equal = true;
  double tol = 1e-10;
  if (vec_1.size() != vec_2.size())
  {
    are_equal = false;
    return are_equal;
  }
  for (unsigned int dof_i = 0; dof_i < vec_1.size(); ++dof_i)
  {
    if (vec_2[dof_i] < vec_1[dof_i] - tol ||
        vec_2[dof_i] > vec_1[dof_i] + tol)
    {
      are_equal = false;
      break;
    }
  }
  return are_equal;
}
