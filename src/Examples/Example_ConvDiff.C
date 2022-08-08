#include <cmath> // e.g. sin
#include <Example.h>
#include <MooNMD_Io.h>

using namespace parmoon;

ParameterDatabase get_default_ConvDiff_Example_parameters(int example)
{
  ParameterDatabase db("default example database");
  db.add("example", example,
         "Choose which example to run. Depending on the example other values "
         "parameters are also taken into acount", -5, 200);
  db.add("diffusion_coefficient", 1.0,
         "The (scalar, constant) coefficient in front of the Laplace term in a "
         "convection-diffusion-reaction problem.", 0., 100.);
  // other parameters may be added depending on the chosen `example`.
  return db;
}


Example Example::ConvDiff(const ParameterDatabase& param_db)
{
  int example = Example::get_example_from_database(param_db);
  auto db = get_default_ConvDiff_Example_parameters(example);
  db.merge(param_db, false);
  double diffusion = db["diffusion_coefficient"];

  // indicate if two or three space dimensions
  bool two_d = true;
#ifdef __3D__
  two_d = false;
#endif

  // function needed to create the PDECoefficients
  std::function<void(const Point&, double, std::vector<double>&)> f;
  // vectors which will contain exactly one entry each.
  std::vector<BoundaryCondition> bc;     // boundary condition
  std::vector<BoundaryData> bd;          // boundary data
  std::vector<AnalyticalFunction> exact; // exact solution (possibly unknown)
  std::string name;                      // to print some information later on

  switch(example)
  {
    case 0: // sine-cosine laplace example
    {
      name = "sine laplace";
      if(two_d)
      {
        // boundary condition:
        auto c = [](unsigned int comp)
        {
          return comp == 1 ? NEUMANN : DIRICHLET;
        };
        bc.push_back(BoundaryCondition(c));
        // exact solution
        auto e = [](const Point& point, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          constexpr double p = M_PI; // 2*M_PI;
          // isnan(point.z()) should be true, because it is a 2D example
          v.set(std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D00);
          v.set(p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D10);
          v.set(p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D01);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D20);
          v.set(p * p * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D11);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(e));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc[0], exact[0]));
        auto g = [diffusion](unsigned int component, double t, double)
        {
          return component == 1 ? -diffusion * M_PI * std::sin(M_PI * t) : 0.;
        };
        bd.push_back(BoundaryData(g));
        // the coefficient function
        f = [diffusion](const Point& point, double, std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          constexpr double p = M_PI; // 2*M_PI;
          // isnan(point.z()) should be true, because it is a 2D example
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = 0;         // convection in x-direction
          coeffs[2] = 0;         // convection in y-direction
          coeffs[3] = 0;         // reaction
          coeffs[4] = diffusion * 2 * p * p * std::sin(p * x) * std::sin(p * y); // rhs
        };
      }
      else // 3D
      {
        // boundary condition:
        bc.push_back(BoundaryCondition(DIRICHLET)); // all Dirichlet
        // exact solution
        auto e = [](const Point& point, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          const double z = point.z;
          constexpr double p = M_PI; // 2*M_PI;
          v.set(std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0, MultiIndex::D000);
          v.set(p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0, MultiIndex::D100);
          v.set(p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0, MultiIndex::D010);
          v.set(p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0, MultiIndex::D001);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D200);
          v.set(p * p * std::cos(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D110);
          v.set(p * p * std::cos(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D101);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D020);
          v.set(p * p * std::sin(p * x) * std::cos(p * y) * std::cos(p * z), 0,
                MultiIndex::D011);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D002);
        };
        exact.push_back(AnalyticalFunction(e));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc[0], exact[0]));
        bd.push_back(BoundaryData(0.0)); // all zero boundary data
        // the coefficient function
        f = [diffusion](const Point& point, double, std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          const double z = point.z;
          constexpr double p = M_PI; // 2*M_PI;
          coeffs[0] = diffusion;     // diffusion
          coeffs[1] = 0;             // convection in x-direction
          coeffs[2] = 0;             // convection in y-direction
          coeffs[3] = 0;             // convection in z-direction
          coeffs[4] = 0;             // reaction
          coeffs[5] = 3 * p * p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z); // rhs
        };
      }
      break;
    }
    case 1: // two interior layers
    {
      name = "two interior layers";
      if(two_d)
      {
        // boundary condition:
        auto c = [](unsigned int comp)
        {
          return comp == 3 ? NEUMANN : DIRICHLET;
        };
        bc.push_back(BoundaryCondition(c));
        // exact solution unknown
        exact.push_back(AnalyticalFunction());
        // boundary data
        auto g = [](unsigned int component, double t, double)
        {
          return (component == 0 && t > 1. / 3. && t < 2. / 3.) ? 1. : 0.;
        };
        bd.push_back(BoundaryData(g));
        // the coefficient function
        f = [diffusion](const parmoon::Point& point, double, std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = -y;        // convection in x-direction
          coeffs[2] = x;         // convection in y-direction
          coeffs[3] = 0.;        // reaction
          coeffs[4] = 0.;        // rhs
          coeffs[5] = std::sqrt(x * x + y * y);
        };
      }
      else // 3D
      {
        ErrThrow("Example Two interior layers: This is a 2D example");
      }
      break;
    }
    case 2: // Hemker
    {
      name = "Hemker";
      // boundary condition:
      auto c = [](unsigned int comp)
      {
        // comp==3 for inflow boundary, comp==4 for cylinder
        return (comp == 3 || comp == 4) ? DIRICHLET : NEUMANN;
      };
      bc.push_back(BoundaryCondition(c));
      // exact solution (unknown)
      exact.push_back(AnalyticalFunction());
      auto g = [](unsigned int component, double, double)
      {
        return component == 4 ? 1. : 0.;
      };
      bd.push_back(BoundaryData(g));
      if(two_d)
      {
        // the coefficient function
        f = [diffusion](const parmoon::Point&, double, std::vector<double>& coeffs)
        {
          const double angle = 0.;
          const double v1 = std::cos(angle);
          const double v2 = std::sin(angle);
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = v1;        // convection in x-direction
          coeffs[2] = v2;        // convection in y-direction
          coeffs[3] = 0.;        // reaction
          coeffs[4] = 0.;        // rhs
          coeffs[5] = std::sqrt(v1 * v1 + v2 * v2);
        };
      }
      else // 3D
      {
        // the coefficient function
        f = [diffusion](const Point&, double, std::vector<double>& coeffs)
        {
          const double angle = 0.;
          const double v1 = std::cos(angle);
          const double v2 = std::sin(angle);
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = v1;        // convection in x-direction
          coeffs[2] = v2;        // convection in y-direction
          coeffs[3] = 0.;        // convection in z-direction
          coeffs[4] = 0.;        // reaction
          coeffs[5] = 0.;        // rhs
        };
      }
      break;
    }
    case 3: // sharp boundary layer
    {
      name = "sharp boundary layer";
      if(two_d)
      {
        // boundary condition:
        auto c = [](unsigned int comp)
        {
          return comp == 2 ? NEUMANN : DIRICHLET;
        };
        bc.push_back(BoundaryCondition(c));
        // exact solution unknown
        exact.push_back(AnalyticalFunction());
        // boundary data
        auto g = [](unsigned int component, double t, double)
        {
          return (component == 3 && t < 0.5) ? 1. : 0.;
        };
        bd.push_back(BoundaryData(g));
        // the coefficient function
        f = [diffusion](const Point&, double, std::vector<double>& coeffs)
        {
          coeffs[0] = diffusion;            // diffusion
          coeffs[1] = std::cos(M_PI / 18.); // convection in x-direction
          coeffs[2] = std::sin(M_PI / 18.); // convection in y-direction
          coeffs[3] = 0.;                   // reaction
          coeffs[4] = 0.;                   // rhs
          coeffs[5] = std::sqrt(coeffs[1] * coeffs[1] + coeffs[2] * coeffs[2]);
        };
      }
      else // 3D
      {
        ErrThrow("Example Two interior layers: This is a 2D example");
      }
      break;
    }
    default:
    {
      ErrThrow("unknown example ", example);
      break;
    }
  }

  PDECoefficients coeffs(f, Problem_type::ConvDiffReac, false, false);

  Output::info<2>("Example", "using the convection-diffusion example: ", name);
  Output::dash<2>("The diffusion constant is ", diffusion);
  return Example(std::move(bc), std::move(bd), std::move(coeffs),
                 std::move(exact));
}
