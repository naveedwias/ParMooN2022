#include <cmath> // e.g. sin
#include <Example_NonStationary.h>
#include <MooNMD_Io.h>
using namespace parmoon;

Example_NonStationary Example_NonStationary::Time_ConvDiff(
  const ParameterDatabase& param_db)
{
  int example = Example::get_example_from_database(param_db);
  double diffusion = param_db["diffusion_coefficient"];
  double initial_time = 0.0; // should this be taken from the database?

  bool lhs_depends_on_time = true;
  bool rhs_depends_on_time = true;
  // indicate if two or three space dimensions
  bool two_d = true;
#ifdef __3D__
  two_d = false;
#endif

  // function needed to create the PDECoefficients
  std::function<void(const Point&, double, std::vector<double>&)> f;
  // vectors which will contain exactly one entry each.
  std::vector<BoundaryCondition> bc;       // boundary condition
  std::vector<BoundaryData> bd;            // boundary data
  std::vector<AnalyticalFunction> exact;   // exact solution (possibly unknown)
  std::vector<AnalyticalFunction> initial; // initial solution
  std::string name;                        // to print some information later on

  switch(example)
  {
    case 0: // sine example
    {
      name = "sine convection-diffusion-reaction with constant coefficients";
      if(two_d)
      {
        // boundary condition:
        bc.push_back(BoundaryCondition(DIRICHLET));
        constexpr double p = 2 * M_PI;
        // exact solution
        auto e = [](const Point& point, double t, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          v.set(std::sin(t) * (std::sin(p * x) * std::sin(p * y) + 1), 0, MultiIndex::D00);
          v.set(std::sin(t) * p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D10);
          v.set(std::sin(t) * p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D01);
          v.set(-std::sin(t) * p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D20);
          v.set(std::sin(t) * p * p * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D11);
          v.set(-std::sin(t) * p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(e));
        // initial condition
        initial.push_back(AnalyticalFunction(exact.back(), initial_time));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc[0], exact[0]));
        auto g = [](const Point& point, double time)
        {
          return std::sin(time) * (std::sin(p * point.x) * std::sin(p * point.y) + 1);
        };
        bd.push_back(BoundaryData(g));
        // the coefficient function
        lhs_depends_on_time = false;
        rhs_depends_on_time = true;
        f = [diffusion](const Point& point, double time,
                           std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = 1.;        // convection in x-direction
          coeffs[2] = 2.;        // convection in y-direction
          coeffs[3] = 1.;        // reaction
          coeffs[4] = std::cos(time)
                      * (std::sin(p * x) * std::sin(p * y) + 1); // time derivative
          coeffs[4] += diffusion * 2 * std::sin(time) * p * p * std::sin(p * x)
                       * std::sin(p * y);
          coeffs[4] += coeffs[1] * std::sin(time) * p * std::cos(p * x) * std::sin(p * y);
          coeffs[4] += coeffs[2] * std::sin(time) * p * std::sin(p * x) * std::cos(p * y);
          coeffs[4] += coeffs[3] * std::sin(time) * (std::sin(p * x) * std::sin(p * y) + 1);
        };
      }
      else // 3D
      {
        ErrThrow("A 3D time dependent convection-diffusion problem is not yet"
                 "implemented");
      }
      break;
    }
    case 1: // linear (in space and time) solution
    {
      name = "linear convection-diffusion-reaction with constant coefficients";
      if(two_d)
      {
        // boundary condition:
        bc.push_back(BoundaryCondition(DIRICHLET));
        // exact solution
        auto e = [](const Point& point, double t, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          v.set(1+2*x+3*t*y, 0, MultiIndex::D00);
          v.set(2., 0, MultiIndex::D10);
          v.set(3*t, 0, MultiIndex::D01);
          v.set(0., 0, MultiIndex::D20);
          v.set(0., 0, MultiIndex::D11);
          v.set(0., 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(e));
        // initial condition
        initial.push_back(AnalyticalFunction(exact.back(), initial_time));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc[0], exact[0]));
        auto g = [](unsigned int component, double t, double time)
        {
          switch(component)
          {
            case 0: return 1*2*t; break;
            case 1: return 3+3*t*time; break;
            case 2: return 1+3*time+2*(1-t); break;
            case 3: return 1+3*time*(1-t); break;
            default: ErrThrow("unknown boundary component ", component);
          }
        };
        bd.push_back(BoundaryData(g));
        // the coefficient function
        lhs_depends_on_time = false;
        rhs_depends_on_time = true;
        f = [diffusion](const Point& point, double time,
                           std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          coeffs[0] = diffusion; // diffusion
          coeffs[1] = 1.;        // convection in x-direction
          coeffs[2] = 2.;        // convection in y-direction
          coeffs[3] = 1.;        // reaction
          coeffs[4] = 3*y; // time derivative
          // Laplace(u) is zero here
          coeffs[4] += coeffs[1] * 2. + coeffs[2] * 3*time; // convection
          coeffs[4] += coeffs[3] * (1 + 2*x + 3*time*y); // reaction
        };
      }
      else // 3D
      {
        ErrThrow("A 3D time dependent convection-diffusion problem is not yet"
                 "implemented");
      }
      break;
    }
    case 3: // rotating bodies
    {
      name = "rotating bodies";
      if(!two_d)
      {
        ErrThrow("the rotating bodies example is a 2D example");
      }
      Output::info<>("Rotating bodies example",
                     "the exact function is known, "
                     "however not its derivatives. In fact it is not in H^1.");
      if(diffusion != 1.0e-20)
      {
        diffusion = 1.0e-20;
        Output::dash<>("diffusion constant set to ", diffusion);
      }
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));
      // exact solution
      auto e = [](const Point& point, double t, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        // rotation
        double x0 = std::cos(t) * x + std::sin(t) * y - 0.5 * std::cos(t) - 0.5 * std::sin(t) + 0.5;
        double y0 = -std::sin(t) * x + std::cos(t) * y + 0.5 * std::sin(t) - 0.5 * std::cos(t)
                    + 0.5;

        double val = 0; // to be set to v
        // cylinder
        if((x0 - 0.5) * (x0 - 0.5) + (y0 - 0.75) * (y0 - 0.75) <= 0.15 * 0.15)
        {
          if(std::abs(x0 - 0.5) >= 0.0225 || y0 >= 0.85)
            val = 1.;
        }

        // conical body
        double conical = 1.0 / 0.15 * std::sqrt((x0 - 0.5) * (x0 - 0.5)
                                           + (y0 - 0.25) * (y0 - 0.25));
        if(conical <= 1.0)
          val = 1. - conical;

        // hump
        double hump = 1.0 / 0.15 * std::sqrt((x0 - 0.25) * (x0 - 0.25)
                                        + (y0 - 0.5) * (y0 - 0.5));
        if(hump <= 1.0)
          val = 0.25 * (1.0 + (std::cos(M_PI * hump)));

        v.set(val, 0, MultiIndex::D00);
        // derivatives are not known (not in L^2)
        v.set(0., 0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(e));
      // initial condition
      initial.push_back(AnalyticalFunction(exact.back(), initial_time));
      // boundary data
      bd.push_back(BoundaryData(0.0));
      // the coefficient function
      lhs_depends_on_time = false;
      rhs_depends_on_time = false;
      f = [diffusion](const Point& point, double, std::vector<double>& coeffs)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        coeffs[0] = diffusion; // diffusion
        coeffs[1] = 0.5 - y;   // convection in x-direction
        coeffs[2] = x - 0.5;   // convection in y-direction
        coeffs[3] = 0.;        // reaction
        coeffs[4] = 0.;        // rhs
      };
      break;
    }
    default:
    {
      ErrThrow("unknown example ", example);
      break;
    }
  }

  PDECoefficients coeffs(f, Problem_type::ConvDiffReac, lhs_depends_on_time, 
                         rhs_depends_on_time);

  Output::info<2>("Example", "using the time dependent convection-diffusion "
                             "example: ",
                  name);
  Output::dash<2>("The diffusion constant is ", diffusion);
  return Example_NonStationary(std::move(bc), std::move(bd), std::move(coeffs),
                               std::move(exact), std::move(initial));
}

