#include <cmath> // e.g. sin
#include <Example_NonStationary.h>
#include <MooNMD_Io.h>
using namespace parmoon;

Example_NonStationary Example_NonStationary::Time_NavierStokes(
  const ParameterDatabase& param_db)
{
  int example = Example::get_example_from_database(param_db);
  double reynolds = param_db["reynolds_number"];
  bool only_stokes = false; // some examples also provide a linear version
  // flag, false means the laplace term is discretized with the gradient, true
  // means it is discretized with the deformation tensor
  //bool deformation_tensor = false;
  double initial_time = 0.0; // should this be taken from the database?

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
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));

      constexpr double p = M_PI; // 2*M_PI
      // exact solution, first velocity component
      if(two_d)
      {
        auto u1 = [](const Point& point, double t, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          v.set(std::sin(t) * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D00);
          v.set(p * std::sin(t) * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D10);
          v.set(p * std::sin(t) * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D01);
          v.set(-p * p * std::sin(t) * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D20);
          v.set(p * p * std::sin(t) * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D11);
          v.set(-p * p * std::sin(t) * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(u1));
        // exact solution, second velocity component
        auto u2 = [](const Point& point, double t, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          v.set(std::sin(t) * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D00);
          v.set(-p * std::sin(t) * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D10);
          v.set(-p * std::sin(t) * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D01);
          v.set(-p * p * std::sin(t) * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D20);
          v.set(p * p * std::sin(t) * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D11);
          v.set(-p * p * std::sin(t) * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(u2));
        // exact solution, pressure component
        auto pr = [](const Point& point, double t, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          v.set(std::sin(t) * (std::sin(p * x) + std::cos(p * y) - 2 / p), 0, MultiIndex::D00);
          v.set(p * std::sin(t) * std::cos(p * x), 0, MultiIndex::D10);
          v.set(-p * std::sin(t) * std::sin(p * y), 0, MultiIndex::D01);
          v.set(-p * p * std::sin(t) * std::sin(p * x), 0, MultiIndex::D20);
          v.set(0., 0, MultiIndex::D11);
          v.set(-p * p * std::sin(t) * std::cos(p * y), 0, MultiIndex::D02);
        };
        exact.push_back(AnalyticalFunction(pr));
        // initial condition
        initial.push_back(AnalyticalFunction(exact[0], initial_time));
        initial.push_back(AnalyticalFunction(exact[1], initial_time));
        initial.push_back(AnalyticalFunction(exact[2], initial_time));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc, exact));
        // Boundary data for first velocity component ((all zero))
        bd.push_back(BoundaryData(0.0));
        // Boundary data for second velocity component
        auto g2 = [](unsigned int component, double t, double time)
        {
          switch(component)
          {
            case 0:
              return std::sin(time) * std::cos(p * t);
              break;
            case 1:
              return -std::sin(time) * std::cos(p * t);
              break;
            case 2:
              return -std::sin(time) * std::cos(p * (1 - t));
              break;
            case 3:
              return std::sin(time) * std::cos(M_PI * (1 - t));
              break;
            default:
              ErrThrow("unknown boundary component");
          }
        };
        bd.push_back(BoundaryData(g2));
        bd.push_back(BoundaryData(0.0)); // all zero for pressure
        // the coefficient function
        f = [reynolds, only_stokes](const Point& point, double time,
                                       std::vector<double>& coeffs)
        {
          const double x = point.x;
          const double y = point.y;
          // isnan(point.z()) should be true, because it is a 2D example
          coeffs[0] = 1. / reynolds; // reynolds number
          // right hand side, first component
          coeffs[1] = std::cos(time) * std::sin(p * x) * std::sin(p * y); // time derivative
          // -(1/reynolds)*laplace(u1)
          coeffs[1] += 2 * p * p * std::sin(time) * std::sin(p * x) * std::sin(p * y)
                       / reynolds;
          coeffs[1] += p * std::sin(time) * std::cos(p * x); // x-derivative of pressure
          // right hand side, second component
          coeffs[2] = std::cos(time) * std::cos(p * x) * std::cos(p * y); // time derivative
          // -(1/reynolds)*laplace(u2)
          coeffs[2] += 2 * p * p * std::sin(time) * std::cos(p * x) * std::cos(p * y)
                       / reynolds;
          coeffs[2] += -p * std::sin(time) * std::sin(p * y); // y-derivative of pressure
          if(!only_stokes)
          {
            // adding nonlinear term
            coeffs[1] += std::sin(time) * std::sin(p * x) * std::sin(p * y) * p * std::sin(time)
                         * std::cos(p * x) * std::sin(p * y);
            coeffs[1] += std::sin(time) * std::cos(p * x) * std::cos(p * y) * p * std::sin(time)
                         * std::sin(p * x) * std::cos(p * y);

            coeffs[2] += std::sin(time) * std::sin(p * x) * std::sin(p * y) * (-p) * std::sin(time)
                         * std::sin(p * x) * std::cos(p * y);
            coeffs[2] += std::sin(time) * std::cos(p * x) * std::cos(p * y) * (-p) * std::sin(time)
                         * std::cos(p * x) * std::sin(p * y);
          }
          coeffs[3] = 0.; // divergence
        };
      }
      else // 3D
      {
        /// @todo implement sine example for TNSE3D, this one is time indepenent
        ErrThrow("implement a sine example for TNSE3D!");
        auto u1 = [](const Point& point, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          const double z = point.z;
          v.set(std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0, MultiIndex::D000);
          v.set(-p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0, MultiIndex::D100);
          v.set(p * std::cos(p * x) * std::cos(p * y) * std::sin(p * z), 0, MultiIndex::D010);
          v.set(p * std::cos(p * x) * std::sin(p * y) * std::cos(p * z), 0, MultiIndex::D001);
          v.set(-p * p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D200);
          v.set(-p * p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D110);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D101);
          v.set(-p * p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D020);
          v.set(p * p * std::cos(p * x) * std::cos(p * y) * std::cos(p * z), 0,
                MultiIndex::D011);
          v.set(-p * p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D002);
        };
        exact.push_back(AnalyticalFunction(u1));
        // exact solution, second velocity component
        auto u2 = [](const Point& point, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          const double z = point.z;
          v.set(std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0, MultiIndex::D000);
          v.set(p * std::cos(p * x) * std::cos(p * y) * std::sin(p * z), 0, MultiIndex::D100);
          v.set(-p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0, MultiIndex::D010);
          v.set(p * std::sin(p * x) * std::cos(p * y) * std::cos(p * z), 0, MultiIndex::D001);
          v.set(-p * p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D200);
          v.set(-p * p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D110);
          v.set(p * p * std::cos(p * x) * std::cos(p * y) * std::cos(p * z), 0,
                MultiIndex::D101);
          v.set(-p * p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D020);
          v.set(-p * p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D011);
          v.set(-p * p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D002);
        };
        exact.push_back(AnalyticalFunction(u2));
        // exact solution, third velocity component
        auto u3 = [](const Point& point, FunctionEvaluation& v)
        {
          const double x = point.x;
          const double y = point.y;
          const double z = point.z;
          v.set(-2 * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0, MultiIndex::D000);
          v.set(-2 * p * std::cos(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D100);
          v.set(-2 * p * std::sin(p * x) * std::cos(p * y) * std::cos(p * z), 0,
                MultiIndex::D010);
          v.set(2 * p * std::sin(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D001);
          v.set(2 * p * p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D200);
          v.set(-2 * p * p * std::cos(p * x) * std::cos(p * y) * std::cos(p * z), 0,
                MultiIndex::D110);
          v.set(2 * p * p * std::cos(p * x) * std::sin(p * y) * std::sin(p * z), 0,
                MultiIndex::D101);
          v.set(2 * p * p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D020);
          v.set(2 * p * p * std::sin(p * x) * std::cos(p * y) * std::sin(p * z), 0,
                MultiIndex::D011);
          v.set(2 * p * p * std::sin(p * x) * std::sin(p * y) * std::cos(p * z), 0,
                MultiIndex::D002);
        };
        exact.push_back(AnalyticalFunction(u3));
        // exact solution, pressure component
        auto pr = [](const Point&, FunctionEvaluation& v)
        {
          v.reset(); // all zero
        };
        exact.push_back(AnalyticalFunction(pr));
        // initial condition
        initial.push_back(AnalyticalFunction(exact[0], initial_time));
        initial.push_back(AnalyticalFunction(exact[1], initial_time));
        initial.push_back(AnalyticalFunction(exact[2], initial_time));
        initial.push_back(AnalyticalFunction(exact[3], initial_time));
        // boundary data, automatically computed from known solution
        // bd.push_back(BoundaryData(bc, exact));
        // Boundary data for first velocity component
        auto g1 = [](const Point& point, double)
        {
          return std::cos(p * point.x) * std::sin(p * point.y) * std::sin(p * point.z);
        };
        bd.push_back(BoundaryData(g1));
        // Boundary data for first velocity component
        auto g2 = [](const Point& point, double)
        {
          return std::sin(p * point.x) * std::cos(p * point.y) * std::sin(p * point.z);
        };
        bd.push_back(BoundaryData(g2));
        // Boundary data for third velocity component
        auto g3 = [](const Point& point, double)
        {
          return -2 * std::sin(p * point.x) * std::sin(p * point.y) * std::cos(p * point.z);
        };
        bd.push_back(BoundaryData(g3));
        bd.push_back(BoundaryData(0.0)); // all zero for pressure
        // the coefficient function
        f = [reynolds, only_stokes, u1, u2, u3, pr](
            const Point& point, double, std::vector<double>& coeffs)
        {
          constexpr MultiIndex D000 = MultiIndex::D000;
          constexpr MultiIndex D100 = MultiIndex::D100;
          constexpr MultiIndex D010 = MultiIndex::D010;
          constexpr MultiIndex D001 = MultiIndex::D001;
          FunctionEvaluation u1_values(3, 1);
          FunctionEvaluation u2_values(3, 1);
          FunctionEvaluation u3_values(3, 1);
          FunctionEvaluation pr_values(3, 1);
          u1(point, u1_values);
          u2(point, u2_values);
          u3(point, u3_values);
          pr(point, pr_values);
          coeffs[0] = 1. / reynolds; // reynolds number
          // Stokes: diffusion and pressure gradient
          // right hand side, first component
          coeffs[1] = -(u1_values(0, MultiIndex::D200)
                        + u1_values(0, MultiIndex::D020)
                        + u1_values(0, MultiIndex::D002)) / reynolds
                      + pr_values(0, D100);
          // right hand side, second component
          coeffs[2] = -(u2_values(0, MultiIndex::D200)
                        + u2_values(0, MultiIndex::D020)
                        + u2_values(0, MultiIndex::D002)) / reynolds
                      + pr_values(0, D010);
          // right hand side, third component
          coeffs[3] = -(u3_values(0, MultiIndex::D200)
                        + u3_values(0, MultiIndex::D020)
                        + u3_values(0, MultiIndex::D002)) / reynolds
                      + pr_values(0, D001);
          if(!only_stokes)
          {
            // adding nonlinear term
            coeffs[1] += u1_values(0, D000) * u1_values(0, D100)
                         + u2_values(0, D000) * u1_values(0, D010)
                         + u3_values(0, D000) * u1_values(0, D001);
            coeffs[2] += u1_values(0, D000) * u2_values(0, D100)
                         + u2_values(0, D000) * u2_values(0, D010)
                         + u3_values(0, D000) * u2_values(0, D001);
            coeffs[3] += u1_values(0, D000) * u3_values(0, D100)
                         + u2_values(0, D000) * u3_values(0, D010)
                         + u3_values(0, D000) * u3_values(0, D001);
          }
          coeffs[4] = 0.; // divergence
        };
      }
      break;
    }
    case 1: // Poisseuille flow
    {
      ErrThrow("Time dependent Poisseuille flow example not yet implemented");
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));
      auto u1 = [](const Point& point, FunctionEvaluation& v)
      {
        const double y = point.y;
        v.set(4 * y * (1 - y), 0, MultiIndex::D00);
        v.set(0., 0, MultiIndex::D10);
        v.set(4 - 8 * y, 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(-8., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component (all zero)
      exact.push_back(AnalyticalFunction(0.0));
      // exact solution, pressure component
      auto pr = [reynolds](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        v.set(8 * reynolds * (0.5 - x), 0, MultiIndex::D00);
        v.set(-8 * reynolds, 0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(pr));
      // initial condition
      initial.push_back(AnalyticalFunction(exact[0], initial_time));
      initial.push_back(AnalyticalFunction(exact[1], initial_time));
      initial.push_back(AnalyticalFunction(exact[2], initial_time));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc, exact));
      // Boundary data for first velocity component
      auto g1 = [](unsigned int component, double t, double)
      {
        return (component == 1 || component == 3) ? 4 * t * (1 - t) : 0.;
      };
      bd.push_back(BoundaryData(g1));
      // Boundary data for first velocity component (all zero)
      bd.push_back(BoundaryData(0.0));
      bd.push_back(BoundaryData(0.0)); // all zero for pressure
      // the coefficient function
      f = [reynolds](const Point&, double,
                     std::vector<double>& coeffs)
      {
        coeffs[0] = 1. / reynolds; // reynolds number
        coeffs[1] = 0.;            // right hand side, first component
        coeffs[2] = 0.;            // right hand side, second component
        coeffs[3] = 0.;            // divergence
      };
      break;
    }
    case 2: // driven cavity
    {
      ErrThrow("Time dependent driven cavity example not yet implemented");
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));

      // exact solution, first velocity component (unknown)
      exact.push_back(AnalyticalFunction());
      // exact solution, second velocity component (unknown)
      exact.push_back(AnalyticalFunction());
      if(!two_d)
      {
        // exact solution, third velocity component (unknown)
        exact.push_back(AnalyticalFunction());
      }
      // exact solution, pressure component (unknown)
      exact.push_back(AnalyticalFunction());
      // initial condition
      initial.push_back(AnalyticalFunction(exact[0], initial_time));
      initial.push_back(AnalyticalFunction(exact[1], initial_time));
      initial.push_back(AnalyticalFunction(exact[2], initial_time));
      if(!two_d)
      {
        initial.push_back(AnalyticalFunction(exact[3], initial_time));
      }
      // boundary data
      // Boundary data for first velocity component
      if(two_d)
      {
        auto g1 = [](unsigned int component, double t, double)
        {
          return (component == 2 && t > 0.00001 && t < 0.99999) ? 1. : 0.;
        };
        bd.push_back(BoundaryData(g1));
      }
      else
      {
        // 3D
        auto g1 = [](const Point& point, double)
        {
          return (1. - point.z < 1e-10) ? 1. : 0.;
        };
        bd.push_back(BoundaryData(g1));
      }
      // Boundary data for second velocity component
      bd.push_back(BoundaryData(0.));
      if(!two_d)
      {
        // Boundary data for third velocity component
        bd.push_back(BoundaryData(0.));
      }
      bd.push_back(BoundaryData(0.0)); // all zero for pressure
      // the coefficient function
      f = [reynolds](const Point& point, double, std::vector<double>& coeffs)
      {
        coeffs[0] = 1. / reynolds; // reynolds number
        coeffs[1] = 0.;            // right hand side, first component
        coeffs[2] = 0.;            // right hand side, second component
        coeffs[3] = 0.; // divergence (2D) / rhs, third component (3D)
        if(point.dimension() == 3)
          coeffs[3] = 0.; // divergence (3D)
      };
      break;
    }
    case 3:
    {
      ErrThrow("NSE Example 3: to be implemented ");
      break;
    }
    default:
      ErrThrow("unknown example ", example);
      break;
  }

  PDECoefficients coeffs(f, Problem_type::ConvDiffReac);

  Output::info<2>("Example", "using the Navier--Stokes example: ", name);
  Output::dash<2>("The Reynolds number is ", reynolds);
  return Example_NonStationary(std::move(bc), std::move(bd), std::move(coeffs),
                               std::move(exact), std::move(initial));
}
