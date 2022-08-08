#include <cmath> // e.g. sin
#include <Example.h>
#include <MooNMD_Io.h>

using namespace parmoon;

ParameterDatabase get_default_Darcy_Example_parameters(int example)
{
  ParameterDatabase db("default example database");
  db.add("example", example,
         "Choose which example to run. Depending on the example other values "
         "parameters are also taken into acount. "
         "0 -- with known sine and cosine solution (e.g. to study convergence "
         "rates),  "
         "1 -- simple benchmark with known solution (linear velocity, zero "
         "pressure),  "
         "2 -- simple benchmark with known solution (cubic velocity, quartic "
         "pressure),  "
         "3 -- obstacle with lower permeability. Set 'obstacle_type' for "
         "various obstacles, and 'permeability_factor' to determine how much "
         "lower its permeability should be."
         "4 -- 5-spot problem, see Arif Masud, Thomas J.R. Hughes: 'A "
         "stabilized mixed finite element method for Darcy flow', chapter 4.2",
         -5, 200);
  db.add("permeability_coefficient", 1.0,
         "Abreviated: K. Its inverse is the (scalar, constant) coefficient in "
         "front of the velocity term in the Darcy problem. Typically K is "
         "small. K^{-1}u + grad p = f", 0., 100.);
  
  switch(example)
  {
    case 3: // obstacle example
    {
      db.add("obstacle_type", 0u, "choose different obstacles:"
             "0 - circle in center, 1 - rectangle in center, "
             "2 - diamond in center, "
             "3 - two boxes, one at the top, one at the bottom, "
             "4 - two half circles, one at the top, one at the bottom, "
             "5 - two boxes, one at the top, one at the bottom, boxes have an "
             "offset, "
             "6 - two half ellipsoids with an offset, one at the top, one at "
             "the bottom, "
             "7 - continuous change in permeability in form of a hump",
             0u, 7u);
      db.add("permeability_factor", 0.1, 
             "The permeability within the obstacle is multiplied by this "
             "factor. So a factor of 1 makes the entire domain uniform. "
             "Smaller factors force the flow around the obstacle to some "
             "extent. Larger factors are possible, but then the 'obstacle' "
             "really becomes a preferred flow path.");
      break;
    }
    case 4:
    {
      db.add("permeability_factor", 0.1,
             "The permeability within the zones of the checkerboard domain is "
             "multiplied by this factor. So a factor of 1 makes the entire "
             "domain uniform. Smaller factors force the flow away from the "
             "diagonal to some extent. Larger factors are possible, too, and "
             "focus the flow towards the diagonal.");
    }
  }
  return db;
}

Example Example::Darcy(const ParameterDatabase& param_db)
{
  int example = Example::get_example_from_database(param_db);
  auto db = get_default_Darcy_Example_parameters(example);
  db.merge(param_db, false);
  double permeability = db["permeability_coefficient"];

#ifdef __3D__
  ErrThrow("Darcy is not yet ported to 3D from MooNMD");
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
    case 0: // sine-cosine example
    {
      // boundary condition:
      auto c = [](unsigned int comp)
      {
        return (comp == 0 || comp == 3) ? NEUMANN : DIRICHLET;
      };
      bc.push_back(BoundaryCondition(c));
      // exact solution, first velocity component
      auto u1 = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        constexpr double p = 2 * M_PI;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set(-p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D00);
        v.set(p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D10);
        v.set(-p * p * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D01);
        v.set(p * p * p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D20);
        v.set(p * p * p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D11);
        v.set(p * p * p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component
      auto u2 = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        constexpr double p = 2 * M_PI;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set(-p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D00);
        v.set(-p * p * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D10);
        v.set(p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D01);
        v.set(p * p * p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D20);
        v.set(p * p * p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D11);
        v.set(p * p * p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u2));
      // exact solution, pressure component
      auto p = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        constexpr double p = 2 * M_PI;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set(std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D00);
        v.set(p * std::cos(p * x) * std::sin(p * y), 0, MultiIndex::D10);
        v.set(p * std::sin(p * x) * std::cos(p * y), 0, MultiIndex::D01);
        v.set(-p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D20);
        v.set(p * p * std::cos(p * x) * std::cos(p * y), 0, MultiIndex::D11);
        v.set(-p * p * std::sin(p * x) * std::sin(p * y), 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(p));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc[0], exact));
      auto g = [](unsigned int component, double t, double)
      {
        if(component == 0 || component == 3) // Neumann
          return 0.0;
        // else: Dirichlet
        if(component == 1)
          return -2 * M_PI * std::sin(2 * M_PI * t);
        // else if(component == 2)
        return -2 * M_PI * std::sin(2 * M_PI * (1 - t));
      };
      bd.push_back(BoundaryData(g));
      // the coefficient function
      f = [permeability](const Point& point, double, 
                         std::vector<double>& coeffs)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        constexpr double p = 2 * M_PI;
        coeffs[0] = 1./permeability; // permeability
        coeffs[1] = 0.;              // right hand side, first component
        coeffs[2] = 0.;              // right hand side, second component
        coeffs[3] = 2 * p * p * std::sin(p * x) * std::sin(p * y); // divergence
      };
      break;
    }
    case 1: // simple benchmark, zero pressure linear velocity
    {
      // boundary condition:
      auto c = [](unsigned int comp)
      {
        return (comp == 0) ? NEUMANN : DIRICHLET;
      };
      bc.push_back(BoundaryCondition(c));
      // exact solution, first velocity component
      auto u1 = [](const Point& point, FunctionEvaluation& v)
      {
        v.set(point.y, 0, MultiIndex::D00);
        v.set(0., 0, MultiIndex::D10);
        v.set(1., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component
      auto u2 = [](const Point& point, FunctionEvaluation& v)
      {
        v.set(point.x, 0, MultiIndex::D00);
        v.set(1., 0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u2));
      // exact solution, pressure component
      auto p = [](const Point&, FunctionEvaluation& v)
      {
        v.set(0., 0, MultiIndex::D00);
        v.set(0., 0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(p));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc[0], exact));
      auto g = [](unsigned int component, double t, double)
      {
        if(component == 0) // Neumann, (x,y) = (t,0), (nx,ny) = (0,-1)
          return -t;
        if(component == 1) // Dirichlet, (x,y) = (1,t), (nx,ny) = (1,0)
          return t;
        if(component == 2) // Dirichlet, (x,y) = (1-t,1), (nx,ny) = (0,1)
          return 1 - t;
        if(component == 3) // Dirichlet, (x,y) = (0,1-t), (nx,ny) = (-1,0)
          return t - 1;
        return std::nan("1"); // only to avoid compiler warning
      };
      bd.push_back(BoundaryData(g));
      // the coefficient function
      f = [permeability](const Point& point, double, 
                         std::vector<double>& coeffs)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        constexpr double p = 2 * M_PI;
        coeffs[0] = 1./permeability; // permeability
        coeffs[1] = 0.;              // right hand side, first component
        coeffs[2] = 0.;              // right hand side, second component
        coeffs[3] = 2 * p * p * std::sin(p * x) * std::sin(p * y); // divergence
      };
      break;
    }
    case 2: // polynomial solution with cubic pressure
    {
      // boundary condition:
      auto c = [](unsigned int comp)
      {
        return (comp == 3) ? DIRICHLET : NEUMANN;
      };
      bc.push_back(BoundaryCondition(c));
      // exact solution, first velocity component
      auto u1 = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set(-x * x * y + y * y * y / 3.0, 0, MultiIndex::D00);
        v.set(-2 * x * y, 0, MultiIndex::D10);
        v.set(-x * x + y * y, 0, MultiIndex::D01);
        v.set(-2 * y, 0, MultiIndex::D20);
        v.set(-2 * x, 0, MultiIndex::D11);
        v.set(2 * y, 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component
      auto u2 = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set(-x * x * x / 3.0 + x * y * y, 0, MultiIndex::D00);
        v.set(-x * x + y * y, 0, MultiIndex::D10);
        v.set(2 * x * y, 0, MultiIndex::D01);
        v.set(-2 * x, 0, MultiIndex::D20);
        v.set(2 * y, 0, MultiIndex::D11);
        v.set(2 * x, 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u2));
      // exact solution, pressure component
      auto p = [](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        v.set((x * x * x * y - x * y * y * y) / 3.0, 0, MultiIndex::D00);
        v.set(x * x * y - y * y * y / 3.0, 0, MultiIndex::D10);
        v.set(x * x * x / 3.0 - x * y * y, 0, MultiIndex::D01);
        v.set(2 * x * y, 0, MultiIndex::D20);
        v.set(x * x - y * y, 0, MultiIndex::D11);
        v.set(-2 * x * y, 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(p));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc[0], exact));
      auto g = [](unsigned int component, double t, double)
      {
        if(component == 0) // Neumann, (x,y) = (t,0), (nx,ny) = (0,-1)
          return 0.;
        if(component == 1) // Neumann, (x,y) = (1,t), (nx,ny) = (1,0)
          return (t - t * t * t) / 3.0;
        if(component == 2) // Neumann, (x,y) = (1-t,1), (nx,ny) = (0,1)
          return ((1 - t) * (1 - t) * (1 - t) - (1 - t)) / 3.0;
        if(component == 3) // Dirichlet, (x,y) = (0,1-t), (nx,ny) = (-1,0)
          return -(1 - t) * (1 - t) * (1 - t) / 3.0;
        return std::nan("1"); // only to avoid compiler warning
      };
      bd.push_back(BoundaryData(g));
      // the coefficient function
      f = [permeability](const Point&, double,
                         std::vector<double>& coeffs)
      {
        coeffs[0] = 1./permeability; // permeability
        coeffs[1] = 0.;              // right hand side, first component
        coeffs[2] = 0.;              // right hand side, second component
        coeffs[3] = 0.;              // divergence
      };
      break;
    }
    case 3: // obstacle with lower permeability than surrounding area
    {
      // 0 - circle in center
      // 1 - rectangle in center
      // 2 - diamond in center
      // 3 - two boxes, one at the top, one at the bottom
      // 4 - two half circles, one at the top, one at the bottom
      // 5 - two boxes, one at the top, one at the bottom, boxes have an offset
      // 6 - two half ellipsoids with an offset, one at the top, one at the
      //     bottom,
      // 7 - continuous change in permeability in form of a hump
      size_t obstacle = db["obstacle_type"];
      // inside the obstacle the porosity sigma is multiplied by this factor
      double factor = db["permeability_factor"];
      name = "obstacle\nThe permeability within the obstacle is "
             + std::to_string(factor) + " times the permeability outside."
             " The obstacle has the shape of";
      switch(obstacle)
      {
        case 0:
          name += " a circle in the center"; break;
        case 1:
          name += " a square in the center"; break;
        case 2:
          name += " a diamond in the center"; break;
        case 3:
          name += " two boxes on one vertical line, one at the top, "
                  "the other at the bottom\n"; break;
        case 4:
          name += " two half circles on one vertical line, one at the "
                  "top, one at the bottom"; break;
        case 5:
          name += " two boxes with a vertical offset, on at the top, one "
                  "at the bottom"; break;
        case 6:
          name += "two half ellipsoids with a vertical offset, on at the "
                  "top, one at the bottom"; break;
        case 7:
          name += "hump in the center (continuous change in permeability)"; 
          break;
        default:
          ErrThrow("Unknown obstacle. Set 'obstacle_type' in the input file to "
                   "0,...,7");
          break;
      }

      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));
      // exact solution, first velocity component (unknown)
      exact.push_back(AnalyticalFunction());
      // exact solution, second velocity component (unknown)
      exact.push_back(AnalyticalFunction());
      // exact solution, pressure component (unknown)
      exact.push_back(AnalyticalFunction());
      // boundary data
      auto g = [](unsigned int component, double, double)
      {
        if(component == 0) // Dirichlet (x,y) = (t,0), (nx,ny) = (0,-1)
          return 0.;
        if(component == 1) // Dirichlet, (x,y) = (1,t), (nx,ny) = (1,0)
          return 1.;
        if(component == 2) // Dirichlet, (x,y) = (1-t,1), (nx,ny) = (0,1)
          return 0.;
        if(component == 3) // Dirichlet, (x,y) = (0,1-t), (nx,ny) = (-1,0)
          return -1.;
        return std::nan("1"); // only to avoid compiler warning
      };
      bd.push_back(BoundaryData(g));
      // the coefficient function
      f = [permeability, obstacle, factor](const Point& point, double,
                                           std::vector<double>& coeffs)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        coeffs[0] = 1./permeability; // permeability
        switch(obstacle)
        {
          case 0:
            if((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.2 * 0.2)
              coeffs[0] /= factor;
            break;
          case 1:
            if(x > 0.45 && x < 0.55 && y > 0.35 && y < 0.65)
              coeffs[0] /= factor;
            break;
          case 2:
            if(x + y > 0.8 && x + y < 1.2 && x - y < 0.2 && x - y > -0.2)
              coeffs[0] /= factor;
            break;
          case 3:
            if(x > 0.35 && x < 0.65 && (y < 0.2 || y > 0.8))
              coeffs[0] /= factor;
            break;
          case 4:
            if((x - 0.5) * (x - 0.5) + y * y < 0.2 * 0.2
               || (x - 0.5) * (x - 0.5) + (y - 1) * (y - 1) < 0.2 * 0.2)
              coeffs[0] /= factor;
            break;
          case 5:
            if((x > 0.65 && x < 0.75 && y < 0.6)
               || (x > 0.25 && x < 0.35 && y > 0.4))
              coeffs[0] /= factor;
            break;
          case 6:
            if((x - 0.75) * (x - 0.75) + 0.1 * y * y < 0.2 * 0.2
               || (x - 0.25) * (x - 0.25) + 0.1 * (y - 1) * (y - 1) < 0.2 * 0.2)
              coeffs[0] /= factor;
            break;
          case 7:
          {
            double r = std::sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
            if(r < 0.2)
              coeffs[0] /= 1 + 0.5 * (factor - 1) * (std::cos(5 * M_PI * r) + 1);
            break;
          }
        }
        coeffs[1] = 0.; // right hand side, first component
        coeffs[2] = 0.; // right hand side, second component
        coeffs[3] = 0.; // divergence
      };
      break;
    }
    case 4: // 5-spot problem
    {
      // see Arif Masud, Thomas J.R. Hughes: "A stabilized mixed finite element
      // method for Darcy flow", chapter 4.2
      // the zones of the checkerboard domain have different permeabilities. The
      // quotient is this factor (1 means uniform permeability).
      double factor = 100;
      // boundary condition:
      bc.push_back(BoundaryCondition(DIRICHLET));
      // helping function
      auto r = [](double x, double y, double x0, double y0)
      {
        return (x - x0) * (x - x0) + (y - y0) * (y - y0);
      };
      // exact solution, first velocity component
      auto u1 = [r](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        double r0 = r(x, y, 0, 0);
        double r1 = r(x, y, 1, 1);
        double r2 = r(x, y, -1, 1);
        double r3 = r(x, y, -1, -1);
        double r4 = r(x, y, 1, -1);
        v.set((x / r0 - (x - 1) / r1 - (x + 1) / r2 - (x + 1) / r3
               - (x - 1) / r4) / 4,
              0, MultiIndex::D00);
        v.set(0., 0, MultiIndex::D10);
        v.set(-(x * y / (r0 * r0) - (x - 1) * (y - 1) / (r1 * r1)
                - (x + 1) * (y - 1) / (r2 * r2) - (x + 1) * (y + 1) / (r3 * r3)
                - (x - 1) * (y + 1) / (r4 * r4)) / 2,
              0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u1));
      // exact solution, second velocity component
      auto u2 = [r](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        double r0 = r(x, y, 0, 0);
        double r1 = r(x, y, 1, 1);
        double r2 = r(x, y, -1, 1);
        double r3 = r(x, y, -1, -1);
        double r4 = r(x, y, 1, -1);
        v.set((y / r0 - (y - 1) / r1 - (y - 1) / r2 - (y + 1) / r3
               - (y + 1) / r4) / 4,
              0, MultiIndex::D00);
        v.set(-(x * y / (r0 * r0) - (x - 1) * (y - 1) / (r1 * r1)
                - (x + 1) * (y - 1) / (r2 * r2) - (x + 1) * (y + 1) / (r3 * r3)
                - (x - 1) * (y + 1) / (r4 * r4)) / 2,
              0, MultiIndex::D10);
        v.set(0., 0, MultiIndex::D01);
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(u2));
      // exact solution, pressure component
      auto p = [r](const Point& point, FunctionEvaluation& v)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        size_t n_plus = 1;
        size_t n_minus = 4;
        double plus[]
            = {0,  0, 2,  0,  2,  2,  0,  2,  -2, 2,  -2, 0,  -2, -2, 0, -2, 2,
               -2, 4, 0,  4,  2,  4,  4,  2,  4,  0,  4,  -2, 4,  -4, 4, -4, 2,
               -4, 0, -4, -2, -4, -4, -2, -4, 0,  -4, 2,  -4, 4,  -4, 4, -2};
        double minus[]
            = {1,  1,  -1, 1,  -1, -1, 1,  -1, 3,  1,  3,  3,  1,  3,  -1,
               3,  -3, 3,  -3, 1,  -3, -1, -3, -3, -1, -3, 1,  -3, 3,  -3,
               3,  -1, 5,  1,  5,  3,  5,  5,  3,  5,  1,  5,  -1, 5,  -3,
               5,  -5, 5,  -5, 3,  -5, 1,  -5, -1, -5, -3, -5, -5, -3, -5,
               -1, -5, 1,  -5, 3,  -5, 5,  -5, 5,  -3, 5,  -1};
        // value
        double value = 0;
        for(size_t i = 0; i < n_plus; i++)
          value -= std::log(std::sqrt(r(x, y, plus[2 * i], plus[2 * i + 1])));
        for(size_t i = 0; i < n_minus; i++)
          value += std::log(std::sqrt(r(x, y, minus[2 * i], minus[2 * i + 1])));
        value *= 0.25;
        v.set(value, 0, MultiIndex::D00);

        // x-derivative
        value = 0;
        for(size_t i = 0; i < n_plus; i++)
          value -= (x - plus[2 * i]) / r(x, y, plus[2 * i], plus[2 * i + 1]);
        for(size_t i = 0; i < n_minus; i++)
          value += (x - minus[2 * i]) / r(x, y, minus[2 * i], minus[2 * i + 1]);
        value *= 0.25;
        v.set(value, 0, MultiIndex::D10);

        // y-derivative
        value = 0;
        for(size_t i = 0; i < n_plus; i++)
          value -= (y - plus[2 * i + 1])
                   / r(x, y, plus[2 * i], plus[2 * i + 1]);
        for(size_t i = 0; i < n_minus; i++)
          value += (y - minus[2 * i + 1])
                   / r(x, y, minus[2 * i], minus[2 * i + 1]);
        value *= 0.25;
        v.set(value, 0, MultiIndex::D01);

        // not computed here:
        v.set(0., 0, MultiIndex::D20);
        v.set(0., 0, MultiIndex::D11);
        v.set(0., 0, MultiIndex::D02);
      };
      exact.push_back(AnalyticalFunction(p));
      // boundary data, automatically computed from known solution
      // bd.push_back(BoundaryData(bc[0], exact));
      auto g = [](unsigned int component, double t, double)
      {
        double h = 1.0 / 32.0; // needs to be the mesh width!
        int order = 1;         // polynomial order (1 or 2 supported)
        switch(component)
        {
          case 0: // Dirichlet, (x,y) = (t,0), (nx,ny) = (0,-1)
            if(t < h)
            {
              if(order == 1)
                return -(1 - t / h) / (4 * h);
              if(order == 2)
                return -(3 / (8 * h)) * (1 - 2 * t / h + t * t / (h * h));
            }
            break;
          case 1:
            if(t > 1 - h)
            {
              if(order == 1)
                return (1 - (1 - t) / h) / (4 * h);
              if(order == 2)
                return (3 / (8 * h))
                       * (1 - 2 * (1 - t) / h + (1 - t) * (1 - t) / (h * h));
            }
            break;
          case 2:
            if(t < h)
            {
              if(order == 1)
                return (1 - t / h) / (4 * h);
              if(order == 2)
                return (3 / (8 * h)) * (1 - 2 * t / h + t * t / (h * h));
            }
            break;
          case 3:
            if(t > 1 - h)
            {
              if(order == 1)
                return -(1 - (1 - t) / h) / (4 * h);
              if(order == 2)
                return -(3 / (8 * h))
                       * (1 - 2 * (1 - t) / h + (1 - t) * (1 - t) / (h * h));
            }
            break;
          default:
            ErrThrow("wrong boundary component ", component);
            break;
        }
        return 0.;
      };
      bd.push_back(BoundaryData(g));
      // the coefficient function
      f = [permeability, factor](const Point& point, double,
                                 std::vector<double>& coeffs)
      {
        const double x = point.x;
        const double y = point.y;
        // isnan(point.z()) should be true, because it is a 2D example
        coeffs[0] = 1./permeability; // permeability
        if((x < 0.5 && y < 0.5) || (x > 0.5 && y > 0.5))
          coeffs[0] /= factor;
        coeffs[1] = 0.; // right hand side, first component
        coeffs[2] = 0.; // right hand side, second component
        coeffs[3] = 0.; // divergence
      };
      break;
    }
    default:
      ErrThrow("unknown example ", example);
      break;
  }

  PDECoefficients coeffs(f, Problem_type::Darcy, false, false);

  Output::info<2>("Example", "using the darcy example: ", name);
  Output::dash<2>("The permeability is ", permeability);
  return Example(std::move(bc), std::move(bd), std::move(coeffs),
                 std::move(exact));
}
