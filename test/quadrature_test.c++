// TODO Add common header here!
#include "MooNMD_Io.h"
#include "QuadFormula.h"
#include <cmath>
#include <functional>
#include <list>
#include <utility> // std::pair, std::make_pair

// function type to be integrated
typedef std::function<double(parmoon::Point p)> function;
// a tuple of a function f and its integral values on two different domains
typedef std::tuple<function, double, double> function_and_integral;
// pair of quadrature weight and corresponding point
typedef std::pair<double, parmoon::Point> quad_pair;
// transformation of quadrature rule from reference cell to some other cell
typedef std::function<quad_pair(quad_pair)> transformation;

bool equal(double a, double b, double eps = 1.e-10)
{
  return std::abs(a - b) < eps;
}

// do a quadrature of a function f:R^3->R using a quadrature formula
double do_quadrature(const TQuadFormula& qf, function f)
{
  double integral = 0.;
  unsigned int n_q_points = qf.GetN_QuadPoints();
  for(unsigned int q = 0; q < n_q_points; ++q)
  {
    double w = qf.get_weight(q);
    parmoon::Point p = qf.get_point(q);
    integral += f(p) * w;
    Output::print<4>(" point (x,y,z) = (", p.x, ",", p.y, ",", p.z,
                     "), f(x,y,z) = ", f(p));
  }
  return integral;
}

TQuadFormula get_transformed_formula(TQuadFormula qf, transformation t)
{
  unsigned int n_points = qf.GetN_QuadPoints();
  for(unsigned int i = 0; i < n_points; ++i)
  {
    auto new_pair = t({qf.get_weight(i), qf.get_point(i)});
    qf.update_pair(i, new_pair);
  }
  return qf;
}

bool check(QuadratureFormula_type qf_type,
           const std::list<function_and_integral>& functions,
           transformation t = [](quad_pair wp){ return wp; })
{
  Output::print<4>();
  Output::print<3>("testing ", qf_type);
  TQuadFormula qf(qf_type);
  TQuadFormula qf_t{get_transformed_formula(qf, t)};
  
  bool ret = true, ret_t = true;
  for(const auto& it : functions)
  {
    double integral = do_quadrature(qf, std::get<0>(it));
    double integral_t = do_quadrature(qf_t, std::get<0>(it));
    Output::print<3>("Computed integral and transformed integral: ", integral,
                     ", ", integral_t);
    ret = equal(integral, std::get<1>(it));
    ret_t = equal(integral_t, std::get<2>(it));
    if(!ret)
    {
      Output::print("wrong result in quadrature ", integral, " while ",
                    std::get<1>(it), " was expected, diff: ",
                    integral - std::get<1>(it));
      qf.info();
      break;
    }
    if(!ret_t)
    {
      Output::print("wrong result in transformed quadrature ", integral_t,
                    " while ", std::get<2>(it), " was expected, diff: ",
                    integral_t - std::get<2>(it));
      qf_t.info();
      break;
    }
  }
  return ret && ret_t;
}

/////////////////////////////////////////////////////////////////////////////
//  1D                                                                     //
/////////////////////////////////////////////////////////////////////////////
bool check_all_1d()
{
  // all functions which shall be testes with various quadrature formulas
  std::list<function_and_integral> functions;
  transformation t = [](quad_pair wp)
                     { return quad_pair{0.5*wp.first, 0.5*(wp.second.x+1)}; };

  // extend the list of function by a polynomial of degree 1
  function f = [](parmoon::Point p) { return 3. * p.x + 1.; };
  functions.push_back({f, 2., 2.5});
  if(!check(QuadratureFormula_type::Gauss1Line, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 3
  f = [](parmoon::Point p)
      { return 3. * p.x * p.x * p.x - 2 * p.x * p.x + p.x + 1.; };
  functions.push_back({f, 2. / 3., 19. / 12.});
  if(!check(QuadratureFormula_type::Gauss2Line, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 5
  f = [](parmoon::Point p)
      { return 3. * p.x * p.x * p.x * p.x * p.x - 2 * p.x * p.x * p.x * p.x; };
  functions.push_back({f, -0.8, 0.1});
  if(!check(QuadratureFormula_type::Gauss3Line, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 7
  f = [](parmoon::Point p)
      { return 3. * p.x * p.x * p.x * p.x * p.x * p.x * p.x
            - 14. * p.x * p.x * p.x * p.x * p.x * p.x; };
  functions.push_back({f, -4., -13. / 8.});
  if(!check(QuadratureFormula_type::Gauss4Line, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 9
  f = [](parmoon::Point p)
      { return 3. * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x
            - 18. * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x; };
  functions.push_back({f, -4., -1.7});
  if(!check(QuadratureFormula_type::Gauss5Line, functions, t))
    return false;

  /// @todo more 1D quadrature tests
  // the rest is tested with the functions above, but these should also be
  // tested for higher order polynomials
  if(!check(QuadratureFormula_type::Gauss6Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss7Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss8Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss9Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss10Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss11Line, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss12Line, functions, t))
    return false;
  // Test integrate f:x^39-500x^38+3x^36
  // f = [](parmoon::Point p)
  //     {
  //       double tx = p.x * p.x * p.x * p.x * p.x * p.x;
  //       tx *= tx * tx * tx * tx * tx; // tx=x^36
  //       return tx * p.x * p.x * p.x - 500 * tx * p.x * p.x + 3 * tx;
  //     };
  // functions.push_back({f, -36766. / 1443.});
  // if(!check(QuadratureFormula_type::Gauss20Line, functions))
  //   return false;
  return true; // success
}


/////////////////////////////////////////////////////////////////////////////
//  2D                                                                     //
/////////////////////////////////////////////////////////////////////////////
bool check_all_2d()
{
  // all functions which shall be testes with various quadrature formulas
  std::list<function_and_integral> functions;
  transformation t = [](quad_pair wp)
                   {
                     parmoon::Point p(0.5*(wp.second.x+1), 0.5*(wp.second.y+1));
                     return quad_pair{0.25*wp.first, p};
                   };
  function f;

  /////////////////////////////////////////////////////////////////////////////
  // tests on quads:
  // extend the list of function by a polynomial of degree 3
  f = [](parmoon::Point p)
      {
        return 9 * (3. * p.x * p.x * p.x - 2 * p.x * p.x + p.x + 4)
              * (p.y * p.y * p.y + p.y * p.y - 2 * p.y + 1);
      };
  functions.push_back({f, 160., 385. / 16});
  if(!check(QuadratureFormula_type::Gauss2Quad, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 5
  f = [](parmoon::Point p)
      {
        return (3. * p.x * p.x * p.x * p.x * p.x - p.x * p.x * p.x * p.x + 1)
            * (p.y * p.y * p.y * p.y * p.y + p.y * p.y * p.y * p.y - p.y - 1);
      };
  functions.push_back({f, -2.56, -221. / 150.});
  if(!check(QuadratureFormula_type::Gauss3Quad, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 7
  f = [](parmoon::Point p)
      { 
        return (3 * p.x * p.x * p.x * p.x - 2 * p.x * p.x + 0.5)
             * (p.x * p.x * p.x + 2 * p.x) * (p.y * p.y * p.y * p.y + 1)
             * (p.y * p.y * p.y + 3); 
      };
  functions.push_back({f, 0., 53. / 20.});
  if(!check(QuadratureFormula_type::Gauss4Quad, functions, t))
    return false;

  /// @todo more 2D quadrature tests on quads
  // the rest is tested with the functions above, but these should also be
  // tested for higher order polynomials
  if(!check(QuadratureFormula_type::Gauss5Quad, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss6Quad, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss7Quad, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss8Quad, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss9Quad, functions, t))
    return false;
  functions.clear();
  /////////////////////////////////////////////////////////////////////////////
  // tests on triangles:
  // extend the list of function by a polynomial of degree 1
  f = [](parmoon::Point p) { return p.x + p.y + 1; };
  functions.push_back({f, 5. / 6., 7. / 24.});
  if(!check(QuadratureFormula_type::BaryCenterTria, functions, t))
    return false;
  if(!check(QuadratureFormula_type::VertexTria, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 2
  f = [](parmoon::Point p)
      { return 3. * p.x * p.x - 2 * p.x * p.y + 4 * p.y * p.y + 0.5; };
  functions.push_back({f, 0.75, 17. / 48.});
  if(!check(QuadratureFormula_type::MidPointTria, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 3
  f = [](parmoon::Point p)
      {
        return 3. * p.x * p.x * p.x - 2 * p.x * p.x * p.y + p.x * p.y * p.y
             + 4. * p.y * p.y * p.y + 0.5;
      };
  functions.push_back({f, 7. / 12., 119. / 384.});
  //if(!check(QuadratureFormula_type::SevenPointTria, functions, t))
  //  return false;
  if(!check(QuadratureFormula_type::CompGauss3Tria, functions, t))
    return false;
  if(!check(QuadratureFormula_type::CompGauss4Tria, functions, t))
    return false;
  // Gauss3Tria should be tested for polynomials of order 5
  if(!check(QuadratureFormula_type::Gauss3Tria, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 8
  f = [](parmoon::Point p)
      {
        return (3 * p.x + p.y) * (p.x * p.x - p.y * p.y)
             * (p.x * p.x * p.x + p.y * p.y * p.y) * (p.x * p.x - 2 * p.y * p.y);
      };
  functions.push_back({f, 253. / 4200., 4819305. / 147835568. });
  if(!check(QuadratureFormula_type::Degree8Tria, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss_Degree8Tria, functions, t))
    return false;

  /// @todo more 2D quadrature tests on triangles
  // the rest is tested with the functions above, but these should also be
  // tested for higher order polynomials
  if(!check(QuadratureFormula_type::Degree9Tria, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Degree11Tria, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Degree19Tria, functions, t))
    return false;

  return true; // success
}


/////////////////////////////////////////////////////////////////////////////
//  3D                                                                     //
/////////////////////////////////////////////////////////////////////////////
bool check_all_3d()
{
  // all functions which shall be testes with various quadrature formulas
  std::list<function_and_integral> functions;
  transformation t = [](quad_pair wp)
  {
    auto op = wp.second;
    parmoon::Point p(0.5*(op.x+1), 0.5*(op.y+1), 0.5*(op.z+1));
    return quad_pair{0.125*wp.first, p};
  };
  function f;

  /////////////////////////////////////////////////////////////////////////////
  // tests on hexas:
  // extend the list of function by a polynomial of degree 3
  f = [](parmoon::Point p)
      { return (3. * p.x - 4) * (p.y + 1) * (3 * p.z - 1.); };
  functions.push_back({f, 32., -15. / 8.});
  if(!check(QuadratureFormula_type::VertexHexa, functions, t))
    return false;
  if(!check(QuadratureFormula_type::VerticesAndOrigin, functions, t))
    return false;
  if(!check(QuadratureFormula_type::VerticesAndOrigin15, functions, t))
    return false;
  if(!check(QuadratureFormula_type::VerticesAndOrigin57, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Degree7_Points38, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 3
  f = [](parmoon::Point p)
      {
        return 1.125 * (3. * p.x * p.x * p.x - 2 * p.x * p.x + p.x + 4)
              * (p.y * p.y * p.y + p.y * p.y - 2 * p.y + 1)
              * (p.z * p.z * p.z + 3 * p.z * p.z - 1.);
      };
  functions.push_back({f, 0., 385. / 512.});
  if(!check(QuadratureFormula_type::Gauss2Hexa, functions, t))
    return false;

  // extend the list of function by a polynomial of degree 5
  f = [](parmoon::Point p)
      {
        return (3. * p.x * p.x * p.x * p.x * p.x 
                - 2 * p.x * p.x * p.x * p.x + p.x + 3)
             * (p.y * p.y * p.y * p.y * p.y + 2 * p.y * p.y * p.y * p.y
                - 2 * p.y + 1)
             * (p.z * p.z * p.z * p.z * p.z + 3 * p.z * p.z * p.z * p.z - 1.);
      };
  functions.push_back({f, -1456. / 125., -119. / 250.});
  if(!check(QuadratureFormula_type::Gauss3Hexa, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss4Hexa, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss5Hexa, functions, t))
    return false;
  // Test integrate f:(12x^11+x^10-4x-6)*(8y^11+4y^10-6)*(-z^11-10z^10+2)
  f = [](parmoon::Point p)
      {
        double t_x = p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x * p.x;
        double t_y = p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y * p.y;
        double t_z = p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z * p.z;
        return (12 * t_x * p.x + t_x - 4 * p.x - 6)
             * (8 * t_y * p.y + 4 * t_y - 6)
             * (-1 * t_z * p.z - 10 * t_z + 2);
      };
  functions.push_back({f, 386880. / 1331., 414428. / 11979.});
  if(!check(QuadratureFormula_type::Gauss6Hexa, functions, t))
    return false;
  /// @todo more 3D quadrature tests on hexahedra
  // the rest is tested with the functions above, but these should also be
  // tested for higher order polynomials
  if(!check(QuadratureFormula_type::Gauss7Hexa, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss8Hexa, functions, t))
    return false;
  if(!check(QuadratureFormula_type::Gauss9Hexa, functions, t))
    return false;
  
  functions.clear();
  /////////////////////////////////////////////////////////////////////////////
  // tests on tetrahedra:
  f = [](parmoon::Point p) { return p.x + p.y + p.z + 1.; };
  functions.push_back({f, 7./24, 23./384});
  if(!check(QuadratureFormula_type::BaryCenterTetra, functions, t))
    return false;
  if(!check(QuadratureFormula_type::VertexTetra, functions, t))
    return false;
  f = [](parmoon::Point p) { return p.x * p.y + p.y * p.z + p.z * p.z + 1.; };
  functions.push_back({f, 1./5., 29./640.});
  if(!check(QuadratureFormula_type::P2Tetra, functions, t))
    return false;
  f = [](parmoon::Point p) { return p.x * p.y * p.y * p.y + p.z * p.z + 1.; };
  functions.push_back({f, 31./168., 1747./53760.});
  if(!check(QuadratureFormula_type::P4Tetra, functions, t))
    return false;
  f = [](parmoon::Point p)
      { return p.x * p.y * p.y * p.y * p.z + p.z * p.z + 1.; };
  functions.push_back({f, 411./2240., 64289203./2061664679.});
  if(!check(QuadratureFormula_type::P5Tetra, functions, t))
    return false;
  f = [](parmoon::Point p)
      { return p.x * p.x * p.y * p.y * p.y * p.z * p.z * p.z + 1.; };
  functions.push_back({f, 92401./554400., 5252943./246459079.});
  if(!check(QuadratureFormula_type::P8Tetra, functions, t))
    return false;
  
  return true; // success
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
#ifdef _MPI
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
#else
int main(int , char* [])
{
#endif // _MPI

  Output::setVerbosity(2);

  if(!check_all_1d())
  {
    return 1; // failure
  }

  if(!check_all_2d())
  {
    return 1; // failure
  }

  if(!check_all_3d())
  {
    return 1; // failure
  }

  Output::print("successful test");
#ifdef _MPI
  MPI_Finalize();
#endif // _MPI
  return 0;
}
