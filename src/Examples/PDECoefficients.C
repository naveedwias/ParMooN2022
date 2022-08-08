#include <cmath> // std::isnan
#include <PDECoefficients.h>
#include <MooNMD_Io.h>

using namespace parmoon;

/* ************************************************************************* */
std::ostream& operator<<(std::ostream& out, const Problem_type value)
{
  const char* s = 0;
#define PROCESS_VAL(p)                                                         \
  case(Problem_type::p):                                                       \
    s = #p;                                                                    \
    break;
  switch(value)
  {
    PROCESS_VAL(ConvDiffReac);
    PROCESS_VAL(Stokes);
    PROCESS_VAL(NavierStokes);
    PROCESS_VAL(Darcy);
    PROCESS_VAL(dummy);
    default:
      s = "unknown Problem_type ";
      break;
  }
#undef PROCESS_VAL
  return out << s;
}

/* ************************************************************************* */
PDECoefficients::PDECoefficients(
    std::function<void(const Point&, double, std::vector<double>&)> f,
    Problem_type t, bool, bool)
  : type(t), coeff(f)
{
}

/* ************************************************************************* */
unsigned int n_coeffs(Problem_type t, unsigned int dim)
{
  if(dim == 2)
  {
    switch(t)
    {
      case Problem_type::ConvDiffReac:
        return 5;
      case Problem_type::Stokes:
      case Problem_type::NavierStokes:
        return 4;
      case Problem_type::Darcy:
        return 7;
      default:
      {
        ErrThrow("unsupported Problem_type ", t);
      }
    }
  }
  else if(dim == 3)
  {
    switch(t)
    {
      case Problem_type::ConvDiffReac:
        return 6;
      case Problem_type::Stokes:
      case Problem_type::NavierStokes:
        return 5;
      case Problem_type::Darcy:
        return 13;
      default:
      {
        ErrThrow("unsupported Problem_type ", t);
      }
    }
  }
  else
  {
    ErrThrow("unknown dimension ", dim);
  }
}

/* ************************************************************************* */
void PDECoefficients::get(const Point& point, double time,
                          std::vector<double>& coeffs) const
{
  // check the size of coeffs
  // if the z-component of the point is nan, it is 2D otherwise 3D
  unsigned int n_c = n_coeffs(this->type, (std::isnan(point.z) ? 2 : 3));
  if(coeffs.size() < n_c)
  {
    coeffs.resize(n_c);
  }
  if(this->coeff) // callable
  {
    this->coeff(point, time, coeffs);
  }
  else
  {
    ErrThrow("no coefficients defined");
  }
}
