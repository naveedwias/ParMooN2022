// ======================================================================
// Pure convection-reaction, no diffusion
// Based on Houston, Schwab, Sueli, SIAM 39(6), 2002, doi: 10.1137/S0036142900374111
// ======================================================================

#include "Constants.h"
#include "MooNMD_Io.h"
#include <cmath>
void ExampleFile()
{
  Output::root_info<1>("Example", "Pure_Convection.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  const double p = M_PI; // Pi;

  double x_hat = 2*x-1; // transform from (0,1) to (-1,1)
  double y_hat = 2*y-1; // transform from (0,1) to (-1,1)
  double z_hat = 2*z-1; // transform from (0,1) to (-1,1)

  double arg = p*(1+x_hat)*(1+y_hat)*(1+z_hat)/8; // argument of sin
  values[0] = 1+std::sin(arg);  //original solution

  double arg_x = p*(1+y_hat)*(1+z_hat)/4; // inner derivative wrt x
  double arg_y = p*(1+x_hat)*(1+z_hat)/4; // inner derivative wrt y
  double arg_z = p*(1+x_hat)*(1+y_hat)/4; // inner derivative wrt z

  values[1] = std::cos(arg)*arg_x; //x-derivative
  values[2] = std::cos(arg)*arg_y; //y-derivative
  values[3] = std::cos(arg)*arg_z; //z-derivative
  double dxx = -std::sin(arg)*arg_x*arg_x;
  double dyy = -std::sin(arg)*arg_y*arg_y;
  double dzz = -std::sin(arg)*arg_z*arg_z;
  values[4] = dxx + dyy + dzz ; //laplacian
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  double values[5];
  Exact(x, y, z, values);
  value = values[0];
}

inline void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *z, const double *const*, double **coeffs)
{
  double values[5];
  for(int i = 0; i < n_points; i++)
  {
    auto coeff = coeffs[i];

    auto x_hat = 2*x[i]-1; // transform from (0,1) to (-1,1)
    auto y_hat = 2*y[i]-1; // transform from (0,1) to (-1,1)
    auto z_hat = 2*z[i]-1; // transform from (0,1) to (-1,1)

    coeff[0] = 0;
    coeff[1] = 2-y_hat*z_hat; //x-coordiante of convection b
    coeff[2] = 2-x_hat*z_hat; //y-coordiante of convection b
    coeff[3] = 2-x_hat*z_hat; //z-coordiante of convection b
    coeff[4] = 1+(1+x_hat)*(1+y_hat)*(1+z_hat); //reaction coeffecient c

    Exact(x[i], y[i], z[i], values);

    coeff[5] = -coeff[0]*values[4]; // diffusion
    coeff[5] += coeff[1]*values[1] + coeff[2]*values[2] + coeff[3]*values[3]; // convection
    coeff[5] += coeff[4]*values[0]; // reaction
  }
}
