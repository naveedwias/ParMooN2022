// ======================================================================
// Polynomial problem for the UnitCube. This MAYBE works on other hexahedrons
// but it is not tested yet.
//
// Idea: We want to have a polynomial of degree at most k mapping from R^3 to R.
// The coefficients are alternately multiplied by -1 to prevent a blow up within
// (-1, 1)^3. So the exact solution is
// coef[l]*x^i * y^j * z^k
// for some coefficient function coef.
// ======================================================================
#include "Constants.h"
#include "MooNMD_Io.h"
#include <cmath>

double PECLET_NUMBER;

unsigned int deg; // degree of polynomial
std::vector<int> coef;
unsigned int n_coef ;

void ExampleFile()
{
  Output::print<1>("Example: Polynomial.h with polynomial degree ", deg);

  // Determine dimension of basis of polynomials mapping from R^3->R of degree
  // at most k, i.e. the number of coefficients. The dimension in general is
  // given by
  // (k+d)! / (k! * (k+d-k)!)
  // where ( )! is the factorial, d is the dimension (here 3) and k the degree.
  // For d = 3 the formula can be manipulated to arrive at
  // (k+1)*(k+2)*(k+3)/6.
  // The polynomial is later given by coef[l]*x^i*y^j*z^k, where i+j+k <= deg
  // and l is a function of i, j and k.
  // The coefficients are stored in coef and are given as described in the
  // beginning of the file.
  auto n_coef = (deg+1) * (deg+2) * (deg+3) / 6;
  coef.resize(n_coef, 0);
  for (unsigned int i = 0; i < n_coef; ++i) {
    coef[i]= i%9 + 1;
    if (i%2 == 1)
    {
      coef[i] *= -1;  // this may lead to a slower blow up of the solution
    } // end if
  } // end for
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  // values[0] = function at (x,y,z)
  unsigned int index = 0;
  values[0] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        values[0] += coef[index] * std::pow(x,i) * std::pow(y,j)
          * std::pow(z,k);
        index++;
      }
    }
  }

  index = 0;
  values[1] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (i == 0)
        {
          index++;
        }
        else
        {
          values[1] += coef[index] * (int) i * std::pow(x,i-1) * std::pow(y,j)
            * std::pow(z,k);
          index++;
        }
      }
    }
  }

  index = 0;
  values[2] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (j == 0)
        {
          index++;
        }
        else
        {
          values[2] += coef[index] * (int) j * std::pow(x,i) * std::pow(y,j-1)
            * std::pow(z,k);
          index++;
        }
      }
    }
  }

  index = 0;
  values[3] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (k == 0)
        {
          index++;
        }
        else
        {
          values[3] += coef[index] * (int) k * std::pow(x,i) * std::pow(y,j)
            * std::pow(z,k-1);
          index++;
        }
      }
    }
  }

  // values[4] = Laplacian of function at (x,y,z)
  index = 0;
  values[4] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. x
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (i <= 1)
        {
          index++;
        }
        else
        {
          values[4] += coef[index] * (int) i * (int) (i-1) * std::pow(x,i-2)
            * std::pow(y,j) * std::pow(z,k);
          index++;
        }
      }
    }
  }

  index = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. y
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (j <= 1)
        {
          index++;
        }
        else
        {
          values[4] += coef[index] * (int) j * (int) (j-1) * std::pow(x,i)
            * std::pow(y,j-2) * std::pow(z,k);
          index++;
        }
      }
    }
  }

  index = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. z
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      for (unsigned int k = 0; k <= deg-i-j; ++k)
      {
        if (k <= 1)
        {
          index++;
        }
        else
        {
          values[4] += coef[index] * (int) k * (int) (k-1) * std::pow(x,i)
            * std::pow(y,j) * std::pow(z,k-2);
          index++;
        }
      }
    }
  }
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

void BilinearCoeffs(int n_points, const double *x, const double *y,
                    const double *z, const double *const*, double **coeffs)
{
  const double eps = PECLET_NUMBER;
  double values[5];
  for(int i = 0; i < n_points; i++)
  {
    Exact(x[i], y[i], z[i], values);
    coeffs[i][0] = eps;
    coeffs[i][1] = 1;
    coeffs[i][2] = 2;
    coeffs[i][3] = 3;
    coeffs[i][4] = 4;
    coeffs[i][5] = -coeffs[i][0] * values[4] + coeffs[i][1]*values[1]
                   + coeffs[i][2]*values[2] + coeffs[i][3]*values[3]
                   + coeffs[i][4]*values[0];
  }
}
