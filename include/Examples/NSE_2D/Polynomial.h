// Navier-Stokes problem with the solution from the BDM space of a given degree.
// The exact velocity is just a polynomial of a given degree and the exact
// pressure is a polynomial of degree -1.
//

#include "Constants.h"
#include "MooNMD_Io.h"
#include <cmath>
#include <vector>
unsigned int deg;
bool is_RT;
double pressure_factor;

// Compute for the given degree the number of coefficients of polynomial and set
// the coefficients as 1, -2, 3, -4, ..., -8, 9, -1, 2, -3, ...
std::vector<int> coef;
unsigned int n_coef ;

void ExampleFile()
{
  if (is_RT)
  {
    Output::print<1>("Example: RT",deg, " polynomial");
  }
  else
  {
    Output::print<1>("Example: 2D polynomial of degree ", deg);
  }
  // number of coeffs of a 2D polynomial of order deg
  auto n_coef = (deg*deg + 3*deg)/2 + 1;
  n_coef += 3;  // add more coefficients to differ between components and pressure
  coef.resize(n_coef, 0);
  for (unsigned int i = 0; i < n_coef; ++i)
  {
    coef[i]= i%9 + 1;
    if (i%2 == 1)
    {
      coef[i] *= -1;  // this may lead to a slower blow up of the solution
    } // end if
  } // end for

}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  // The polynomial is given by a sum of coefficient times powers of x and y.
  // The coefficients come from the previously defined coef vector.

  // values[0] = function at (x,y)
  unsigned int index = 0;
  values[0] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      values[0] += coef[index] * std::pow(x,i) * std::pow(y,j);
      values[0] += (is_RT) ?
        coef[index + 2] * std::pow(x,i) * std::pow(y,j) * x : 0;
      index++;
    }
  }

  // values[1] = derivative wrt. x of function at (x,y)
  index = 0;
  values[1] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      values[1] += (is_RT) ?
        coef[index + 2] * std::pow(x,i) * std::pow(y,j) : 0;
      if (i == 0)
      {
        index++;
      }
      else
      {
        values[1] += coef[index] * (int) i * std::pow(x,i-1) * std::pow(y,j);
        values[1] += (is_RT) ?
          coef[index + 2] * (int) i * std::pow(x,i-1) * std::pow(y,j) * x : 0;
        index++;
      }
    }
  }


  // values[2] = derivative wrt. y of function at (x,y)
  index = 0;
  values[2] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (j == 0)
      {
        index++;
      }
      else
      {
        values[2] += coef[index] * (int) j * std::pow(x,i) * std::pow(y,j-1);
        values[2] += (is_RT) ?
          coef[index + 2] * (int) j * std::pow(x,i) * std::pow(y,j-1) * x : 0;
        index++;
      }
    }
  }

  // values[3] = Laplacian of function at (x,y)
  index = 0;
  values[3] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. x
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (i == 0)
      {
        index++;
      }
      else
      {
        if (i > 1)
        {
          values[3] += coef[index] * (int) i * (int) (i-1) * std::pow(x,i-2) * std::pow(y,j);
          values[3] += (is_RT) ?
            coef[index + 2] * (int) i * (int) (i-1) * std::pow(x,i-2) * std::pow(y,j) * x : 0;
        }
        values[3] += (is_RT) ?
          2 * coef[index + 2] * (int) i * std::pow(x,i-1) * std::pow(y,j): 0;
        index++;
      }
    }
  }
  index = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. y
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (j <= 1)
      {
        index++;
      }
      else
      {
        values[3] += coef[index] * (int) j * (int) (j-1) * std::pow(x,i) * std::pow(y,j-2);
        values[3] += (is_RT) ?
          coef[index + 2] * (int) j * (int) (j-1) * std::pow(x,i) * std::pow(y,j-2) : 0;
        index++;
      }
    }
  }
}

void ExactU2(double x, double y, double *values)
{
  // The polynomial is given by a sum of coefficient times powers of x and y.
  // The coefficients come from the previously defined coef vector times -1.

  // values[0] = function at (x,y)
  unsigned int index = 0;
  values[0] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      values[0] += coef[index + 1] * std::pow(x,i) * std::pow(y,j);
      values[0] += (is_RT) ?
        coef[index + 2] * std::pow(x,i) * std::pow(y,j) * y : 0;
      index++;
    }
  }

  // values[1] = derivative wrt. x of function at (x,y)
  index = 0;
  values[1] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (i == 0)
      {
        index++;
      }
      else
      {
        values[1] += coef[index + 1] * (int) i * std::pow(x,i-1) * std::pow(y,j);
        values[1] += (is_RT) ?
          coef[index + 2] * (int) i * std::pow(x,i-1) * std::pow(y,j) * y : 0;
        index++;
      }
    }
  }


  // values[2] = derivative wrt. y of function at (x,y)
  index = 0;
  values[2] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  {
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      values[2] += (is_RT) ?
        coef[index + 2] * std::pow(x,i) * std::pow(y,j) : 0;
      if (j == 0)
      {
        index++;
      }
      else
      {
        values[2] += coef[index + 1] * (int) j * std::pow(x,i) * std::pow(y,j-1);
        values[2] += (is_RT) ?
          coef[index + 2] * (int) j * std::pow(x,i) * std::pow(y,j-1) * y : 0;
        index++;
      }
    }
  }

  // values[3] = Laplacian of function at (x,y)
  index = 0;
  values[3] = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. x
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (i <= 1)
      {
        index++;
      }
      else
      {
        values[3] += coef[index + 1] * (int) i * (int) (i-1) * std::pow(x,i-2) * std::pow(y,j);
        values[3] += (is_RT) ?
          coef[index + 2] * (int) i * (int) (i-1) * std::pow(x,i-2) * std::pow(y,j) * y : 0;
        index++;
      }
    }
  }
  index = 0;
  for (unsigned int i = 0; i <= deg; ++i)
  { // 2nd derivative wrt. y
    for (unsigned int j = 0; j <= deg-i; ++j)
    {
      if (j == 0)
      {
        index++;
      }
      else
      {
        if (j > 1)
        {
        values[3] += coef[index + 1] * (int) j * (int) (j-1) * std::pow(x,i) * std::pow(y,j-2); values[3] += (is_RT) ?
          coef[index + 2] * (int) j * (int) (j-1) * std::pow(x,i) * std::pow(y,j-2) * y : 0;
        }
        values[3] += (is_RT) ?
          2 * coef[index + 2] * (int) j * std::pow(x,i) * std::pow(y,j-1) : 0;
        index++;
      }
    }
  }
}

void ExactP(double x, double y, double *values)
{
  // The polynomial is given by a sum of coefficient times powers of x and
  // y such that the degree is deg - 1.
  // The coefficients come from the previously used coef vector times a factor
  // of +-1.

  auto degree = (is_RT) ? deg + 1 : deg;

  // values[0] = function at (x,y)
  unsigned int index = 0;
  values[0] = 0;
  for (unsigned int i = 0; i < degree; ++i)
  {
    for (unsigned int j = 0; j < degree-i; ++j)
    {
      auto factor = ((i + j) % 2 == 0) ? 1 : -1;
      double integral = 1. / (i+1) / (j+1);
      values[0] += pressure_factor * factor * coef[index + 2] * (std::pow(x,i) * std::pow(y,j) - integral);
      index++;
    }
  }

  // values[1] = derivative wrt. x of function at (x,y)
  index = 0;
  values[1] = 0;
  for (unsigned int i = 0; i < degree; ++i)
  {
    for (unsigned int j = 0; j < degree-i; ++j)
    {
      if (i == 0)
      {
        index++;
      }
      else
      {
        auto factor = ((i + j) % 2 == 0) ? 1 : -1;
        values[1] += pressure_factor * factor * coef[index + 2] * (int) i * std::pow(x,i-1)
          * std::pow(y,j);
        index++;
      }
    }
  }


  // values[2] = derivative wrt. y of function at (x,y)
  index = 0;
  values[2] = 0;
  for (unsigned int i = 0; i < degree; ++i)
  {
    for (unsigned int j = 0; j < degree-i; ++j)
    {
      if (j == 0)
      {
        index++;
      }
      else
      {
        auto factor = ((i + j) % 2 == 0) ? 1 : -1;
        values[2] += pressure_factor * factor * coef[index + 2] * (int) j * std::pow(x,i)
          * std::pow(y,j-1);
        index++;
      }
    }
  }

  // values[3] = Laplacian of function at (x,y)
  index = 0;
  values[3] = 0;
  for (unsigned int i = 0; i < degree; ++i)
  { // 2nd derivative wrt. x
    for (unsigned int j = 0; j < degree-i; ++j)
    {
      if (i <= 1)
      {
        index++;
      }
      else
      {
        auto factor = ((i + j) % 2 == 0) ? 1 : -1;
        values[3] += pressure_factor * factor * coef[index + 2] * (int) i * (int) (i-1)
          * std::pow(x,i-2) * std::pow(y,j);
        index++;
      }
    }
  }
  index = 0;
  for (unsigned int i = 0; i < degree; ++i)
  { // 2nd derivative wrt. y
    for (unsigned int j = 0; j < degree-i; ++j)
    {
      if (j <= 1)
      {
        index++;
      }
      else
      {
        auto factor = ((i + j) % 2 == 0) ? 1 : -1;
        values[3] += pressure_factor * factor * coef[index + 2] * (int) j * (int) (j-1)
          * std::pow(x,i) * std::pow(y,j-2);
        index++;
      }
    }
  }
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void transform(const int BdComp, const double Param, double& x, double& y,
               double& nx, double& ny)
{
  switch(BdComp)
  {
    case 0:
      x = Param;
      y = 0.;
      nx = 0.;
      ny = -1.;
      break;
    case 1:
      x = 1.;
      y = Param;
      nx = 1.;
      ny = 0.;
      break;
    case 2:
      x = 1. - Param;
      y = 1.;
      nx = 0.;
      ny = 1.;
      break;
    case 3:
      x = 0.;
      y = 1. - Param;
      nx = -1.;
      ny = 0.;
      break;
    default:
      ErrThrow("wrong boundary part number", BdComp);
      break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  if(cond == DIRICHLET)
    value = u1[0];
  else
  {
    // NEUMANN
    value = nu * (nx * u1[1] + ny * u1[2]) - p[0] * nx;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * nu * (nx * u1[1] + ny * u2[1]);
    }
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  if(cond == DIRICHLET)
    value = u2[0];
  else
  {
    // NEUMANN
    value = nu * (nx * u2[1] + ny * u2[2]) - p[0] * ny;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * nu * (nx * u1[2] + ny * u2[2]);
    }
  }
  return;
}

// u times normal vector
void UNormalValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  
  if(cond == DIRICHLET)
    value = nx * u1[0] + ny * u2[0];
  else
  {
    // NEUMANN
    value = nu * (nx * (nx*u1[1] + ny*u1[2]) + ny * (nx*u2[1] + ny*u2[2])) - p[0];
    if(lt == 1)
    {
      ErrThrow("LAPLACETYPE=1 not yet supported");
      //value *= 0.5;
      //value += 0.5 * nu * (nx * u1[1] + ny * u2[1]);
    }
  }
  return;
}

// u times tangential vector
void UTangValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
  const double nu = 1./TDatabase::ParamDB->RE_NR;
  // find out boundary condition at the evaluation point on the boundary
  BoundCond cond;
  BoundCondition(BdComp, Param, cond);
  // find coordinates and normal of evaluation point on the boundary
  double x, y, nx, ny;
  transform(BdComp, Param, x, y, nx, ny);
  // evaluate the exact solution at the given point
  double u1[4];
  double u2[4];
  double p[4];
  ExactU1(x, y, u1);
  ExactU2(x, y, u2);
  ExactP(x, y, p);
  double tx = -ny;
  double ty = nx;
  
  if(cond == DIRICHLET)
    value = tx * u1[0] + ty * u2[0];
  else
  {
    // NEUMANN
    value = nu * (tx * (nx*u1[1] + ny*u1[2]) + ty * (nx*u2[1] + ny*u2[2]));
    if(lt == 1)
    {
      ErrThrow("LAPLACETYPE=1 not yet supported");
      //value *= 0.5;
      //value += 0.5 * nu * (nx * u1[2] + ny * u2[2]);
    }
  }
  return;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  const double nu=1./TDatabase::ParamDB->RE_NR;
  double val1[4];
  double val2[4];
  double val3[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;  // viscosity

    ExactU1(X[i], Y[i], val1);
    ExactU2(X[i], Y[i], val2);
    ExactP(X[i], Y[i], val3);

    coeffs[i][1] = -nu*val1[3] + val3[1]; // f1
    coeffs[i][2] = -nu*val2[3] + val3[2]; // f2

    if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5) // Navier-Stokes (3 means Stokes)
    {
      coeffs[i][1] += val1[0]*val1[1] + val2[0]*val1[2]; // f1
      coeffs[i][2] += val1[0]*val2[1] + val2[0]*val2[2]; // f2
    }
    coeffs[i][3] = val1[1] + val2[2]; // g (divergence)

    // additional coefficient (used only in the Brinkman problem)
    coeffs[i][4] = 0.;
  }
}
