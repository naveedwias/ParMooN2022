// Navier-Stokes problem with sine and cosine functions
// 

#include "Constants.h"
#include "MooNMD_Io.h"

double pressure_factor;


void ExampleFile()
{
  bool ns = (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == 5);
  Output::print<1>("Example: 2d ", ns ? "Navier-" : "", "Stokes, ",
                   "polynomial_solution.h.");
  Output::print<1>("         This is Example D.3 in Volker John's book "
                   "\"Finite Element Methods for Incompressible Flow Problems\""
                   ", 2016");
}

double DIMENSIONLESS_VISCOSITY = 1.;

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double a = x - 1;
  double b = y - 1;
  values[0] = 1000 * x*x * y*y * a*a*a*a * b * (5*y - 3);
  values[1] = 2000 * x * y*y * a*a*a * (3*x-1) * b * (5*y-3);
  values[2] = 1000 * x*x * y * a*a*a*a * (5*y*b + y*(5*y-3) + 2*b*(5*y-3));
  values[3] = 2000*a*a*(30*x*x*y*y*a*a - 24*x*x*y*a*a + 3*x*x*a*a
                          + y*y*b*(5*y - 3)*(6*x*x + 8*x*a + a*a));
}

void ExactU2(double x, double y, double *values)
{
  double a = x - 1;
  double b = y - 1;
  values[0] = -2000 * x * y*y*y * a*a*a * (3*x-1) * b*b;
  values[1] = 2000 * y*y*y * a*a * b*b * (-3*x*a - 3*x*(3*x-1) + (-3*x+1)*a);
  values[2] = 2000 * x * y*y * a*a*a * (3*x-1) * (-5*y + 3) * b;
  values[3] = -4000 * y * a * (x*a*a*(3*x - 1)*(y*y + 6*y*b + 3*b*b)
                               + 3*y*y*b*b*(3*x*a + x*(3*x-1) + a*a + a
                                 *(3*x-1)));
}

void ExactP(double x, double y, double *values)
{
  values[0] = pressure_factor * M_PI*M_PI*x*y*(-x*std::sin(2*M_PI*x*y) + y*y*std::cos(2*M_PI*x*x*y)) + 0.125;
  values[1] = pressure_factor * M_PI*M_PI*y*(-4*M_PI*x*x*y*y*y*std::sin(2*M_PI*x*x*y)
                       -2*M_PI*x*x*y*std::cos(2*M_PI*x*y) - 2*x*std::sin(2*M_PI*x*y)
                       +y*y*std::cos(2*M_PI*x*x*y));
  values[2] = pressure_factor * M_PI*M_PI*x*(-2*M_PI*x*x*y*y*y*std::sin(2*M_PI*x*x*y)
                       -2*M_PI*x*x*y*std::cos(2*M_PI*x*y) - x*std::sin(2*M_PI*x*y)
                       +3*y*y*std::cos(2*M_PI*x*x*y));
  values[3] = pressure_factor * 2*M_PI*M_PI*(-2*M_PI*M_PI*x*x*x*x*x*y*y*y*std::cos(2*M_PI*x*x*y)
                       +2*M_PI*M_PI*x*x*x*x*y*std::sin(2*M_PI*x*y)
                       -8*M_PI*M_PI*x*x*x*y*y*y*y*y*std::cos(2*M_PI*x*x*y)
                       -6*M_PI*x*x*x*y*y*std::sin(2*M_PI*x*x*y)
                       -2*M_PI*x*x*x*std::cos(2*M_PI*x*y)
                       +2*M_PI*M_PI*x*x*y*y*y*std::sin(2*M_PI*x*y)
                       -6*M_PI*x*y*y*y*y*std::sin(2*M_PI*x*x*y)
                       -4*M_PI*x*y*y*std::cos(2*M_PI*x*y) + 3*x*y*std::cos(2*M_PI*x*x*y)
                       -y*std::sin(2*M_PI*x*y));
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
    value = DIMENSIONLESS_VISCOSITY * (nx * u1[1] + ny * u1[2]) - p[0] * nx;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * DIMENSIONLESS_VISCOSITY * (nx * u1[1] + ny * u2[1]);
    }
  }
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  // Neumann boundary means setting T.n where n is outer normal and T is either
  // 2nu D(u) - p  (for LAPLACETYPE==1) or nu grad(u) - p   (for LAPLACETYPE==0)
  const int lt = TDatabase::ParamDB->LAPLACETYPE;
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
    value = DIMENSIONLESS_VISCOSITY * (nx * u2[1] + ny * u2[2]) - p[0] * ny;
    if(lt == 1)
    {
      value *= 0.5;
      value += 0.5 * DIMENSIONLESS_VISCOSITY * (nx * u1[2] + ny * u2[2]);
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
  double val1[4];
  double val2[4];
  double val3[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    
    ExactU1(X[i], Y[i], val1);
    ExactU2(X[i], Y[i], val2);
    ExactP(X[i], Y[i], val3);
    
    coeffs[i][1] = -DIMENSIONLESS_VISCOSITY*val1[3] + val3[1]; // f1
    coeffs[i][2] = -DIMENSIONLESS_VISCOSITY*val2[3] + val3[2]; // f2
    
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
