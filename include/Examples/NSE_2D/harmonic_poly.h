// Navier-Stokes problem, Harmonic polynomial
// 
// h(x,y) = h = x^3 - 3 x y^2
// p(x,y) = -0.5 (\nabla h \cdot \nabla h) +  14/5

void ExampleFile()
{
  Output::print<1>("Example: Harmonic_poly.h");
}

double DIMENSIONLESS_VISCOSITY;

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double u = 3.*x*x-3.*y*y;
  double ux = 6.*x;
  double uy = -6.*y;
  values[0] = u;
  values[1] = ux;
  values[2] = uy;
  values[3] = 0.0;
}

void ExactU2(double x, double y, double *values)
{
  double u = -6.*x*y;
  double ux = -6.*y;
  double uy = -6.*x;
  values[0] = u;
  values[1] = ux;
  values[2] = uy;
  values[3] = 0.0;
}

void ExactP(double x, double y, double *values)
{
  double p  = -0.5* ((3.*x*x-3.*y*y)*( 3.*x*x-3.*y*y) + 36.*x*x*y*y)+ 14./5.;
  double px = -0.5* (2.*(3.*x*x-3.*y*y)*( 6.*x) + 72.*x*y*y);
  double py = -0.5* (2.*(3.*x*x-3.*y*y)*(-6.*y) + 72.*x*x*y);
  values[0] = p;
  values[1] = px;//-0.5* (2. *  (3.*x*x-3.*y*y)*( 6.*x) + 72.*x*y*y);
  values[2] = py;//-0.5* (2.*(3.*x*x-3.*y*y)*(-6.*y) + 72.*x*x*y);
  values[3] = 0;
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
  const double nu = DIMENSIONLESS_VISCOSITY;
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
  const double nu = DIMENSIONLESS_VISCOSITY;
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

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  static double nu = DIMENSIONLESS_VISCOSITY;
  double sigma = TDatabase::ParamDB->P0;
  double u1[4], u2[4], p[4];
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);

    if(TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR)
    {
      coeffs[i][1] = -nu*u1[3] + p[1] + u1[0]*u1[1]+u2[0]*u1[2] + sigma*u1[0];
      coeffs[i][2] = -nu*u2[3] + p[2] + u1[0]*u2[1]+u2[0]*u2[2] + sigma*u2[0];      
      coeffs[i][4] = u1[0];
      coeffs[i][5] = u2[0];
      coeffs[i][6] = sigma;
    }
    else
    {
      coeffs[i][1] = -nu*u1[3] + p[1] + u1[0]*u1[1]+u2[0]*u1[2];
      coeffs[i][2] = -nu*u2[3] + p[2] + u1[0]*u2[1]+u2[0]*u2[2];
    }
    coeffs[i][3] = u1[1] + u2[2]; // g (divergence)
  }//exit(0);
}

