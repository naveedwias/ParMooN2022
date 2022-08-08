// Navier-Stokes problem, Harmonic polynomial
// 
// h(x,y) = h = x^3 - 3 x y^2
// p(x,y) = -0.5 (\nabla h \cdot \nabla h) +  14/5

void ExampleFile()
{
  Output::print<1>("Example: Harmonic_poly_degree.h");
  TDatabase::ParamDB->INTERNAL_QUAD_RULE=99;
}

double DIMENSIONLESS_VISCOSITY;

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double u1 = 3.*x*x*y-y*y*y; 
  double u1x = 6.*x*y; 
  double u1y=3*x*x-3.*y*y;
  values[0] = u1;
  values[1] = u1x;
  values[2] = u1y;
  values[3] = 0.;
}

void ExactU2(double x, double y, double *values)
{
  double u2 = x*x*x-3*y*y*x;  
  double u2x = 3.*x*x-3.*y*y; 
  double u2y=-6.*y*x;
  values[0] = u2;
  values[1] = u2x;
  values[2] = u2y;
  values[3] = 0.;
}

void ExactP(double x, double y, double *values)
{
  double r2 = x*x + y*y;
  double px = -3.*x*r2*r2; 
  double py = -3.*y*r2*r2;
  values[0] = -r2*r2*r2/2.0+1.0/8.0;
  values[1] = px;
  values[2] = py;
  values[3] = 0.;
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
  int sigma = TDatabase::ParamDB->P0;
  for(int i = 0; i < n_points; i++)
  {
    coeffs[i][0] = nu;    
    const double x=X[i]; const double y=Y[i];
    double u1  = 3.*x*x*y-y*y*y; 
    double u1x = 6.*x*y; 
    double u1y = 3*x*x-3.*y*y;
    
    double u2  = x*x*x-3.*y*y*x;  
    double u2x = 3.*x*x-3.*y*y; 
    double u2y = -6.*y*x;
    
    double px = -3.*x*(x*x+y*y)*(x*x+y*y); 
    double py = -3.*y*(x*x+y*y)*(x*x+y*y);
    
    if(TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR)
    {
      coeffs[i][1] = px + u1*u1x + u2*u1y + sigma * u1;
      coeffs[i][2] = py + u1*u2x + u2*u2y + sigma * u2;
      // cout<<px + u1*u1x + u2*u1y<< "  " << py + u1*u2x + u2*u2y<<endl;
      coeffs[i][4] = u1;
      coeffs[i][5] = u2;
    }
    coeffs[i][3] = u1x + u2y; // g (divergence)
  }
}

