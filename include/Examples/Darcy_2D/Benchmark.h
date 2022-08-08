// Bechmark for Darcy problem, exact solution is
// 
// u(x,y) = (y,x)
// p(x,y) = 0

void ExampleFile()
{
  Output::print<1>("Example: Benchmark.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double y, double *values)
{
  values[0] = y;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
}

void ExactU2(double x, double, double *values)
{
  values[0] = x;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int bdComp, double, BoundCond &cond)  
{
  cond = (bdComp == 0) ? NEUMANN : DIRICHLET;
}

// u \cdot n
void UNBoundValue(int BdComp, double t, double &value)
{
  double x, y, nx, ny;
  switch(BdComp)
  {
    case 0:
      x = t; y = 0.;
      nx = 0.; ny = -1.;
      break;
    case 1: 
      x = 1.; y = t;
      nx = 1.; ny = 0.;
      break;
    case 2:
      x = 1. - t; y = 1.;
      nx = 0.; ny = 1.;
      break;
    case 3:
      x = 0; y = 1. - t;
      nx = -1.; ny = 0.;
      break;
    default: 
      ErrThrow("wrong boundary part number ", BdComp);
      break;
  }
  BoundCond cond;
  BoundCondition(BdComp, t, cond);
  if(cond == DIRICHLET)
  {
    double u1[4], u2[4];
    ExactU1(x, y, u1);
    ExactU2(x, y, u2);
    value = nx * u1[0] + ny * u2[0];
  }
  else // Neumann
  {
    double p[4];
    ExactP(x, y, p);
    value = p[0];
  }
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y,
               const double *const*, double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->SIGMA_PERM;
  double u1[4], u2[4], p[4];
  for(int i = 0; i < n_points; i++)
  {
    ExactU1(X[i], Y[i], u1);
    ExactU2(X[i], Y[i], u2);
    ExactP(X[i], Y[i], p);
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = u1[0] + p[1]; // f1
    coeffs[i][2] = u2[0] + p[2]; // f2
    coeffs[i][3] = u1[1] + u2[2];// g
  }
}

