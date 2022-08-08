// Bechmark for Darcy problem, exact solution is
// 
// u(x,y) = (y,z,x)
// p(x,y) = 0

constexpr double tol = 1e-10;

void ExampleFile()
{
  Output::print<1>("Example: Benchmark.h");
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double, double z, double *values)
{
  values[0] = x-1;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double, double y, double, double *values)
{
  values[0] = y;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double, double, double z, double *values)
{
  values[0] = z;
  values[1] = 0;
  values[2] = 0;
  values[3] = 1;
  values[4] = 0;
}
void ExactP(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int, double x, double, double, BoundCond &cond)  
{
  if(std::abs(x) <= tol) // (x == 0.0)
    cond = NEUMANN;
  else
    cond = DIRICHLET;
}

// u \cdot n
void UNBoundValue(int, double x, double y, double z, double &value)
{
  // check if (x,y,z) is on the boundary, determine normal vector
  double nx = 0.;
  double ny = 0.;
  double nz = 0.;
  if(std::abs(x) <= tol) // (x == 0.0)
  {
    nx = -1.;
  }
  else if(std::abs(x-1.0) <= tol) //(x == 1.0)
  {
    nx = 1.;
  }
  else if(std::abs(y) <= tol) // (y == 0.0)
  {
    ny = -1.;
  }
  else if(std::abs(y-1.0) <= tol) // (y == 1.0)
  {
    ny = 1.;
  }
  else if(std::abs(z) <= tol) // (z == 0.0)
  {
    nz = -1.;
  }
  else if(std::abs(z-1.0) <= tol) // (z == 1.0)
  {
    nz = 1.;
  }
  else
  {
    ErrThrow("Evaluating boundary condition at (", x, ",", y, ",", z, ") which "
             "is not on the boundary!");
  }
  BoundCond cond;
  BoundCondition(-1, x, y, z, cond);
  if(cond == DIRICHLET)
  {
    double u1[5], u2[5], u3[5];
    ExactU1(x, y, z, u1);
    ExactU2(x, y, z, u2);
    ExactU3(x, y, z, u3);
    value = nx*u1[0] + ny*u2[0] + nz*u3[0];
  }
  else
  {
    double p[5];
    ExactP(x, y, z, p);
    value = p[0];
  }
  //Output::print("Boundary Data: ", value, " at (", x, ",", y, ",", z, ")");
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y, const double *Z,
               const double *const*, double **coeffs)
{
  const double eps = 1./TDatabase::ParamDB->SIGMA_PERM;
  double u1[5], u2[5], u3[5], p[5];
  for(int i = 0; i < n_points; i++)
  {
    ExactU1(X[i], Y[i], Z[i], u1);
    ExactU2(X[i], Y[i], Z[i], u2);
    ExactU3(X[i], Y[i], Z[i], u3);
    ExactP(X[i], Y[i], Z[i], p);
    
    coeffs[i][0] = eps;
    // RHS for exact solution
    coeffs[i][1] = eps*u1[0] + p[1]; // f1
    coeffs[i][2] = eps*u2[0] + p[2]; // f2
    coeffs[i][3] = eps*u3[0] + p[3]; // f2
    coeffs[i][4] = u1[1] + u2[2] + u3[3];// g
  }
}


