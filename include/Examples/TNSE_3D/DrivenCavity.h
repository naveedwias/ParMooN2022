// Navier-Stokes problem, Driven Cavity
// 
// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;
double pi = 3.14159265358979;

void ExampleFile()
{
  Output::print("Example: DrivenCavity.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialU2(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialU3(double, double, double, double *values)
{
  values[0] = 0;
}

void InitialP(double, double, double, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP = unknown_solution_3d;

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double x, double y, double z, double &value)
{
  double eps = 1e-8;

  if (std::abs(z-1)<eps)
  {
    if ((std::abs(x)>eps)&&(std::abs(1-x)>eps)&&(std::abs(y)>eps)&&(std::abs(1-y)>eps))
      value = 1.0;
    else
      value = 0.0;
  }
  else
    value =0.0 ;
}

// value of boundary condition
void U2BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int, double, double, double, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *, const double *,
               const double *const*, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = 0; // f1
    coeffs[i][2] = 0; // f2
    coeffs[i][3] = 0; // f3
    coeffs[i][4] = 0; // f3
  }
}


