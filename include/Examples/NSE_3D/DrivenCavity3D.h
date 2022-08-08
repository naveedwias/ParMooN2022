/**
 * @file The driven cavity benchmark example in 3D.
 *
 * The boundary data is adapted to the [0.1]^3 unit cube example, which is
 * availabe as default geometry in ParMooN. It will throw an error if you
 * try running it on any other domain - just to make you aware of that fact.
 */

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("EXAMPLE"," DrivenCavity3D.h");
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_3d;
auto& ExactU2 = unknown_solution_3d;
auto& ExactU3 = unknown_solution_3d;
auto& ExactP = unknown_solution_3d;

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double, double, double z, double &value)
{
  double tol = 1e-10;

  if (std::abs(1-z)<tol)
  {
    value = 1.0;
  }
  else
    value = 0.0 ;
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
  static double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3
    coeff[4] = 0; // g
    coeff[5] = 0; // sigma
  }
}
