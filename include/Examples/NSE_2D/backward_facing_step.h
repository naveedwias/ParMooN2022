// Navier-Stokes problem, backward facing step
//
//  __________________________     
//  |                         |     y=1
//  |____                     |     y=0
//       |____________________|     y=-0.2
//
//x=0   x=1                   x=5

//WIll be manipulated by the example class.
double DIMENSIONLESS_VISCOSITY = 10;

void ExampleFile()
{
  Output::print<1>("Example: backward_facing_step.h");
}

// ========================================================================
// exact solution
// ========================================================================
auto& ExactU1 = unknown_solution_2d;
auto& ExactU2 = unknown_solution_2d;
auto& ExactP = unknown_solution_2d;

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double, BoundCond &cond)
{
  cond = i == 3 ? NEUMANN : DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  value = BdComp == 5 ? 4.*Param*(1.-Param) : 0.;
}

void U2BoundValue(int, double, double &value)
{
  value = 0.;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *, const double *,
               const double *const*, double **coeffs)
{
  for(auto i = 0; i < n_points; i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = 0.; // f1
    coeffs[i][2] = 0.; // f2
    coeffs[i][3] = 0.; // div(u)
    // additional coefficient (used only in the Brinkman problem)
    coeffs[i][4] = 0.;

  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void NonLinCoeffs(int n_points, const double *, const double *,
                  double *const*parameters, double **coeffs)
{
  for(auto i = 0; i < n_points; i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = parameters[i][0]; // u1
    coeffs[i][2] = parameters[i][1]; // u2

  }
}
