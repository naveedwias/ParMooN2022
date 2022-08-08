// Bechmark for Darcy problem, exact solution is in a BDM
// (Brezzi-Douglas-Marini) or RT (Raviart-Thomas) space.
// This is essentially taken from a Navier-Stokes example, but the right-hand
// side is computed differently.

#include "NSE_2D/Polynomial.h"

void ExampleFileDarcy()
{
  ExampleFile();
}

// ========================================================================
// coefficients for Darcy
// ========================================================================
void LinCoeffsDarcy(int n_points, const double *X, const double *Y,
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

