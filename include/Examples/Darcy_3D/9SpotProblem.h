// Bechmark for Darcy problem
// This is the extension of the 5SpotProblem to 3D, for the 2D version see
// Arif Masud, Thomas J.R. Hughes: "A stabilized mixed finite element method
// for Darcy flow", chapter 4.2, and also the file 5SpotProblem.h.
//
// The analytical solution in 2D (and supposely also here in 3D) exists, but is
// not easy to write. Therefore we set it to zero (unknown) here.

double mesh_width_at_corners = 1./16.;
double factor = 1.;

void ExampleFile()
{
  Output::print<1>("Example: 9SpotProblem.h with mesh width ",
                   mesh_width_at_corners, ". The hydraulic conductivity ", 
                   (factor == 1. ? "is uniform" : " has a checkerboard pattern"),
                   "(", (factor == 1. ? 1.0 : factor), ").");
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
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void FluxBoundValue(int, double x, double y, double z, double &value)
{
  double h = mesh_width_at_corners;
  value = 0;
  if(std::max({x, y, z}) <= h)
  {
    value = -1./(3*8*h);
  }
  else if(std::min({x, y, z}) >= 1-h)
  {
    value = 1./(3*8*h);
  }
}

// coefficients in the pde
void LinCoeffs(int n_points, const double* X, const double* Y, const double* Z,
               const double *const*, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    if(std::max({X[i], Y[i], Z[i]}) < 0.5 || std::min({X[i], Y[i], Z[i]}) > 0.5)
      coeffs[i][0] = 1.0;
    else 
      coeffs[i][0] = 1.0 * factor;
    coeffs[i][1] = 0;   // f1, forces
    coeffs[i][2] = 0;   // f2, forces
    coeffs[i][3] = 0;   // f3, forces
    coeffs[i][4] = 0;   // source terms
  }
}


