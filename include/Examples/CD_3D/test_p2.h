/**
 * Simple example for Convection diffusion 3D, with analytic solution,
 * used for development, debugging and testing.
 *
 * Quadratic concentration solution, c = x^2+y*z.
 * Diffusion coefficient 1, convection (4,3,2), reaction coefficient 1.
 * Dirichlet boundary conditions.
 *
 * @author ???, Clemens Bartsch imported this from MooNMD.
 * @date 2016/03/30 Import to ParMooN.
 *
 */

void ExampleFile()
{
  Output::root_info<1>("Example", "test_p2.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = x*x+y*z;
  values[1] = 2*x;
  values[2] = z;
  values[3] = y;
  values[4] = 2;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  value = x*x+y*z;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *Z, const double *const*, double **coeffs)
{
  double eps = 1;
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    double x = X[i];
    double y = Y[i];
    double z = Z[i];

    coeff[0] = eps; //diffusion coefficient
    coeff[1] = 4;   //ux
    coeff[2] = 3;   //uy
    coeff[3] = 2;   //uy
    coeff[4] = 1;   //reaction coefficient
    coeff[5] = -eps * 2 + 4 * 2*x + 3 * z + 2 * y + 1 * (x*x+y*z); //rhs
  }
}

