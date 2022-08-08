/**
 * Simple example for Convection diffusion 3D, with analytic solution,
 * used for development, debugging and testing.
 *
 * Constant concentration c = 1.
 * Diffusion coefficient 1, convection coefficients 1, reaction coefficient 1.
 * Dirichlet boundary conditions 1.
 *
 * @author ???, Clemens Bartsch imported this from MooNMD.
 * @date 2016/03/30 Import to ParMooN.
 *
 */

void ExampleFile()
{
  Output::root_info<1>("Example", "test_p0_zero.h");
}

// exact solution
void Exact(double, double, double, double *values)
{
  values[0] = 1;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double, double, double, double &value)
{
  value = 1;
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *, const double *const*, double **coeffs)
{
  double eps=1;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps; //diffusion coefficient
    coeff[1] = 1;   //ux
    coeff[2] = 1;   //uy
    coeff[3] = 1;   //uz
    coeff[4] = 1;   //reaction coefficient
    coeff[5] = 1;   //rhs
  }
}

