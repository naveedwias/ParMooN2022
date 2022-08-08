/**
 * Simple example for Convection diffusion 3D, with analytic solution,
 * used for development, debugging and testing.
 *
 * Linear concentration solution, c = 1+2*x+3*y+4*z.
 * Diffusion coefficient 1, convection (1,2,3), no reaction.
 * Dirichlet boundary conditions.
 *
 * @author ???, Clemens Bartsch imported this from MooNMD.
 * @date 2016/03/30 Import to ParMooN.
 *
 */

void ExampleFile()
{
  Output::root_info<1>("Example", "test_p1.h");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 1+2*x+3*y+4*z;
  values[1] = 2;
  values[2] = 3;
  values[3] = 4;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  value = 1+2*x+3*y+4*z;
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *, const double *const*, double **coeffs)
{
  double eps = 1; //diffusion coefficient
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1; // 1;
    coeff[2] = 2; // 2;
    coeff[3] = 3; // 3;
    coeff[4] = 0; // 0;
    coeff[5] = 1*2+2*3+3*4; // 1*2+2*3+3*4;
  }
}

