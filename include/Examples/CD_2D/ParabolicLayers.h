// ======================================================================
// smooth solution problem
// ======================================================================
#define __PARABOLIC_LAYERS__
#define __STORE_ANALYTIC__

double DIFFUSION;

void ExampleFile()
{
  Output::print<1>("Example: PrabolicLayers.h");
}

// exact solution
void Exact(double x, double y, double *values)
{
    values[0] = x;
    values[1] = 1;
    values[2] = 0;
    values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *const*, double **coeffs)
{
  double eps=DIFFUSION;
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 1; //x-coordiante of convection b
    coeff[2] = 0; //y-coordiante of convection b
    coeff[3] = 0; //reaction coeffecient c
    coeff[4] = 1;
  }
}
