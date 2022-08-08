// ======================================================================
// Sine problem 3D
// ======================================================================
// #include <ConvDiff3D.h>

void ExampleFile()
{
  Output::root_info<1>("Example", "Aorta_Laplace.h");
}

double PECLET_NUMBER;
// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[1] = M_PI*std::cos(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
  values[2] = M_PI*std::sin(M_PI*x)*std::cos(M_PI*y)*std::sin(M_PI*z);
  values[3] = M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::cos(M_PI*z);
  values[4] = -3*M_PI*M_PI*std::sin(M_PI*x)*std::sin(M_PI*y)*std::sin(M_PI*z);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int reference, double, double, double, BoundCond &cond)
{
  if (reference == 1) {
    cond = NEUMANN;
  } else
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int reference, double, double, double, double &value)
{
  switch (reference) {
  case 1:
    value = 0;
    break;
  case 2:
    value = 0;
    break;
  case 3:
    value = 10;
    break;
  case 4:
    value = 20;
    break;
  case 5:
    value = 30;
    break;
  case 6:
    value = 40;
    break;
  case 7:
    value = 100;
    break;
  default:
    Output::print(" ** ERROR: no condition for bd. reference = ",reference);
  }
    
}

void BilinearCoeffs(int n_points, const double *, const double *,
                    const double *, const double *const*, double **coeffs)
{
  const double eps = PECLET_NUMBER;
  for(int i = 0; i < n_points; i++)
  {

    coeffs[i][0] = eps;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
    coeffs[i][3] = 0;
    coeffs[i][4] = 0;
    coeffs[i][5] = 0;
  }
}

