// ==========================================================================
// instationary problem
// ==========================================================================
void ExampleFile()
{
  Output::print<1>("Example: quadratic_space_time.h\n");
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*(x*x+y*z);
  values[1] = 2*x*t*t;
  values[2] = z*t*t;
  values[3] = y*t*t;
  values[4] = 2*t*t;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*t*(x*x+y*z);
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*t*(x*x+y*z);
}

void BilinearCoeffs(int n_points, const double *X, const double *Y,
                    const double *Z, const double *const*, double **coeffs)
{
  double eps = 1;
  double *coeff;
  double t=TDatabase::TimeDB->CURRENTTIME;
  
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
    coeff[5] = 2*t*(x*x+y*z) -eps * 2*t*t + 4 * 2*x * t*t+ 3 * z *t*t+ 2 * y*t*t + 1 * t*t*(x*x+y*z); //rhs
  }
}

