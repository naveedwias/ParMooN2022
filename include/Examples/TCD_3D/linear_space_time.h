void ExampleFile()
{
  Output::print<1>("Example: linear_space_time.h\n") ;
}

void Exact(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = 1+2*x+3*t*y+4*t*z;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 4*t;
  values[4] = 0;
}

void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
    values[0] = 1+2*x+3*t*y+4*t*z;
}

void BoundValue(int, double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = 1+2*x+3*t*y+4*t*z;
}

void BilinearCoeffs(int n_points, const double *, const double *y,
                    const double *z, const double *const*, double **coeffs)
{
  double eps = 1; 
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff;

  for(int i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1;
    coeff[2] = 2;
    coeff[3] = 3;
    coeff[4] = 0;
    
    coeff[5] = 3*y[i] + 4*z[i] + 2 + 6*t + 12*t; 
  }
}

