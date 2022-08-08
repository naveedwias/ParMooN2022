// time dependent Navier-Stokes problem 3D, ansatz
// 
// u(t,x,y,z) = t*(y+z, x-z, 2*x+y)
// p(t,x,y,z) = 0

// This is also called nu, or eps, it is equal
// to 1/Reynolds_number and is dimensionless
double DIMENSIONLESS_VISCOSITY;

void ExampleFile()
{
  Output::root_info<1>("Example", "AnsatzLinConst.h");
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(y+z);
}

void InitialU2(double x, double, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(x-z);
}

void InitialU3(double x, double y, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(2*x+y);

}

void InitialP(double, double, double, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(y+z);
  values[1] = 0;
  values[2] = t;
  values[3] = t;
  values[4] = 0;
}

void ExactU2(double x, double, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(x-z);
  values[1] = t;
  values[2] = 0;
  values[3] = -t;
  values[4] = 0;
}

void ExactU3(double x, double y, double, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*(2*x+y);
  values[1] = 2*t;
  values[2] = t;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double, double, double, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;

}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(int, double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(int, double, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(y+z);
}

// value of boundary condition
void U2BoundValue(int, double x, double, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(x-z);
}

// value of boundary condition
void U3BoundValue(int, double x, double y, double, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = t*(2*x+y);    
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, const double *X, const double *Y, const double *Z,
               const double *const*, double **coeffs)
{
  static double eps = DIMENSIONLESS_VISCOSITY;
  int i;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
    coeff[1] = (y+z) + t*t*(x-z) + t*t*(2*x+y); // f1
    coeff[2] = (x-z) + t*t*(y+z) - t*t*(2*x+y); // f2
    coeff[3] =  (2*x+y) + 2*t*t*(y+z) + t*t*(x-z);// f3
  }  
}


